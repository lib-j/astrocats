"""
"""
from glob import glob
import os
import sys

from astrocats.catalog.utils import (compress_gz, is_integer, pbar,
                                     read_json_dict, repo_priority,
                                     uncompress_gz, uniq_cdl)


class Paths:
    """Store and control catalog file-structure information.

    Individual catalogs must provide the below file structure.
    -   `repos.json`
    -   `tasks.json`

    Attributes
    ----------
    catalog : `astrocats.catalog.catalog.Catalog` (sub)class object
    catalog_dir : str
    tasks_dir : str
    PATH_BASE : str
    PATH_INPUT : str
    PATH_OUTPUT : str
    REPOS_LIST : str
    TASK_LIST : str
    repos_dict : dict
        Dictionary of 'repo-types: repo-lists' key-value pairs.
        Loaded from `REPOS_LIST` file.

    Methods
    -------
    get_all_repo_folders : get a list of paths for all data repositories
    get_repo_boneyard : get the path of the boneyard repository
    get_repo_input_folders : get the paths of all input data repositories
    get_repo_output_file_list : get the paths of all files in output repos
    get_repo_output_folders : get the paths of all input data repositories

    """

    def __init__(self, catalog):
        self.catalog = catalog
        this_file = sys.modules[self.__module__].__file__
        self.catalog_dir = os.path.dirname(this_file)
        self.tasks_dir = os.path.join(self.catalog_dir, 'tasks')
        self.PATH_BASE = os.path.join(
            catalog.args.base_path, self.catalog_dir, '')
        self.PATH_INPUT = os.path.join(self.PATH_BASE, 'input', '')
        self.PATH_OUTPUT = os.path.join(self.PATH_BASE, 'output', '')
        # critical datafiles
        self.REPOS_LIST = os.path.join(self.PATH_INPUT, 'repos.json')
        self.TASK_LIST = os.path.join(self.PATH_INPUT, 'tasks.json')
        self.repos_dict = read_json_dict(self.REPOS_LIST)
        return

    def _get_repo_file_list(self, repo_folders, normal=True, bones=True):
        """Get filenames for files in each repository, `boneyard` optional.
        """
        # repo_folders = get_repo_output_folders()
        files = []
        for rep in repo_folders:
            if 'boneyard' not in rep and not normal:
                continue
            if not bones and 'boneyard' in rep:
                continue
            these_files = glob(rep + "/*.json") + glob(rep + "/*.json.gz")
            self.catalog.log.debug("Found {} files in '{}'".format(
                len(these_files), rep))
            files += these_files

        return files

    def get_all_repo_folders(self, boneyard=True):
        """Get the full paths of all data repositories.
        """
        all_repos = self.get_repo_input_folders()
        all_repos.extend(self.get_repo_output_folders(bones=boneyard))
        return all_repos

    def get_repo_boneyard(self):
        bone_path = self.repos_dict['boneyard']
        try:
            bone_path = bone_path[0]
        except TypeError:
            pass
        bone_path = os.path.join(self.PATH_OUTPUT, bone_path, '')
        return bone_path

    def get_repo_input_folders(self):
        """Get the full paths of the input data repositories.
        """
        repo_folders = []
        repo_folders += self.repos_dict['external']
        repo_folders += self.repos_dict['internal']
        repo_folders = list(sorted(set(repo_folders)))
        repo_folders = [os.path.join(self.PATH_INPUT, rf)
                        for rf in repo_folders if len(rf)]
        return repo_folders

    def get_repo_output_file_list(self, normal=True, bones=True):
        """Get a list of all existing output files.

        These are the files deleted in the `delete_old_entry_files` task.
        """
        repo_folders = self.get_repo_output_folders()
        return self._get_repo_file_list(
            repo_folders, normal=normal, bones=bones)

    def get_repo_output_folders(self, bones=True):
        """Get the full paths of the output data repositories.
        """
        repo_folders = []
        repo_folders += self.repos_dict['output']
        if bones:
            repo_folders += self.repos_dict['boneyard']
        repo_folders = list(sorted(list(set(repo_folders)),
                                   key=lambda key: repo_priority(key)))
        repo_folders = [os.path.join(self.PATH_OUTPUT, rf)
                        for rf in repo_folders if len(rf)]
        return repo_folders
