"""
"""
from glob import glob
import os
import sys

from astrocats import main, _ROOT_BASE_PATH
from astrocats.catalog.utils import read_json_dict, repo_priority

_STR_DEPRECATION_ERROR = "WARNING: using deprecated `{}`, use `{}` instead."


class Paths:
    """Store and control catalog file-structure information.

    Individual catalogs must provide the below file structure.
    -   `repos.json`
    -   `tasks.json`

    Attributes
    ----------
    catalog_dir : str
    tasks_dir : str
    catalog_dir : str
    input : str
    output : str
    repos_list_fname : str
    task_list_fname : str
    repos_dict : dict
        Dictionary of 'repo-types: repo-lists' key-value pairs.
        Loaded from `repos_list_fname` file.

    Methods
    -------
    get_all_repo_folders : get a list of paths for all data repositories
    get_repo_boneyard : get the path of the boneyard repository
    get_repo_input_folders : get the paths of all input data repositories
    get_repo_output_file_list : get the paths of all files in output repos
    get_repo_output_folders : get the paths of all input data repositories

    """

    def __init__(self, catalog=None, log=None):
        # New parameter for the overall root directory for astrocats as a whole
        self.root = str(_ROOT_BASE_PATH)

        this_file = None
        if catalog is not None:
            try:
                this_file = sys.modules[catalog.__module__].__file__
            except AttributeError:
                try_dirs = [os.path.join(self.root, 'astrocats'), self.root, '']
                for td in try_dirs:
                    this_file = os.path.join(td, catalog, '')
                    if os.path.isdir(this_file):
                        break

                    import warnings
                    warnings.warn("Couldnt find '{}' in '{}'...".format(catalog, td))
                    this_file = None
            else:
                log = catalog.log

        else:
            this_file = sys.modules[self.__module__].__file__

        if this_file is None:
            raise RuntimeError("`this_file` is None!  catalog = '{}'.".format(catalog))

        # If no log is given (and no catalog), construct a new one
        if log is None:
            log = main.load_log()

        self.log = log
        # Catalog Directories
        self.catalog_dir = os.path.dirname(this_file)
        self.catalog_parent = os.path.abspath(os.path.join(self.catalog_dir, os.pardir))
        self.tasks_dir = os.path.join(self.catalog_dir, 'tasks')
        self.input = os.path.join(self.catalog_dir, 'input', '')
        self.output = os.path.join(self.catalog_dir, 'output', '')
        self.cache = os.path.join(self.output, 'cache')

        # Catalog data files
        self.repos_list_fname = os.path.join(self.input, 'repos.json')
        self.task_list_fname = os.path.join(self.input, 'tasks.json')
        self.repos_dict = read_json_dict(self.repos_list_fname)
        self.md5_file = os.path.join(self.output, 'md5s.json')
        self.host_imgs_file = os.path.join(self.output, 'hostimgs.json')
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
            self.log.debug("Found {} files in '{}'".format(
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
        bone_path = os.path.join(self.output, bone_path, '')
        return bone_path

    def get_repo_input_folders(self):
        """Get the full paths of the input data repositories.
        """
        repo_folders = []
        repo_folders += self.repos_dict['external']
        repo_folders += self.repos_dict['internal']
        repo_folders = list(sorted(set(repo_folders)))
        repo_folders = [os.path.join(self.input, rf)
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
        repo_folders = [os.path.join(self.output, rf)
                        for rf in repo_folders if len(rf)]
        return repo_folders

    @property
    def PATH_INPUT(self):
        self.log.error(_STR_DEPRECATION_ERROR.format("PATH_INPUT", "input"))
        return self.input

    @property
    def PATH_OUTPUT(self):
        self.log.error(_STR_DEPRECATION_ERROR.format("PATH_OUTPUT", "output"))
        return self.output
