"""Handle user arguments when running AstroCats
"""
import argparse


class ArgsHandler:

    def __init__(self, log):
        self.log = log
        parser = self._setup_argparse()
        self.parser = parser
        return

    def run_subcommand(self, args, catalog):
        # Data Import
        # -----------
        if args.subcommand == 'import':
            self.log.info("Running 'import'.")
            catalog.import_data()

        # Producer
        # --------
        elif args.subcommand == 'produce':
            self.log.info("Running 'produce'.")
            from astrocats.producer import webcat
            webcat.main(catalog)

        # Git Subcommands
        # ---------------
        elif args.subcommand == 'git-clone':
            self.log.info("Running 'git clone'.")
            catalog.git_clone_all_repos()
        elif args.subcommand == 'git-push':
            self.log.info("Running 'git push'.")
            catalog.git_add_commit_push_all_repos()
        elif args.subcommand == 'git-pull':
            self.log.info("Running 'git pull'.")
            # catalog.git_pull_all_repos()
            raise RuntimeError("NOT YET IMPLEMENTED!")
        elif args.subcommand == 'git-reset-local':
            self.log.info("Running 'git reset' using the local HEAD.")
            catalog.git_reset_all_repos(hard=True, origin=False, clean=True)
        elif args.subcommand == 'git-reset-origin':
            self.log.info("Running 'git reset' using 'origin/master'.")
            catalog.git_reset_all_repos(hard=True, origin=True, clean=True)
        elif args.subcommand == 'git-status':
            self.log.info("Running 'git status'.")
            catalog.git_status_all_repos()

        # Analyze Catalogs
        # ----------------
        elif args.subcommand == 'analyze':
            self.log.info("Running 'analyze'.")
            from .analysis import Analysis
            # Create an `Analysis` instance
            lysis = Analysis(catalog, self.log)
            # Pass the command-line arguments to run.
            lysis.analyze(args)

        return

    def load_args(self, args, clargs):
        """Parse arguments and return configuration settings.
        """
        # Parse All Arguments
        args = self.parser.parse_args(args=clargs, namespace=args)

        # Print the help information if no subcommand is given
        # subcommand is required for operation
        if args.subcommand is None:
            self.parser.print_help()
            args = None

        return args

    def _setup_argparse(self):
        """Create `argparse` instance, and setup with appropriate parameters.
        """
        parser = argparse.ArgumentParser(
            prog='catalog', description='Parent Catalog class for astrocats.')

        subparsers = parser.add_subparsers(
            description='valid subcommands', dest='subcommand')

        # Data Import
        # -----------
        # Add the 'import' command, and related arguments
        self._add_parser_arguments_import(subparsers)

        # Producer
        # --------
        # Add the 'import' command, and related arguments
        self._add_parser_arguments_produce(subparsers)

        # Git Subcommands
        # ---------------
        self._add_parser_arguments_git(subparsers)

        # Analyze Catalogs
        # ----------------
        # Add the 'analyze' command, and related arguments
        self._add_parser_arguments_analyze(subparsers)

        return parser

    def _add_parser_arguments_import(self, subparsers):
        """Create parser for 'import' subcommand, and associated arguments.
        """
        import_pars = subparsers.add_parser(
            "import", help="Import data.")

        import_pars.add_argument(
            '--update', '-u', dest='update',
            default=False, action='store_true',
            help='Only update catalog using live sources.')
        import_pars.add_argument(
            '--archived', '-a', dest='archived',
            default=False, action='store_true',
            help='Always use task caches.')

        # Control which 'tasks' are executed
        # ----------------------------------
        import_pars.add_argument(
            '--tasks', dest='args_task_list', nargs='*', default=None,
            help='space delimited list of tasks to perform.')
        import_pars.add_argument(
            '--yes', dest='yes_task_list', nargs='+', default=None,
            help='space delimited list of tasks to turn on.')
        import_pars.add_argument(
            '--no', dest='no_task_list', nargs='+', default=None,
            help='space delimited list of tasks to turn off.')
        import_pars.add_argument(
            '--min-task-priority', dest='min_task_priority',
            default=None,
            help='minimum priority for a task to run')
        import_pars.add_argument(
            '--max-task-priority', dest='max_task_priority',
            default=None,
            help='maximum priority for a task to run')
        import_pars.add_argument(
            '--task-groups', dest='task_groups',
            default=None,
            help='predefined group(s) of tasks to run.')

        return import_pars

    def _add_parser_arguments_produce(self, subparsers):
        """Create parser for 'produce' subcommand, and associated arguments.
        """
        produce_pars = subparsers.add_parser(
            "produce", help="Generate a catalog JSON file and plot HTML files from SNE data.")

        produce_pars.add_argument(
            '--no-write-catalog', '-nwc',
            dest='writecatalog', default=True, action='store_false',
            help="Don't write catalog file")
        produce_pars.add_argument(
            '--no-write-html', '-nwh',
            dest='writehtml', default=True, action='store_false',
            help='Don\'t write html plot files')
        produce_pars.add_argument(
            '--no-collect-hosts', '-nch',
            dest='collecthosts', default=True, action='store_false',
            help='Don\'t collect host galaxy images')
        produce_pars.add_argument(
            '--force-html', '-fh',
            dest='forcehtml', default=False,
            help='Force write html plot files',
            action='store_true')
        produce_pars.add_argument(
            '--event-list', '-el',
            dest='eventlist', default=[], type=str, nargs='+',
            help='Process a list of events')
        produce_pars.add_argument(
            '--test', '-te',
            dest='test', default=False, action='store_true',
            help='Test this script')
        produce_pars.add_argument(
            '--travis', '-tr',
            dest='travis', default=False, action='store_true',
            help='Set some options when using Travis')
        produce_pars.add_argument(
            '--boneyard', '-by',
            dest='boneyard', default=False, action='store_true',
            help='Make "boneyard" catalog')
        produce_pars.add_argument(
            '--delete-orphans', '-do',
            dest='deleteorphans', default=False, action='store_true',
            help='Delete orphan JSON files')

        return produce_pars

    def _add_parser_arguments_git(self, subparsers):
        """Create a sub-parsers for git subcommands.
        """
        subparsers.add_parser(
            "git-clone",
            help="Clone all defined data repositories if they dont exist.")

        subparsers.add_parser(
            "git-push",
            help="Add all files to data repositories, commit, and push.")

        subparsers.add_parser(
            "git-pull",
            help="'Pull' all data repositories.")

        subparsers.add_parser(
            "git-reset-local",
            help="Hard reset all data repositories using local 'HEAD'.")

        subparsers.add_parser(
            "git-reset-origin",
            help="Hard reset all data repositories using 'origin/master'.")

        subparsers.add_parser(
            "git-status",
            help="Get the 'git status' of all data repositories.")

        return

    def _add_parser_arguments_analyze(self, subparsers):
        """Create a parser for the 'analyze' subcommand.
        """
        lyze_pars = subparsers.add_parser(
            "analyze",
            help="Perform basic analysis on this catalog.")

        lyze_pars.add_argument(
            '--count', '-c', dest='count',
            default=False, action='store_true',
            help='Determine counts of entries, files, etc.')

        return lyze_pars
