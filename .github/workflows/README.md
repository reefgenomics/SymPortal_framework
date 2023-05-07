# Symportal CI Workflows

We have two main workflows for Continuous Integration (CI) with GItHub Actions:

* [Create and Cache Conda Environment](https://github.com/greenjune-ship-it/SymPortal/actions/workflows/setup-conda-env.yml)
* [Run Tests](https://github.com/greenjune-ship-it/SymPortal/actions/workflows/run-tests.yml)

# Create and Cache Conda Environment

This is a workflow for creating and caching a Conda environment named `symportal`. This caching step can help speed up future builds by reusing the cached environment if the `environment.yml` file has not changed.

The workflow is triggered when there are changes to either the `environment.yml` file or any file that matches the pattern `**/setup-conda-env.yml` in the repository.

The workflow defines an environment variable named `ENV_NAME`, which is the name of the Conda environment to be created. Another environment variable named `CACHE_PATH`, which is the location where the Conda environment will be cached.

The job consists of three steps:

1) The first step checks out the repository code.
2) The second step sets up the Conda environment and installs dependencies listed in the `environment.yml` file. 
3) The third step caches the Conda environment at the specified path with a key that is unique to the environment name and the hash of the `environment.yml` file.

# Run Tests

The workflow helps automate the testing process for a Python project using GitHub Actions.

This workflow that is triggered when a push event occurs to the repository's `environment.yml` or any file that matches the pattern `**/run-tests.yml` file paths.

The workflow defines an environment variables called `ENV_NAME` and `CACHE_PATH` that are used by the job.

The job has several steps that are executed sequentially:

1) The first step checks out the repository code.
2) The "Set up Conda environment" step sets up the Conda environment with the name specified in the `ENV_NAME` environment variable.
3) The "Restore Conda environment from cache" step restores the Conda environment from the cache specified in the `CACHE_PATH` environment variable.
4) The "Set up settings.py" step sets up the settings.py file by renaming the settings_blank.py file to settings.py and replacing the `SECRET_KEY` with a randomly generated value.
5) The "Configure sp_config.py" step sets up the sp_config.py file by renaming the sp_config_blank.py file to sp_config.py.
6) The "Create SQLite database" step creates an SQLite database for Django application by running the `python manage.py migrate` command.
7) The "Populate local database" step populates the local database by running the `python populate_db_ref_seqs.py` command.
8) The "Run tests" step runs the tests by running the `python -m tests.tests` command.

For more details regarding [SymPortal](https://github.com/didillysquat/SymPortal_framework/wiki/0_3_18_19_SymPortal-setup) setup section of project [Wiki](https://github.com/didillysquat/SymPortal_framework/wiki).

# Additional Information

Workflow is set up by Yulia Iakovleva [yulia_iakovleva@proton.me]().
