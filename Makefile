MAKEFLAGS += --warn-undefined-variables
.SHELLFLAGS := -eu -o pipefail -c
.DEFAULT_GOAL := all

# Root directory - directory of the Makefile
root_dir := $(dir $(realpath $(firstword $(MAKEFILE_LIST))))
# Directory of the Methpype package
package_dir := $(cur_dir)methpype
# Directory for the docs
docs_dir := $(cur_dir)docs

# Path to the conda environment file
conda_env_path := $(root_dir)conda-env.yml
# Path to the requirements file
requirements_path := $(root_dir)requirements.txt
# Local vars
venv_name := methpype

.PHONY: all
all: init activate #help: inits and activates the environment

.PHONY: init
init: #help: Creates the virtual environment
	@echo 'Creating conda environment:'
	@echo '$(venv_name)'
	@echo ''
	@echo 'Command:'
	conda env create --file $(conda_env_path)

.PHONY: update
update: #help: Updates the virtual environment
	@echo 'Updating conda environment:'
	@echo '$(venv_name)'
	@echo ''
	@echo 'Command:'
	conda env update --file $(conda_env_path)

.PHONY: test
.SILENT: test
test: #help: Runs the unit test suite
	( \
	   source activate $(venv_name); \
	   python setup.py test; \
	)

clean: #help: Cleans python cache files
	@echo 'Cleaning python cache files'
	( \
		find . -name '*.pyc' -type f -exec rm -rf {} + ; \
		find . -name '__pycache__' -type d -exec rm -rf {} + ; \
	)
