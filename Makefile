
.PHONY: environment remove-env install uninstall test # .PHONY is something we can add when our target dependencies are not files.

MODULE=mcsim
ENVIRONMENT=chem274A_lab1

environment: remove-env
	conda create -n $(ENVIRONMENT) "python=3.10" --yes
	conda install -c conda-forge notebook numpy matplotlib pytest-cov --name $(ENVIRONMENT) --yes

remove-env:
	conda remove --name $(ENVIRONMENT) --all --yes

install: uninstall ## install the package to the active Python's site-packages
	pip install .

uninstall: ## uninstall the package
	pip uninstall --yes $(MODULE)

test: 
	pytest -v
