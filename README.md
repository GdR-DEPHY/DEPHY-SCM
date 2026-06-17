# DEPHY-SCM: Single-Column Model standards and case drivers

Detailed information about the DEPHY-SCM Standards can be found in DEPHY-SCM\_CommonStandards\_v1.0.pdf. 
Latest version of the standards is open for discussion under <a href="https://docs.google.com/document/d/1eAWY-ELL5Ua6a9WIsv4ODHmLXvfgla5TNQAuAwNASo0" target="_blank">here</a>.

The DEPHY-SCM tool is using Python 3

# pip installation

```shell
mkdir -p $HOME/.local/pyenvs/
python -m venv $HOME/.local/pyenvs/dephy-scm
source $HOME/.local/pyenvs/dephy-scm/bin/activate
cd dephy-scm # si on est pas déjà dedans...
pip install -e .
```

# Activate environment before running driver scripts!

```shell
source $HOME/.local/pyenvs/dephy-scm/bin/activate
```
# Required Python packages

They are listed in pyproject.toml and will be installed automatically with pip
install. 

  * netCDF4
  * numpy
  * scipy
  * xarray
  * matplotlib

# Using the tool

The tools is provided as a Python module named dephycf. 
To use it just activate your environment with

```shell
source $HOME/.local/pyenvs/dephy-scm/bin/activate
```
