Installation
============
This part shows how to install BioKowloon in a virtual environment.
***
It is recommended to install BioKowloon in a virtual environment, such as `conda`. This will isolate the package from your system Python installation and prevent any conflicts with other packages.

## Create a virtual environment
- Install `conda` if you don't already have it. For the usage of `conda`, please refer to [conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
- Create a new virtual environment named deside and activate it:

```shell
conda create -n biokowloon python=3.10
conda activate biokowloon
```

## Update pip 
Since we are using `pyproject.toml` to manage dependencies, the minimum required version of `pip` should be 10.0 or higher.
```shell
python3 -m pip install --upgrade pip
```

## Install DeSide

```shell
pip install biokowloon
```
