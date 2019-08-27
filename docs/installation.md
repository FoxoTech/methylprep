# Installation

methylprep maintains configuration files for your Python package manager of choice: [conda](https://conda.io), [pipenv](https://pipenv.readthedocs.io/en/latest/), and [pip](https://pip.pypa.io/en/stable/).

---

## conda

```Shell
conda env create --file conda-env.yml
```

---

## pipenv

When using the library as-is

```Shell
pipenv install
```

When working on the library

```Shell
pipenv install -d
```

---

## pip

```Shell
pip install -r requirements.txt
```

---

## make

```Shell
make init
```
