# pz-rail-lephare

[![Template](https://img.shields.io/badge/Template-LINCC%20Frameworks%20Python%20Project%20Template-brightgreen)](https://lincc-ppt.readthedocs.io/en/latest/)
[![codecov](https://codecov.io/gh/LSSTDESC/rail_lephare/graph/badge.svg?token=I597zffkaM)](https://codecov.io/gh/LSSTDESC/rail_lephare)
[![PyPI](https://img.shields.io/pypi/v/pz-rail-lephare?color=blue&logo=pypi&logoColor=white)](https://pypi.org/project/pz-rail-lephare/)

## RAIL LePHARE

This package allows users to run the [LePHARE code](https://github.com/lephare-photoz/lephare/) 
through [RAIL](https://github.com/LSSTDESC/RAIL). 
It can be installed via pip with the following command. 
Full documentation is available at the main [LePHARE documentation](https://lephare.readthedocs.io/en/latest/index.html).

```console
pip install pz-rail-lephare
```

## RAIL: Redshift Assessment Infrastructure Layers

This package is part of the larger ecosystem of Photometric Redshifts
in [RAIL](https://github.com/LSSTDESC/RAIL). This package in particular
is concerned with wrapping the [LePHARE](https://gitlab.lam.fr/Galaxies/LEPHARE/) template fitting code within RAIL.

### Developer installation

This package can be installed from source using the following commands to allow working on the code.
```console
git clone https://github.com/LSSTDESC/rail_lephare.git
cd rail_lephare
conda install -c conda-forge pytables
pip install -e '.[dev]'
```

### Citing RAIL

This code, while public on GitHub, has not yet been released by DESC and is
still under active development. Our release of v1.0 will be accompanied by a
journal paper describing the development and validation of RAIL.

If you make use of the ideas or software in RAIL, please cite the repository 
<https://github.com/LSSTDESC/RAIL>. You are welcome to re-use the code, which
is open source and available under terms consistent with the MIT license.

External contributors and DESC members wishing to use RAIL for non-DESC projects
should consult with the Photometric Redshifts (PZ) Working Group conveners,
ideally before the work has started, but definitely before any publication or 
posting of the work to the arXiv.

### Citing this package

If you use this package, you should also cite the appropriate papers for each
code used.  A list of such codes is included in the 
[Citing RAIL](https://lsstdescrail.readthedocs.io/en/stable/source/citing.html)
section of the main RAIL Read The Docs page.

