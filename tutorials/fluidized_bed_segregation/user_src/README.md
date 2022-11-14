# Fluidized Bed Segregation Tutorial

## CFD_user.c, CFD_user_ext.c

- Converts frames of CFD-DEM based field information into Fluent *.ip format
- The base version just considers the available restricted set of CFD-DEM data

## rCFD_user, rCFD_user_ext, rCFD_user_batch

- Controls rCFD_prep.c and rCFD_run.c
- The base version just considers the available restricted set of CFD-DEM data
- The _batch version is used for batch run
- The _ext version considers reference data which are available upon request

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

- This software is distributed under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.html).
- Copyright © 2021 JKU Linz
