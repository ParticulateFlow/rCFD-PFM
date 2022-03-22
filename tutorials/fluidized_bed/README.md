# Fluidized Bed Tutorial

## Contents

- Prepares rCFD database on a 65k grid (C2C patterns for gas and solid phase).
- Runs rCFD considering (i) secondary gas species, (ii) drifting solid fractions

## Prepare rCFD

- start Ansys/fluent by e.g. *fluent 3ddp -t4 &* from this folder
- read FB_tutorial.scm
- type into fluent console (rcfd_p1), followed by (rcfd_p1), ..., (rcfd_p11); or simply (rcfd_prep)

- c2c patterns will be stored into ./data/c2c
- vof fields will be stored into ./data/vof
- jump vector will be stored into ./rec

- close Ansys/fluent

## Run rCFD

- start Ansys/fluent with the same number of threads from this folder.
- read FB_tutorial.scm
- type into fluent console (rcfd_r1), followed by (rcfd_r2), ..., (rcfd_r7); or simply (rcfd_run)

- balances will be stored into ./post
- run transcript will be stored into ./post
- animation snap-shots will be stored into ./post

- close Ansys/fluent

## Reference CFD simulation

- start Ansys/fluent by e.g. *fluent 3ddp -t4 &* from this folder
- read FB_tutorial.scm
- type into fluent console (cfd_r1), followed by (cfd_r2), ..., (cfd_r6); or simply (cfd_ref)

- monitors of *secondary gas mass* and *solid mixing index* will be stored into ./ref
- animation snap-shots of secondary gas and solid species concentration will be stored into ./ref

- close Ansys/fluent

## Known Limitations

- (rcfd_run) works only after (rcfd_prep)
- missing cross-partition diffusion

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

- This software is distributed under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.html).
- Copyright © 2021 JKU Linz
