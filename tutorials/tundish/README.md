# Tundish Tutorial

## Contents

- Prepares rCFD database on a 1.2M grid.
- Runs rCFD considering (i) incoming drifting particles, (ii) change in incoming steel grade 
and (iii) change in incoming temperature

## Prepare rCFD

- start Ansys/fluent by e.g. *fluent 3ddp -t12 &* from this folder
- read tundish_tutorial.scm
- type into fluent console (rcfd_p1), followed by (rcfd_p1), ..., (rcfd_p11); or simply (rcfd_prep)

- c2c patterns will be stored into ./data/c2c
- jump vector will be stored into ./rec

- close Ansys/fluent

## Run rCFD

- start Ansys/fluent with the same number of threads from this folder.
- read tundish_tutorial.scm
- type into fluent console (rcfd_r1), followed by (rcfd_r1), ..., (rcfd_r7); or simply (rcfd_run)

- balances will be stored into ./post
- run transcript will be stored into ./post
- animation snap-shots will be stored into ./post

- close Ansys/fluent

## Known Limitations

- (rcfd_run) works only after (rcfd_prep)
- deviating global balances
- missing cross-partition diffusion

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

- This software is distributed under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.html).
- Copyright © 2021 JKU Linz
