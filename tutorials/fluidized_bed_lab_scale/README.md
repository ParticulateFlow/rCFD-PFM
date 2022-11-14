# Fluidized Bed Tutorial - Lab Scale

## Contents

- Prepares rCFD database on a 65k grid (C2C patterns for gas and solid phase).
- Runs rCFD considering (i) secondary gas species, (ii) solid mixing

## Retrieve Ansys/Fluent Files

- switch to /ansys_fluent Folder
- follow local README.md instructions to retrieve cas and dat Files
- switch back to main folder

## Prepare rCFD

- start Ansys/fluent by e.g. *fluent 3ddp -t4 &* from this folder
- read FB_tutorial.scm
- type into fluent console (rcfd_p1), followed by (rcfd_p2), ..., (rcfd_p11); or simply (rcfd_prep)

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
- monitors will be stored into ./post

- additional post-processing can be performed with matlab script ./post/monitor.m

- close Ansys/fluent

## Known Limitations

- (rcfd_run) works only after (rcfd_prep)

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

- This software is distributed under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.html).
- Copyright © 2021 JKU Linz
