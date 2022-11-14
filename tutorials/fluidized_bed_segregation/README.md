# Fluidized Bed Segregation Tutorial

## Contents

- Converts frames of CFD-DEM based field information into Fluent *.ip format
- Generates rCFD database
- Runs rCFD in replay mode, considering segregation of two particle fractions (small/large)

## Retrieve Ansys/Fluent cas and dat Files

- switch to folder ./ansys_fluent
- retrieve files by following the links given in the local README file
- switch back

## Retrieve CFD-DEM data

- switch to folder ./data/csv
- retrieve compressed zip-file by following the link given in the local README file
- extract zip-file into this folder
- switch back

## Convert CFD-DEM based field data

- start Ansys/fluent by *fluent 3ddp -t4 &* from this folder
- from the file menu read the scheme tutorial.scm
- type into fluent console (cfd_conv1), followed by (cfd_conv2), (cfd_conv3); or simply (cfd_conv)

- close ANSYS/fluent

## Prepare rCFD

- start Ansys/fluent by *fluent 3ddp -t4 &* from this folder
- from the file menu read the scheme tutorial.scm
- type into fluent console (rcfd_p1), followed by (rcfd_p2), ..., (rcfd_p11); or simply (rcfd_prep)

- c2c patterns will be stored into ./data/c2c
- jump vector will be stored into ./rec

- close Ansys/fluent

## Run rCFD in GUI mode

- start Ansys/fluent by *fluent 3ddp -t4 &* from this folder
- from the file menu read the scheme tutorial.scm
- type into fluent console (rcfd_r1), followed by (rcfd_r2), ..., (rcfd_r7); or simply (rcfd_run)

- balances will be stored into ./post
- monitors will be stored into ./post
- run transcript will be stored into ./post
- animation snap-shots will be stored into ./post

- close Ansys/fluent

- ./post/monitor.m post-processes monitor_rCFD.out
- ./post/reference contains corresponding monitors of full CFD-DEM simulations


## Run rCFD in batch mode

- type *fluent 3ddp -t4 -g  < run_batch.scm >& run_batch.trn* into your C-shell

- balances will be stored into ./post
- monitors will be stored into ./post
- transcripts will be stored into ./post

## Known Limitations

- (rcfd_run) works only after (rcfd_prep)
- Ansys/fluent has to be closed between *prepare rCFD* and *run rCFD*

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

- This software is distributed under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.html).
- Copyright © 2021 JKU Linz
