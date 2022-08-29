# CFD_DEM based Field data

Request data files by developers or generate them by yourself

*Filename format:* 

	number.txt (number: 1000 to 145000 by 1000; number stands for CFD-DEM time-step)

*Data Format:*

	80.000 lines (one line for each cell)
	10 columns of ascii scientific notation numbers stating
	x, y, z, alpha, Ux, Uy, Uz, Usx, Usy, and Usz, 

	where
	
	(x,y,z) are the coordinate of the cell centers, 
	(alpha) is the void fraction, 
	(Ux, Uy, Uz) the fluid velocity and 
	(Usx, Usy, Usz) is the mass-averaged solid velocity

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

- This software is distributed under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.html).
- Copyright © 2021 JKU Linz
