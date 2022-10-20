#include "udf.h"

static short 	_monitor_file_initiated	= 0;

static int		_monitor_step = 0;

static double	_average_solid_conc = 0.0;

#define			_index_of_gas_phase		0
#define			_index_of_primary_gas   0
#define 		_index_of_secondary_gas 1

#define			_index_of_solid_phase	1
#define			_index_of_mixing_solid	0

#define			_monitor_reset_interval ((int)(15./0.0005))

/*************************************************************************************/
DEFINE_EXECUTE_AT_END(CFD_reset)
/*************************************************************************************/
{
	
#if 1	/* local vars */

	int  	i_cell;
	
	double	x[3], mass, mass_conc;

	Domain  *d=Get_Domain(1);
	
	Thread  *t = NULL, *tg = NULL, *ts = NULL;
	
#endif	

	if((_monitor_step % _monitor_reset_interval) == 0){
		
		/* 1. initialize conc fields */
		
		thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){
			
			tg = THREAD_SUB_THREAD(t, _index_of_gas_phase);
			
			ts = THREAD_SUB_THREAD(t, _index_of_solid_phase);
			
			begin_c_loop_int(i_cell, t){
			
				C_YI(i_cell, tg, _index_of_secondary_gas) = 0.0;
				
				C_YI(i_cell, tg, _index_of_primary_gas) = 1.0;
				
				C_CENTROID(x, i_cell, t);
				
				if(x[1] > 0.0){
					
					C_YI(i_cell, ts, _index_of_mixing_solid) = 1.0;
				}
				else{
					
					C_YI(i_cell, ts, _index_of_mixing_solid) = 0.0;
				}					

			}end_c_loop_int(i_cell, t)
		
		}}
		
		/* 2. initialize _average_solid_conc */
		{
			mass = 0.0;
			
			mass_conc = 0.0;
			
			thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){
				
				ts = THREAD_SUB_THREAD(t, _index_of_solid_phase);
				
				begin_c_loop_int(i_cell, t){
					
					mass += C_R(i_cell, ts) * C_VOF(i_cell, ts) * C_VOLUME(i_cell, t);
					
					mass_conc += C_YI(i_cell, ts, _index_of_mixing_solid) * 
					
						C_R(i_cell, ts) * C_VOF(i_cell, ts) * C_VOLUME(i_cell, t);
					
				}end_c_loop_int(i_cell, t)
				
			}}
			
			mass = PRF_GRSUM1(mass);
			
			mass_conc = PRF_GRSUM1(mass_conc);				
			
			if(mass > 0.0){
				
				_average_solid_conc = mass_conc / mass;
			}
			else{
				
				_average_solid_conc = 0.0;
			}
		}
		
		Message0("\n\nCFD_reset\n");
	}
}

/*************************************************************************************/
DEFINE_EXECUTE_AT_END(CFD_monitor)
/*************************************************************************************/
{
	
#if 1	/* local vars */

	int  	i_cell;
	
	double	mass_of_secondary_gas, mixing_nom, mixing_denom, solid_mixing_index;

	Domain  *d=Get_Domain(1);
	
	Thread  *t = NULL, *tg = NULL, *ts = NULL;

	FILE	*f_out = NULL;	
	
#endif	

	/* 1. calculate mass_of_secondary_gas */
	{
		mass_of_secondary_gas = 0.0;
		
		thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){
			
			tg = THREAD_SUB_THREAD(t, _index_of_gas_phase);
			
			begin_c_loop_int(i_cell, t){
			
				mass_of_secondary_gas += C_YI(i_cell, tg, _index_of_secondary_gas) * 
				
					C_R(i_cell, tg) * C_VOF(i_cell, tg) * C_VOLUME(i_cell, t);
			
			}end_c_loop_int(i_cell, t)
			
		}}

		mass_of_secondary_gas = PRF_GRSUM1(mass_of_secondary_gas);
	}

	/* 2. calculate solid_mixing_index */
	{
		
		mixing_nom = 0.0;
		
		mixing_denom = 0.0;
		
		thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){
			
			ts = THREAD_SUB_THREAD(t, _index_of_solid_phase);
			
			begin_c_loop_int(i_cell, t){
				
				mixing_nom += fabs(C_YI(i_cell, ts, _index_of_mixing_solid) - _average_solid_conc) * 
				
					C_R(i_cell, ts) * C_VOF(i_cell, ts) * C_VOLUME(i_cell, t);
					
				mixing_denom += C_R(i_cell, ts) * C_VOF(i_cell, ts) * C_VOLUME(i_cell, t);

			}end_c_loop_int(i_cell, t)
			
		}}
		
		mixing_nom = PRF_GRSUM1(mixing_nom);
		
		mixing_denom = PRF_GRSUM1(mixing_denom);
		
		if(mixing_denom > 0.0){
			
			solid_mixing_index = 1.0 - 2.0 * mixing_nom / mixing_denom;
			
			/* mixing index should be within the value of 0 and 1 */
		}
		else{
			
			solid_mixing_index = 0.0;
		}			
	}
	
	/* 3. write monitor file */
	{
		if(myid == 0){
			
			if(_monitor_file_initiated == 0){
				
				f_out = fopen("./monitor_CFD.out", "w");
				
				_monitor_file_initiated = 1;
			}
			else{

				f_out = fopen("./monitor_CFD.out", "a");	
			}
			
			if(!f_out){
				
				Message("\nERROR CFD_monitor: Could not open ./monitor_CFD.out");
				
				return;
			}
			
			fprintf(f_out, "%e %e %e\n", CURRENT_TIME, mass_of_secondary_gas, solid_mixing_index); 
			
			fclose(f_out);
		}
	}
	
	_monitor_step ++;
	
    Message0("\n\nCFD_monitor\n");
}

