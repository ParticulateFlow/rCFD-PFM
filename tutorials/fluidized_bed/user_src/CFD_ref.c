#include <udf.h>

/* (C)  2022 
    Stefan Pirker
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/  

#if RP_NODE
static short	ref_species_initiated = 0;

static double	mean_solid_A_conc = 0.0;

enum{	
	gas_A,
	gas_B
};

enum{	
	solid_A,
	solid_B
};
#endif

/*************************************************************************************/
DEFINE_EXECUTE_AT_END(CFD_ref_monitors)
/*************************************************************************************/
{
#if RP_NODE

#if 1	/* local vars */ 

    Domain      *d = Get_Domain(1);
    Thread      *t, *t_gas, *t_solid;
    cell_t      c;
  
    double      x[3], mass_gas_B, mixing_nom, mixing_denom, mixing_index; 
	
#endif	

	/* R.1. initialize conc. at start */
	if(ref_species_initiated == 0){
			
		/* R.1.1. init gas species */
		thread_loop_c(t,d){
			
			t_gas = THREAD_SUB_THREAD(t,0);
			
			begin_c_loop_int(c,t){
			
				C_CENTROID(x,c,t);

				if(x[1] >= -0.074375){
					
					C_YI(c,t_gas, gas_A) = 1.0;
					C_YI(c,t_gas, gas_B) = 0.0;
				}
				else{
					
					C_YI(c,t_gas, gas_A) = 0.0;
					C_YI(c,t_gas, gas_B) = 1.0;
				}
			
		}end_c_loop_int(c,t)}			
		
		/* R.1.2. init solid species */
		thread_loop_c(t,d){
			
			t_solid = THREAD_SUB_THREAD(t,1);
			
			begin_c_loop_int(c,t){
			
				C_CENTROID(x,c,t);

				if(x[1] >= 0.0){
					
					C_YI(c,t_solid, solid_A) = 1.0;
					C_YI(c,t_solid, solid_B) = 0.0;
				}
				else{
					
					C_YI(c,t_solid, solid_A) = 0.0;
					C_YI(c,t_solid, solid_B) = 1.0;
				}
			
		}end_c_loop_int(c,t)}			

		/* R.1.3. calc mean_solid_A_conc */
		{
			mixing_nom = 0.0;
			
			mixing_denom = 0.0;

			thread_loop_c(t,d){
				
				t_solid = THREAD_SUB_THREAD(t,1);
				
				begin_c_loop_int(c,t){
				
					mixing_nom += C_YI(c,t_solid, solid_A) * C_VOLUME(c,t) * C_VOF(c,t_solid) * C_R(c, t_solid);

					mixing_denom += C_VOLUME(c,t) * C_VOF(c,t_solid) * C_R(c, t_solid);
					
			}end_c_loop_int(c,t)}			

			mixing_nom = PRF_GRSUM1(mixing_nom);
			
			mixing_denom = PRF_GRSUM1(mixing_denom);
			
			if(mixing_denom > 0.0){
				
				mean_solid_A_conc = mixing_nom / mixing_denom;
			}
		}
	}
	
	/* R.2. calc m_gas_B, mixing_solid */
	if(ref_species_initiated){
		
		mass_gas_B = 0.0;

		thread_loop_c(t,d){

			t_gas = THREAD_SUB_THREAD(t,0);
			
			begin_c_loop_int(c,t){
			
				mass_gas_B += C_YI(c,t_gas, gas_B) * C_VOLUME(c,t) * C_VOF(c,t_gas) * C_R(c, t_gas);
			
		}end_c_loop_int(c,t)}			
		
		mass_gas_B = PRF_GRSUM1(mass_gas_B);
		
		mixing_nom = 0.0;
		
		mixing_denom = 0.0;
		
		thread_loop_c(t,d){

			t_solid = THREAD_SUB_THREAD(t,1);
			
			begin_c_loop_int(c,t){
			
				mixing_nom += fabs(C_YI(c,t_solid, solid_A) - mean_solid_A_conc) * C_VOLUME(c,t) * C_VOF(c,t_solid) * C_R(c, t_solid);
				
				mixing_denom += C_VOLUME(c,t) * C_VOF(c,t_solid) * C_R(c, t_solid);				
			
		}end_c_loop_int(c,t)}		

		mixing_nom = PRF_GRSUM1(mixing_nom);
		
		mixing_denom = PRF_GRSUM1(mixing_denom);
		
		if(mixing_denom > 0.0){
			
			mixing_index = 1.0 - 2.0 * mixing_nom / mixing_denom;
		}
		else{
			
			mixing_index = 0.0;
		}
	}
	
	/* R.3. node-0 writes values to ref file */
	if((ref_species_initiated) && (myid == 0)){

		FILE	*f_out = NULL;
		
		if(ref_species_initiated == 0){
			
			f_out = fopen("./monitor.out","w");
		}
		else{
			f_out = fopen("./monitor.out","a");
		}
		
		if(f_out == NULL){
			
			Message0("\nERROR: Could not open ref monitor file\n");
			
			return;
		}
		
		fprintf(f_out, "%e %e %e\n", CURRENT_TIME, mass_gas_B, mixing_index);
		
		fclose(f_out);
	}
		
	/* R.4. write species conc to UDMI */
	if(ref_species_initiated){
		
		thread_loop_c(t,d){
			
			t_gas = THREAD_SUB_THREAD(t,0);			
			
			begin_c_loop_int(c,t){
			
				C_UDMI(c,t, 1) = C_YI(c,t_gas, gas_B) * C_VOF(c, t_gas);
			
		}end_c_loop_int(c,t)}	
		
		thread_loop_c(t,d){
			
			t_solid = THREAD_SUB_THREAD(t,1);			
			
			begin_c_loop_int(c,t){
			
				C_UDMI(c,t, 2) = C_YI(c,t_solid, solid_A) * C_VOF(c, t_solid);
			
		}end_c_loop_int(c,t)}	
		
	}
	
	ref_species_initiated = 1;
	
	Message0("\n\n...CFD_ref_monitors\n");
	
#endif
}
