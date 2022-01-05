#ifndef RCFD_USER
#define RCFD_USER

#include "rCFD_types.h"
#include "rCFD_globals.h"
#include "rCFD_defaults.h"
#include "rCFD_macros.h"

/* (C)  2021 
    Stefan Pirker
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/  

    
    /* user Names */
    /*************************************************************************************/ 

#if 1
        /* names of phases */
        enum{
            air,
            number_of_phase_names
        };

        /* names of phase user vars - phase_0 */
        enum{
            number_of_phase0_user_vars_names
        };

        /* names of data - phase_0 */
        enum{
            c_drift_0p0,
            c_drift_0p001,
            c_drift_0p002,
            number_of_phase0_data_names
        };

#endif
    
    /* user set Dicts */
    /*************************************************************************************/
    
#if 1
    
    void rCFD_user_set_Solver_Dict(void)
    {
        Solver_Dict.max_number_of_cells_per_time_step =    30;

        Solver_Dict.global_time_step =                     0.005;   
        
        Solver_Dict.time_steps_per_monitoring_interval =   25;
        
        Solver_Dict.number_of_frames =                     10;         

        Solver_Dict.number_of_runs =                       1; 
        
        Solver_Dict.data_drifting_on =                     1; 

        Solver_Dict.control_conc_sum_on =                  0;      
    }   

    void rCFD_user_set_File_Dict(void)
    {

    }   
    
    void rCFD_user_set_Phase_Dict(void)
    {
#if RP_NODE     
        int i_phase;
        
        loop_phases{
            
            if(i_phase == air){
                
                Phase_Dict[i_phase].number_of_data = 3;

                Phase_Dict[i_phase].time_step = Solver_Dict.global_time_step;
                
                Phase_Dict[i_phase].density = 1.25;
            }
        }
#endif    
    }

    void rCFD_user_set_Tracer_Dict(void)
    {

    }

    void rCFD_user_set_Norm_Dict(void)
    {

    }

    void rCFD_user_set_Rec_Dict(void)
    {

    }
    
    void rCFD_user_set_Data_Dict(void)
    {
#if RP_NODE     
        int i_phase, i_data;
        
        loop_phases{
            
            loop_data{
                
                if((i_phase == air) && (i_data == c_drift_0p0)){
                    
                    Data_Dict[i_phase][i_data].type = concentration_data;
                }

                if((i_phase == air) && (i_data == c_drift_0p001)){
                    
                    Data_Dict[i_phase][i_data].type = concentration_data;

                    Data_Dict[i_phase][i_data].drifting_data = 1;
                    
                    Data_Dict[i_phase][i_data].drift_velocity = (double*)malloc( 3 * sizeof(double));
                    
                    Data_Dict[i_phase][i_data].drift_velocity[0] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[1] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[2] = -0.015;
                }

                if((i_phase == air) && (i_data == c_drift_0p002)){
                    
                    Data_Dict[i_phase][i_data].type = concentration_data;

                    Data_Dict[i_phase][i_data].drifting_data = 1;
                    
                    Data_Dict[i_phase][i_data].drift_velocity = (double*)malloc( 3 * sizeof(double));
                    
                    Data_Dict[i_phase][i_data].drift_velocity[0] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[1] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[2] = -0.03;
                }
            }
        }
#endif    
    }
    
    void rCFD_user_set_Balance_Dict(void)
    {
#if RP_NODE     
       int i_phase, i_data;

        loop_phases{

            loop_data{
            
                /*Balance_Dict[i_phase][i_data].max_correction_loops = 5;*/

                Balance_Dict[i_phase][i_data].write_balance_to_file = 1;

                Balance_Dict[i_phase][i_data].write_balance_to_file_interval = Solver_Dict.number_of_runs;
            }
        }
#endif    
    }

    void rCFD_user_set_Topo_Dict(void)
    {
        
    }

    void rCFD_user_set_Cell_Dict(const short i_layer)
    {

    }
    
    void rCFD_user_set_Face_Dict(const short i_layer)
    {
        
    }
#endif

    /* user init */
    /*************************************************************************************/
    
#if 1   
    void rCFD_user_pre_proc(void)
    {
#if RP_NODE     
        Domain  *d=Get_Domain(1);
        Thread  *t;

        int i_phase, i_cell, i_UDMI, i_layer;
        
        i_layer = 0;
        
        thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){begin_c_loop_int(i_cell, t){
            
            i_UDMI = 0;     /* start index */
            
            loop_phases{
                
                C_UDMI(i_cell, t, i_UDMI) = _C.crossing_time[i_phase][i_cell];

                C_UDMI(i_cell, t, (i_UDMI+1)) = _C.average_velocity[i_phase][i_cell];
                
                i_UDMI += 2;
            }
            
        }end_c_loop_int(i_cell, t)}}        
#endif    
    }

    void rCFD_user_init_Data(const short i_layer)
    {
#if RP_NODE     
        int     i_phase, i_cell, i_data;

        loop_phases{
        
            loop_cells{
                             
                loop_data{
            
                    _C.data[_i_data] = 0.0; 
                }
            }
        }
        
        Message0("\n\n...rCFD_user_init_Data\n");               
#endif    
    }

#endif  
    
    /* user accesses Data */
    /*************************************************************************************/

#if 1
    
    void rCFD_user_access_data_before_shift(const short i_phase, const short i_layer)
    {
#if RP_NODE
        int     i_cell, i_data, i_dim;

        double  x[3];
        
        double  radius_in, tracer_gas_flowrate, V_in, V_in_global; 

        V_in = 0.0;
        
        tracer_gas_flowrate = 1.0e-6;   /* (kg/s) */            
        
        /* inflow bc volume */
        loop_int_cells{
            
            loop_dim{
                
                x[i_dim] = _C.x[i_cell][i_dim];
            }
            
            radius_in = sqrt(x[0]*x[0] + x[1]*x[1]);
            
            if((radius_in <= 0.01) && (x[2] < 0.075)){

                V_in += _C.volume[i_cell];
            }
        }
                
        V_in_global = PRF_GRSUM1(V_in);
        
        loop_int_cells{
            
            /* coord's */
            loop_dim{
                
                x[i_dim] = _C.x[i_cell][i_dim];
            }
            
            /* inflow bc */
            {                  
                radius_in = sqrt(x[0]*x[0] + x[1]*x[1]);
                
                if((radius_in <= 0.01) && (x[2] < 0.075)){

                    loop_data{
                        
                        _C.data[_i_data] = (_C.data[_i_data] * _C.volume[i_cell] * Phase_Dict[i_phase].density +
                        
                            _C.volume[i_cell]/V_in_global * tracer_gas_flowrate * Phase_Dict[i_phase].time_step) /
                            
                            (_C.volume[i_cell] * Phase_Dict[i_phase].density);  /* (.) */

                        Balance[i_phase][i_data].mass_source +=  _C.volume[i_cell]/V_in_global * tracer_gas_flowrate * Phase_Dict[i_phase].time_step;  /* (kg) */
                    }
                }   
            }
            
            /* outflow bc */
            {
                if(x[0] > 0.6){

                    loop_data{
                        
                        Balance[_i_balance].mass_source -=  _C.data[_i_data] * _C.volume[i_cell] * Phase_Dict[i_phase].density;  /* (kg) */

                        _C.data[_i_data] = 0.0;  /* (kg) */
                    }
                }   
            }
        }           
#endif
    }

    void rCFD_user_access_data_after_swap(const short i_phase, const short i_layer)
    {

    }

    void rCFD_user_post(const short i_layer)
    {
#if RP_NODE     
        Domain  *d=Get_Domain(1);
        Thread  *t;

        int i_phase, i_cell, i_data, i_UDMI;
        
        thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){begin_c_loop_int(i_cell, t){
            
            i_UDMI = 0;     /* start index */
            
            loop_phases{
                
                loop_data{
                
                    C_UDMI(i_cell, t, i_UDMI) = _C.data[_i_data];
                    
                    i_UDMI++;
                }
            }
            
        }end_c_loop_int(i_cell, t)}}        
#endif    
    }
#endif

    /* user sub models */
    /*************************************************************************************/

#if 1
    
    double rCFD_user_set_random_walk_velocity(void)
    {
        return 0.0;
    }   

    void rCFD_user_set_Norm(void)
    {
        
    }
    
#endif  

#endif