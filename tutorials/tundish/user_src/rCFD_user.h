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
            steel,
            number_of_phase_names
        };

        /* names of data - phase_0 */
        enum{
            c_drift_0p0,
            c_drift_0p001,
            c_drift_0p002,
            c_steel_grade_A,
            c_steel_grade_B,
            c_steel_temp,
            number_of_phase0_data_names
        };
#endif
    
    /* user set Dicts */
    /*************************************************************************************/
    
#if 1
    
    void rCFD_user_set_Solver_Dict(void)
    {
        Solver_Dict.max_number_of_cells_per_time_step =    30;

        Solver_Dict.global_time_step =                     0.02;
        
        Solver_Dict.time_steps_per_monitoring_interval =   20; 
        
        Solver_Dict.number_of_frames =                     10;         

        Solver_Dict.number_of_runs =                       50;
        
        Solver_Dict.data_drifting_on =                     1;   

        Solver_Dict.balance_correction_on =                1;
        
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
            
            if(i_phase == steel){
                
                Phase_Dict[i_phase].number_of_data = 6;

                Phase_Dict[i_phase].time_step = Solver_Dict.global_time_step;
                
                Phase_Dict[i_phase].density = 1000.0;

                Phase_Dict[i_phase].heat_capacity = 4182.0;
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
                
                if((i_phase == steel) && (i_data == c_drift_0p0)){
                    
                    Data_Dict[i_phase][i_data].type = concentration_data;
                }

                if((i_phase == steel) && (i_data == c_drift_0p001)){
                    
                    Data_Dict[i_phase][i_data].type = concentration_data;

                    Data_Dict[i_phase][i_data].drifting_data = 1;
                    
                    Data_Dict[i_phase][i_data].drift_velocity = (double*)malloc( 3 * sizeof(double));
                    
                    Data_Dict[i_phase][i_data].drift_velocity[0] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[1] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[2] = 0.001;
                }

                if((i_phase == steel) && (i_data == c_drift_0p002)){
                    
                    Data_Dict[i_phase][i_data].type = concentration_data;

                    Data_Dict[i_phase][i_data].drifting_data = 1;
                    
                    Data_Dict[i_phase][i_data].drift_velocity = (double*)malloc( 3 * sizeof(double));
                    
                    Data_Dict[i_phase][i_data].drift_velocity[0] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[1] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[2] = 0.002;
                }

                if((i_phase == steel) && (i_data == c_steel_grade_A)){
                    
                    Data_Dict[i_phase][i_data].type = concentration_data;
                }

                if((i_phase == steel) && (i_data == c_steel_grade_B)){
                    
                    Data_Dict[i_phase][i_data].type = concentration_data;
                }

                if((i_phase == steel) && (i_data == c_steel_temp)){
                    
                    Data_Dict[i_phase][i_data].type = temperature_data;
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
        int     i_phase, i_cell, i_data, i_dim;

        double  x[3], radius;

        i_phase = steel;
        
        loop_cells{
            
            /* coord's */
            loop_dim{
                
                x[i_dim] = _C.x[i_cell][i_dim];
            }
            
            radius = sqrt((x[0] - 0.935)*(x[0] - 0.935) + x[1]*x[1]);
            
            /* Field Data */
            loop_data{
        
                /* global initialization */
                if(i_data == c_steel_grade_A){
                    
                    _C.data[_i_data] = 1.0;
                }
                else{
                    
                    _C.data[_i_data] = 0.0;
                }                   
                
                /* inflow initialization */
                if((radius <= 0.0115) && (x[2] > 0.2)){
                    
                    switch (i_data){
                            
                        case c_drift_0p0:
                        {
                            _C.data[_i_data] = 1.0e-3;
                            
                            break;
                        }
                        case c_drift_0p001:
                        {
                            _C.data[_i_data] = 1.0e-3;
                            
                            break;
                        }
                        case c_drift_0p002:
                        {
                            _C.data[_i_data] = 1.0e-3;
                            
                            break;
                        }
                        case c_steel_grade_A:
                        {
                            _C.data[_i_data] = 0.0;
                            
                            break;
                        }
                        case c_steel_grade_B:
                        {
                            _C.data[_i_data] = 1.0;
                            
                            break;
                        }
                        case c_steel_temp:
                        {
                            _C.data[_i_data] = 23;
                            
                            break;
                        }                       
                    }               
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

        double  x[3], radius_in, radius_out, flowrate;
        
        double  mean_value_in[Phase_Dict[i_phase].number_of_data], V_in, V_in_global;
        
        double  mean_value_out[Phase_Dict[i_phase].number_of_data], V_out, V_out_global; 

        V_in = 0.0; V_out = 0.0; 
        
        flowrate = 0.415;   /* (kg/s) */

        /* init mean values */
        loop_data{
            
            mean_value_in[i_data] = 0.0;
            
            mean_value_out[i_data] = 0.0;
            
            Balance[i_phase][i_data].mass_source = 0.0;
        }
        
        /* volumes, constant and mean values at bc */ 
        {
            loop_int_cells{
                
                loop_dim{
                    
                    x[i_dim] = _C.x[i_cell][i_dim];
                }

                radius_in = sqrt((x[0] - 0.935)*(x[0] - 0.935) + x[1]*x[1]);

                if((radius_in <= 0.0115) && (x[2] > 0.2)){

                    V_in += _C.volume[i_cell];
                     
                    /* constant concentrations */
                    loop_data{
                        
                        switch (i_data){
                            
                            case c_drift_0p0:
                            
                                _C.data[_i_data] = 1.0e-3; break;

                            case c_drift_0p001:
                            
                                _C.data[_i_data] = 1.0e-3; break;

                            case c_drift_0p002:
                            
                                _C.data[_i_data] = 1.0e-3; break;
                                
                            case c_steel_grade_A:
                            
                                _C.data[_i_data] = 0.0; break;
                            
                            case c_steel_grade_B:
                            
                                _C.data[_i_data] = 1.0; break;
                            
                            case c_steel_temp:
                            
                                _C.data[_i_data] = 23.0; break;
                            
                            default: break;
                        }

                        mean_value_in[i_data] += _C.data[_i_data] * _C.volume[i_cell];                                      
                    }
                }

                radius_out = sqrt((x[0] - 0.085)*(x[0] - 0.085) + x[1]*x[1]);
                
                if((radius_out <= 0.015) && (x[2] < 0.025)){

                    V_out += _C.volume[i_cell];

                    /* mean values */
                    loop_data{
                        
                        mean_value_out[i_data] += _C.data[_i_data] * _C.volume[i_cell];                 
                    }
                    
                }           
            }
            
            V_in_global = PRF_GRSUM1(V_in);
            
            V_out_global = PRF_GRSUM1(V_out);
            
            if(V_in > 0.0){
                
                loop_data{
                    
                    mean_value_in[i_data] /= V_in;
                }
            }

            if(V_out > 0.0){
                
                loop_data{
                    
                    mean_value_out[i_data] /= V_out;
                }
            }
        }
                
        /* set balance.sources */
        {    
            if(V_in > 0){
                
                loop_data{
        
                    if(i_data < c_steel_temp){
                        
                        Balance[_i_balance].mass_source += V_in/V_in_global * mean_value_in[i_data] * 
                        
                            flowrate * Phase_Dict[i_phase].time_step;  /* (kg) */
                    }
                    
                    if(i_data == c_steel_temp){
                        
                        Balance[_i_balance].mass_source += V_in/V_in_global * mean_value_in[i_data] * 
                        
                            flowrate * Phase_Dict[i_phase].heat_capacity * Phase_Dict[i_phase].time_step;  /* (J) */
                    }
                }
            }           

            if(V_out > 0){
                
                loop_data{
                            
                    if(i_data < c_steel_temp){
                        
                        Balance[_i_balance].mass_source -= V_out/V_out_global * mean_value_out[i_data] * 
                        
                            flowrate * Phase_Dict[i_phase].time_step;  /* (kg) */
                    }
                    
                    if(i_data == c_steel_temp){
                        
                        Balance[_i_balance].mass_source -= V_out/V_out_global * mean_value_out[i_data] * 
                        
                            flowrate * Phase_Dict[i_phase].heat_capacity * Phase_Dict[i_phase].time_step;  /* (J) */
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
                    
                    if(i_data == c_steel_temp){
                        
                        C_UDMI(i_cell, t, i_UDMI) += Phase_Dict[i_phase].reference_temperature;
                    }
                    
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