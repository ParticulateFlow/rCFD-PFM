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

        /* names of phase user vars - phase_0 */
        enum{
            number_of_phase0_user_vars_names
        };

        /* names of data - phase_0 */
        enum{
            c_drift_0p0,
            c_drift_0p001,
            c_drift_0p002,
            c_steel_grade_A,
            c_steel_grade_B,
            c_temp,
            number_of_phase0_data_names
        };

#endif
    
    /* user set Dicts */
    /*************************************************************************************/
    
#if 1
    
    void rCFD_user_set_Solver_Dict(Solver_Dict_type *Solver_Dict)
    {
        Solver_Dict->max_number_of_cells_per_time_step =    30;
        
        Solver_Dict->number_of_frames =                     10;         

        Solver_Dict->number_of_runs =                       50; 
        
        Solver_Dict->data_drifting_on =                     1;      
    }   

    void rCFD_user_set_File_Dict(Solver_Dict_type *Solver_Dict, File_Dict_type *File_Dict)
    {

    }   
    
    void rCFD_user_set_Phase_Dict(Solver_Dict_type *Solver_Dict, Phase_Dict_type* Phase_Dict)
    {
        int i_phase;
        
        loop_phases_ptr{
            
            if(i_phase == steel){
                
                Phase_Dict[i_phase].number_of_data = 6;

                Phase_Dict[i_phase].time_step = 0.02;
                
                Phase_Dict[i_phase].density = 1000.0;

                Phase_Dict[i_phase].heat_capacity = 4182.0;
            }
        }
    }

    void rCFD_user_set_Tracer_Dict(Solver_Dict_type *Solver_Dict, Tracer_Dict_type *Tracer_Dict)
    {

    }

    void rCFD_user_set_Norm_Dict(Solver_Dict_type *Solver_Dict, Norm_Dict_type *Norm_Dict)
    {

    }

    void rCFD_user_set_Rec_Dict(Solver_Dict_type *Solver_Dict, Rec_Dict_type *Rec_Dict)
    {

    }
    
    void rCFD_user_set_Data_Dict(Solver_Dict_type *Solver_Dict, Phase_Dict_type *Phase_Dict, Data_Dict_type **Data_Dict)
    {
        int i_phase, i_data;
        
        loop_phases_ptr{
            
            loop_data{
                
                if((i_phase == steel) && (i_data == c_drift_0p0)){
                    
                    Data_Dict[i_phase][i_data].type = generic_data;
                }

                if((i_phase == steel) && (i_data == c_drift_0p001)){
                    
                    Data_Dict[i_phase][i_data].type = generic_data;

                    Data_Dict[i_phase][i_data].drifting_data = 1;
                    
                    Data_Dict[i_phase][i_data].drift_velocity = (double*)malloc( 3 * sizeof(double));
                    
                    Data_Dict[i_phase][i_data].drift_velocity[0] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[1] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[2] = 0.001;
                }

                if((i_phase == steel) && (i_data == c_drift_0p002)){
                    
                    Data_Dict[i_phase][i_data].type = generic_data;

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

                if((i_phase == steel) && (i_data == c_temp)){
                    
                    Data_Dict[i_phase][i_data].type = temperature_data;
                }

            }
        }
    }
    
    void rCFD_user_set_Balance_Dict(Solver_Dict_type *Solver_Dict, Phase_Dict_type *Phase_Dict, Balance_Dict_type **Balance_Dict)
    {
       int i_phase, i_data;

        loop_phases_ptr{

            loop_data{
            
                Balance_Dict[i_phase][i_data].write_balance_to_file = 1;

                Balance_Dict[i_phase][i_data].write_balance_to_file_interval = Solver_Dict->number_of_runs;
            }
        }
    }

    void rCFD_user_set_Topo_Dict(Solver_Dict_type *Solver_Dict, Topo_Dict_type *Topo_Dict)
    {
        
    }

    void rCFD_user_set_Cell_Dict(Solver_Dict_type *Solver_Dict, Cell_Dict_type *Cell_Dict, const short i_layer)
    {

    }
    
    void rCFD_user_set_Face_Dict(Solver_Dict_type *Solver_Dict, Face_Dict_type *Face_Dict, const short i_layer)
    {
        
    }
#endif

    /* user init */
    /*************************************************************************************/
    
#if 1   
    void rCFD_user_pre_proc(Solver_Dict_type *Solver_Dict, Phase_Dict_type *Phase_Dict, Topo_Dict_type *Topo_Dict, Cell_type *C)
    {
        Domain  *d=Get_Domain(1);
        Thread  *t;

        int i_phase, i_cell, i_UDMI;
        
        thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){begin_c_loop_int(i_cell, t){
            
            i_UDMI = 0;     /* start index */
            
            loop_phases_ptr{
                
                C_UDMI(i_cell, t, i_UDMI) = C->crossing_time[i_phase][i_cell];

                C_UDMI(i_cell, t, (i_UDMI+1)) = C->average_velocity[i_phase][i_cell];
                
                i_UDMI += 2;
            }
            
        }end_c_loop_int(i_cell, t)}}        
    }

    void rCFD_user_set_recurrence_time_step(Solver_Dict_type *Solver_Dict, Phase_Dict_type* Phase_Dict)
    {
        
    }
        
    void rCFD_user_init_Data(Solver_Dict_type *Solver_Dict, Balance_type** Balance, Phase_Dict_type* Phase_Dict, Data_Dict_type** Data_Dict, 
        Topo_Dict_type *Topo_Dict, Cell_type *C, const short i_layer)
    {
        int     i_phase, i_cell, i_data, i_dim;

        double  x[3], radius;

        loop_phases_ptr{
        
            loop_cells_ptr{
                
                /* coord's */
                loop_dim{
                    
                    x[i_dim] = C->x[i_cell][i_dim];
                }
                
                radius = sqrt((x[0] - 0.935)*(x[0] - 0.935) + x[1]*x[1]);
                
                /* Field Data */
                loop_data{
            
                    if((i_phase == steel) && (i_data == c_drift_0p0)){
                    
                        C->data[i_phase][i_cell][i_data] = 0.0;
                        
                        /* mark incoming tube */
                        if((radius <= 0.0115) && (x[2] > 0.2)){
                            
                            C->data[i_phase][i_cell][i_data] = 0.5;
                        }
                    }

                    if((i_phase == steel) && (i_data == c_drift_0p002)){
                    
                        C->data[i_phase][i_cell][i_data] = 0.0;
                        
                        /* mark incoming tube */
                        if((radius <= 0.0115) && (x[2] > 0.2)){
                            
                            C->data[i_phase][i_cell][i_data] = 0.5;
                        }
                    }
                    
                    if((i_phase == steel) && (i_data == c_drift_0p001)){
                    
                        C->data[i_phase][i_cell][i_data] = 0.0;
                        
                        /* mark incoming tube */
                        if((radius <= 0.0115) && (x[2] > 0.2)){
                            
                            C->data[i_phase][i_cell][i_data] = 0.5;
                        }
                    }
                    
                    if((i_phase == steel) && (i_data == c_steel_grade_A)){
                    
                        C->data[i_phase][i_cell][i_data] = 1.0;
                        
                        /* mark incoming tube */
                        if((radius <= 0.0115) && (x[2] > 0.2)){
                            
                            C->data[i_phase][i_cell][i_data] = 0.0;
                        }
                    }
                    
                    if((i_phase == steel) && (i_data == c_steel_grade_B)){
                    
                        C->data[i_phase][i_cell][i_data] = 0.0;
                        
                        /* mark incoming tube */
                        if((radius <= 0.0115) && (x[2] > 0.2)){
                            
                            C->data[i_phase][i_cell][i_data] = 1.0;
                        }
                    }
                    
                    if((i_phase == steel) && (i_data == c_temp)){
                    
                        C->data[i_phase][i_cell][i_data] = 300.0;
                        
                        /* mark incoming tube */
                        if((radius <= 0.0115) && (x[2] > 0.2)){
                            
                            C->data[i_phase][i_cell][i_data] = 323.0;
                        }
                    }                   
                }
            }
        }
        
        Message0("\n\n...rCFD_user_init_Data\n");               
    }

#endif  
    
    /* user accesses Data */
    /*************************************************************************************/

#if 1
    
    void rCFD_user_access_data_before_shift(Balance_type** Balance, Phase_Dict_type* Phase_Dict, 
        Topo_Dict_type *Topo_Dict, Cell_type *C, Rec_type *Rec, 
        const short i_phase, const short i_layer)
    {
#if RP_NODE
        int     i_cell, i_data, i_dim;

        double  x[3], radius_in, radius_out, flowrate;
        
        double  mean_value_in[Phase_Dict->number_of_data], V_in, V_in_global;
        
        double  mean_value_out[Phase_Dict->number_of_data], V_out, V_out_global; 

        V_in = 0.0; 
        V_out = 0.0;
        
        flowrate = 0.415;   /* (kg/s) */

        /* init mean values */
        loop_data{
            
            mean_value_in[i_data] = 0.0;
            
            mean_value_out[i_data] = 0.0;
        }
        
        loop_cells_ptr{
            
            /* coord's */
            loop_dim{
                
                x[i_dim] = C->x[i_cell][i_dim];
            }
            
            /* inflow bc */
            {
                radius_in = sqrt((x[0] - 0.935)*(x[0] - 0.935) + x[1]*x[1]);
                
                if((radius_in <= 0.0115) && (x[2] > 0.2)){
                    
                    loop_data{
                
                        /* set inflow bc value for i_data */
                        switch (i_data){
                            
                            case c_drift_0p0:
                            
                                C->data[i_phase][i_cell][i_data] = 0.5; break;
                            
                            case c_drift_0p001:
                            
                                C->data[i_phase][i_cell][i_data] = 0.5; break;
                            
                            case c_drift_0p002:
                            
                                C->data[i_phase][i_cell][i_data] = 0.5; break;
                            
                            case c_steel_grade_A:
                            
                                C->data[i_phase][i_cell][i_data] = 0.0; break;
                            
                            case c_steel_grade_B:
                            
                                C->data[i_phase][i_cell][i_data] = 1.0; break;
                            
                            case c_temp:
                            
                                C->data[i_phase][i_cell][i_data] = 323.0; break;
                            
                            default: break;
                        }
                        
                        /* calc mean bc values (accum. values and volume */
                        mean_value_in[i_data] += C->data[i_phase][i_cell][i_data] * C->volume[i_cell];
                    }
                    
                    V_in += C->volume[i_cell];
                }   
            }
            
            /* outflow bc */
            {
                radius_out = sqrt((x[0] - 0.085)*(x[0] - 0.085) + x[1]*x[1]);
                
                if((radius_out <= 0.015) && (x[2] < 0.025)){
                        
                    loop_data{
                
                        /* calc mean bc values (accum. values and volume */
                        mean_value_out[i_data] += C->data[i_phase][i_cell][i_data] * C->volume[i_cell];
                    }
                    
                    V_out += C->volume[i_cell];
                }
            }
            
            /* heat losses, TODO */
            {
            }
        }

        V_in_global = PRF_GRSUM1(V_in); 
        
        V_out_global = PRF_GRSUM1(V_out);
        
        if(V_in > 0){
            
            loop_data{
                        
                mean_value_in[i_data] /= V_in;
                    
                /* per partition source */  
                if(i_data < c_temp){
                    
                    Balance[i_phase][i_data].mass_source += V_in/V_in_global * flowrate * mean_value_in[i_phase] * 
                    
                        Phase_Dict[i_phase].time_step;  /* (kg) */
                }
                else{
                    
                    Balance[i_phase][i_data].mass_source += V_in/V_in_global * flowrate * mean_value_in[i_phase] * 
                    
                        Phase_Dict[i_phase].heat_capacity * Phase_Dict[i_phase].time_step;  /* (J) */
                }
            }
        }           

        if(V_out > 0){
            
            loop_data{
                        
                mean_value_out[i_data] /= V_out;
                    
                /* per partition source */  
                if(i_data < c_temp){
                    
                    Balance[i_phase][i_data].mass_source += V_out/V_out_global * flowrate * mean_value_out[i_phase] * 
                    
                        Phase_Dict[i_phase].time_step;  /* (kg) */
                }
                else{
                    
                    Balance[i_phase][i_data].mass_source += V_out/V_out_global * flowrate * mean_value_out[i_phase] * 
                    
                        Phase_Dict[i_phase].heat_capacity * Phase_Dict[i_phase].time_step;  /* (J) */
                }
            }
        }           
#endif
    }

    void rCFD_user_access_data_after_swap(Balance_type** Balance, Phase_Dict_type* Phase_Dict, 
        Data_Dict_type **Data_Dict, Topo_Dict_type *Topo_Dict, Cell_type *C, 
        const short i_phase, const short i_layer)
    {

    }

    void rCFD_user_post(Solver_Dict_type *Solver_Dict, Phase_Dict_type *Phase_Dict, Topo_Dict_type *Topo_Dict, Cell_type *C)
    {
        Domain  *d=Get_Domain(1);
        Thread  *t;

        int i_phase, i_cell, i_data, i_UDMI;
        
        thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){begin_c_loop_int(i_cell, t){
            
            i_UDMI = 0;     /* start index */
            
            loop_phases_ptr{
                
                loop_data{
                
                    C_UDMI(i_cell, t, i_UDMI) = C->data[i_phase][i_cell][i_data];
                    
                    i_UDMI++;
                }
            }
            
        }end_c_loop_int(i_cell, t)}}        
    }
#endif

    /* user sub models */
    /*************************************************************************************/

#if 1
    
    double rCFD_user_set_random_walk_velocity(Solver_Dict_type *Solver_Dict, Thread *t, int i_cell)
    {
        return 0.0;
    }   

    void rCFD_user_set_Norm(Solver_Dict_type *Solver_Dict, Norm_type *Norms)
    {
        
    }
    
#endif  

#endif