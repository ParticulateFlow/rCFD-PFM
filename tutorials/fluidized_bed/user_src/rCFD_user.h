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
            gas,
            solid,
            number_of_phase_names
        };

        /* names of data - phase_0 */
        enum{
            c_gas_A,
            c_gas_B,
            number_of_phase0_data_names
        };

        /* names of data - phase_1 */
        enum{
            c_drift_0,
            c_drift_z_pos,
            c_drift_z_neg,
            number_of_phase1_data_names
        };

#endif
    
    /* user set Dicts */
    /*************************************************************************************/
    
#if 1
    
    void rCFD_user_set_Solver_Dict(Solver_Dict_type *Solver_Dict)
    {
        Solver_Dict->number_of_phases =                     2;
        
        Solver_Dict->max_number_of_cells_per_time_step =    30;

        Solver_Dict->global_time_step =                     0.005;   
        
        Solver_Dict->time_steps_per_monitoring_interval =   10;
        
        Solver_Dict->number_of_frames =                     300;         

        Solver_Dict->number_of_runs =                       1; 
        
        Solver_Dict->data_drifting_on =                     1;

        Solver_Dict->face_diffusion_on =                    0; 
        
        Solver_Dict->recurrence_process_on =                1;  

        Solver_Dict->face_swap_max_per_loop =               (1./8.);        
        
    }   

    void rCFD_user_set_File_Dict(Solver_Dict_type *Solver_Dict, File_Dict_type *File_Dict)
    {

    }   
    
    void rCFD_user_set_Phase_Dict(Solver_Dict_type *Solver_Dict, Phase_Dict_type* Phase_Dict)
    {
        int i_phase;
        
        loop_phases_ptr{
            
            if(i_phase == gas){
                
                Phase_Dict[i_phase].number_of_data = 2;

                Phase_Dict[i_phase].time_step = Solver_Dict->global_time_step;
                
                Phase_Dict[i_phase].density = 1.0;
            }

            if(i_phase == solid){
                
                Phase_Dict[i_phase].number_of_data = 3;

                Phase_Dict[i_phase].time_step = Solver_Dict->global_time_step;
                
                Phase_Dict[i_phase].density = 1000.0;
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
                
                if((i_phase == gas) && (i_data == c_gas_A)){
                    
                    Data_Dict[i_phase][i_data].type = concentration_data;
                }

                if((i_phase == gas) && (i_data == c_gas_B)){
                    
                    Data_Dict[i_phase][i_data].type = concentration_data;
                }
                
                if((i_phase == solid) && (i_data == c_drift_0)){
                    
                    Data_Dict[i_phase][i_data].type = concentration_data;
                }

                if((i_phase == solid) && (i_data == c_drift_z_pos)){
                    
                    Data_Dict[i_phase][i_data].type = concentration_data;

                    Data_Dict[i_phase][i_data].drifting_data = 1;
                    
                    Data_Dict[i_phase][i_data].drift_velocity = (double*)malloc( 3 * sizeof(double));
                    
                    Data_Dict[i_phase][i_data].drift_velocity[0] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[1] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[2] = 0.01;
                    
                }

                if((i_phase == solid) && (i_data == c_drift_z_neg)){
                    
                    Data_Dict[i_phase][i_data].type = concentration_data;

                    Data_Dict[i_phase][i_data].drifting_data = 1;
                    
                    Data_Dict[i_phase][i_data].drift_velocity = (double*)malloc( 3 * sizeof(double));
                    
                    Data_Dict[i_phase][i_data].drift_velocity[0] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[1] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[2] = -0.01;                   
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
        
        double  x[3];

        loop_phases_ptr{
        
            loop_cells_ptr{
                
                loop_dim{
                    
                    x[i_dim] = C->x[i_cell][i_dim];
                }
                             
                loop_data{
            
                    if(i_phase == solid){
                        
                        C->data[_i_data] = 1./3.;
#if 0                       
                        if(i_data == c_drift_0){
                            
                            if(x[0] > 0.0){
                                
                                C->data[i_phase][i_cell][i_data] = 1.;
                            }
                            else{
                                
                                C->data[i_phase][i_cell][i_data] = 0.;
                            }
                        }
                        
                        if(i_data == c_drift_z_pos){
                            
                            if(x[1] > 0.0){
                                
                                C->data[i_phase][i_cell][i_data] = 1.;
                            }
                            else{
                                
                                C->data[i_phase][i_cell][i_data] = 0.;
                            }
                        }
                        
                        if(i_data == c_drift_z_neg){
                            
                            if(x[2] < 0.05){
                                
                                C->data[i_phase][i_cell][i_data] = 1.;
                            }
                            else{
                                
                                C->data[i_phase][i_cell][i_data] = 0.;
                            }
                        }
#endif                      
                    }
                    if(i_phase == gas){
                        
                        if(x[1] < -0.074){
                            
                            C->data[i_phase][i_cell][c_gas_A] = 1.0; 
                            C->data[i_phase][i_cell][c_gas_B] = 0.0; 
                        }
                        else{
                            
                            C->data[i_phase][i_cell][c_gas_A] = 0.0; 
                            C->data[i_phase][i_cell][c_gas_B] = 1.0; 
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


        if(i_phase == gas){
            
            int     i_cell, i_frame;
            
            double  mean_value[Phase_Dict[i_phase].number_of_data];

            double  volume_in, volume_in_global, volume_in2, volume_in2_global, volume_out, volume_out_global;
            
            double  primary_inflow, secondary_inflow;
            
            mean_value[c_gas_A] = 0.0; mean_value[c_gas_B] = 0.0;
            
            volume_in = 0.0; volume_in2 = 0.0; volume_out = 0.0;
            
            primary_inflow = 3.28e-3;       /* (kg/s) */
            
            secondary_inflow = 7.14e-5;
            
            loop_cells_ptr{
                
                i_frame = Rec->global_frame[C->island_id[i_cell]];
                
                if(C->x[i_cell][2] < 0.02){
                    
                    C->data[i_phase][i_cell][c_gas_A] = 0.0;
                    
                    C->data[i_phase][i_cell][c_gas_B] = 1.0;
                    
                    volume_in += C->volume[i_cell] * C->vof[i_frame][i_cell][i_phase];
                }

                if(C->x[i_cell][1] < -0.075){
                    
                    C->data[i_phase][i_cell][c_gas_A] = 1.0;
                    
                    C->data[i_phase][i_cell][c_gas_B] = 0.0;
                    
                    volume_in2 += C->volume[i_cell] * C->vof[i_frame][i_cell][i_phase];
                }
                
                if(C->x[i_cell][2] > 0.4){
                    
                    mean_value[c_gas_A] += C->data[i_phase][i_cell][c_gas_A] * C->volume[i_cell] * C->vof[i_frame][i_cell][i_phase];
                    
                    mean_value[c_gas_B] += C->data[i_phase][i_cell][c_gas_B] * C->volume[i_cell] * C->vof[i_frame][i_cell][i_phase];
                    
                    volume_out += C->volume[i_cell] * C->vof[i_frame][i_cell][i_phase];
                }
            }
            
            volume_in_global =  PRF_GRSUM1(volume_in);
            volume_in2_global = PRF_GRSUM1(volume_in2);
            volume_out_global = PRF_GRSUM1(volume_out);
            
            if(volume_in > 0.0){
                
                Balance[gas][c_gas_B].mass_source += primary_inflow * Phase_Dict[gas].time_step * volume_in / volume_in_global;
            }
            
            if(volume_in2 > 0.0){
                
                Balance[gas][c_gas_A].mass_source += secondary_inflow * Phase_Dict[gas].time_step * volume_in2 / volume_in2_global;
            }

            if(volume_out > 0.0){
                
                mean_value[c_gas_A] /= volume_out;
                mean_value[c_gas_B] /= volume_out;
                            
                Balance[gas][c_gas_A].mass_source -= mean_value[c_gas_A] * (primary_inflow + secondary_inflow) * Phase_Dict[gas].time_step * volume_out / volume_out_global;                
                Balance[gas][c_gas_B].mass_source -= mean_value[c_gas_B] * (primary_inflow + secondary_inflow) * Phase_Dict[gas].time_step * volume_out / volume_out_global;
            }
        }

        if(Solver_Dict.verbal){
            
            Message0("\n\nUser before Shift: i_phase %d, i_layer %d", i_phase, i_layer);
        }
#endif
    }

    void rCFD_user_access_data_after_swap(Balance_type** Balance, Phase_Dict_type* Phase_Dict, 
        Data_Dict_type **Data_Dict, Topo_Dict_type *Topo_Dict, Cell_type *C, 
        const short i_phase, const short i_layer)
    {

    }

    void rCFD_user_post(Solver_Dict_type *Solver_Dict, Phase_Dict_type *Phase_Dict, Topo_Dict_type *Topo_Dict, Cell_type *C, Rec_type *Rec)
    {
        Domain  *d=Get_Domain(1);
        Thread  *t;

        int i_phase, i_cell, i_data, i_UDMI, i_frame;
        
        thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){begin_c_loop_int(i_cell, t){
            
            i_UDMI = 0;     /* start index */
            
            i_frame = Rec->global_frame[C->island_id[i_cell]];
            
            C_UDMI(i_cell, t, i_UDMI) = C->vof[i_frame][i_cell][i_phase];
            
            i_UDMI++;

            i_phase = gas;
            
            loop_data{
            
                i_frame = Rec->global_frame[C->island_id[i_cell]];
                
                C_UDMI(i_cell, t, i_UDMI) = C->data[i_phase][i_cell][i_data];

                i_UDMI++;
            }

            i_phase = solid;
            
            loop_data{
            
                i_frame = Rec->global_frame[C->island_id[i_cell]];
                
                C_UDMI(i_cell, t, i_UDMI) = C->data[i_phase][i_cell][i_data] * C->vof[i_frame][i_cell][i_phase];;

                i_UDMI++;
            }
            
        }end_c_loop_int(i_cell, t)}}  

        /* eval segregation */
        if(i_phase == solid){
            
            int     i_layer = 0;
            
            double  data_z_0, data_z_pos, data_z_neg;
            double  vol_z_0, vol_z_pos, vol_z_neg;
            
            data_z_0 = 0.0; data_z_pos = 0.0; data_z_neg = 0.0;         
            vol_z_0 =  0.0; vol_z_pos =  0.0; vol_z_neg =  0.0;
            
            loop_cells_ptr{
                
                i_frame = Rec->global_frame[C->island_id[i_cell]];
                
                i_data = c_drift_0;
                
                data_z_0 += C->data[_i_data] * C->volume[i_cell] * C->vof[_i_vof] * C->x[i_cell][2];
                
                vol_z_0 += C->data[_i_data] * C->volume[i_cell] * C->vof[_i_vof];
                
                i_data = c_drift_z_pos;
                
                data_z_pos += C->data[_i_data] * C->volume[i_cell] * C->vof[_i_vof] * C->x[i_cell][2];
                
                vol_z_pos += C->data[_i_data] * C->volume[i_cell] * C->vof[_i_vof];
                
                i_data = c_drift_z_neg;
                
                data_z_neg += C->data[_i_data] * C->volume[i_cell] * C->vof[_i_vof] * C->x[i_cell][2];
                
                vol_z_neg += C->data[_i_data] * C->volume[i_cell] * C->vof[_i_vof];
            }
            
            data_z_0 =      PRF_GRSUM1(data_z_0);
            data_z_pos =    PRF_GRSUM1(data_z_pos);
            data_z_neg =    PRF_GRSUM1(data_z_neg);
            
            vol_z_0 =       PRF_GRSUM1(vol_z_0);
            vol_z_pos =     PRF_GRSUM1(vol_z_pos);
            vol_z_neg =     PRF_GRSUM1(vol_z_neg);
            
            if(vol_z_0 > 0.0){
                
                data_z_0 /=     vol_z_0;
                
                Message0("\n\nPost: i_phase %d data_z_0 %e ", i_phase, data_z_0);
            }
            
            if(vol_z_pos > 0.0){
                
                data_z_pos /=   vol_z_pos;
                
                Message0("data_z_pos %e ", data_z_pos);
            }
            
            if(vol_z_neg > 0.0){
                
                data_z_neg /=   vol_z_neg;
                
                Message0("data_z_neg %e ", data_z_neg);
            }
        }
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