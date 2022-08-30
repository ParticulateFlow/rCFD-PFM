#ifndef RCFD_USER
#define RCFD_USER

#include "rCFD_types.h"
#include "rCFD_globals.h"
#include "rCFD_defaults.h"
#include "rCFD_macros.h"
#include "rCFD_layer.h"

/* (C)  2021 
    Stefan Pirker
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/  


    /* user Names */
    /*************************************************************************************/ 

#if 1

#if RP_NODE
        static short    ref_species_initiated = 0;

        static double   mean_c_drift_0_conc = 0.0;
#endif

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
    
    void rCFD_user_set_Solver_Dict(void)
    {
        Solver_Dict.number_of_phases =                     2;
        
        Solver_Dict.number_of_layers =                     1;

        Solver_Dict.max_number_of_cells_per_time_step =    30;

        Solver_Dict.global_time_step =                     0.005;   
        
        Solver_Dict.time_steps_per_monitoring_interval =   10;
        
        Solver_Dict.number_of_frames =                     50;         

        Solver_Dict.number_of_runs =                       1; 
        
        Solver_Dict.data_drifting_on =                     1;

        Solver_Dict.face_diffusion_on =                    0; 
        
        Solver_Dict.recurrence_process_on =                1;  

        Solver_Dict.face_swap_max_per_loop =               (1./8.); 

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
            
            if(i_phase == gas){
                
                Phase_Dict[i_phase].number_of_data = 2;

                Phase_Dict[i_phase].time_step = Solver_Dict.global_time_step;
                
                Phase_Dict[i_phase].density = 1.0;
            }

            if(i_phase == solid){
                
                Phase_Dict[i_phase].number_of_data = 3;

                Phase_Dict[i_phase].time_step = Solver_Dict.global_time_step;
                
                Phase_Dict[i_phase].density = 1000.0;
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
        
        int     i_phase, i_cell, i_data, i_dim, i_frame;

        i_data = c_drift_0;
        
        double  x[3], mixing_nom, mixing_denom;

        loop_phases{
        
            loop_cells{
                    
                loop_dim{
                    
                    x[i_dim] = _C.x[i_cell][i_dim];
                }
                             
                loop_data{
            
                    if(i_phase == solid){
                        
                        _C.data[_i_data] = 1./3.;
                        
                        if(i_data == c_drift_0){
                            
                            if(x[1] > 0.0){
                                
                                _C.data[_i_data] = 1.0;                            /* mixing condition for solid phase */
                            }
                            else{
                                
                                _C.data[_i_data] = 0.0;
                            }
                        }       
                    }
                    
                    if(i_phase == gas){
                        
                        if(x[1] < -0.074){
                            
                            _C.data[i_phase][i_cell][c_gas_A] = 1.0;              /* mixing condition for gas phase */
                            _C.data[i_phase][i_cell][c_gas_B] = 0.0; 
                        }
                        else{
                            
                            _C.data[i_phase][i_cell][c_gas_A] = 0.0; 
                            _C.data[i_phase][i_cell][c_gas_B] = 1.0; 
                        }
                    }
                }
            }
        }
        
        /* initialize mean_solid_A_conc (using +=) */   
        {       
            mixing_nom = 0.0;
            
            mixing_denom = 0.0;
                
            i_phase = solid;
            
            i_data = c_drift_0;
            
            loop_cells{
                    
                i_frame = Rec.global_frame[_C.island_id[i_cell]];
                
                mixing_nom += _C.data[_i_data] * _C.volume[i_cell] * _C.vof[i_frame][i_cell][i_phase];
                
                mixing_denom += _C.volume[i_cell] * _C.vof[i_frame][i_cell][i_phase];
            }
        
            mixing_nom = PRF_GRSUM1(mixing_nom);            /* sum of mixing_nom */
            
            mixing_denom = PRF_GRSUM1(mixing_denom);        /* sum of mixing_denom */
            
            if(mixing_denom > 0.0){
                
                mean_c_drift_0_conc = mixing_nom / mixing_denom;              /* mean solid concentration */
                
            }
            else{
                
                mean_c_drift_0_conc = 0.0;          
            }
        }

        Message0("\n\n...rCFD_user_init_Data\n");               
#endif
    }

#endif
    
    /* user accesses Data */
    /*************************************************************************************/

#if 1

    short rCFD_user_set_layer(const short i_layer)
    {
        return 0;
    }
    
    short rCFD_user_phase_switch(const short i_phase)
    {
#if RP_NODE

        short   i_phase_new = i_phase;
        
        if(i_phase == ph_gas){
            
            i_phase_new = ph_solid;
        }
        
        if((i_phase_new < 0) || (i_phase_new >= Solver_Dict.number_of_phases)){
            
            i_phase_new = i_phase;
        }
        
        return i_phase_new;
#else

        return 0;   /* just to avoid compiler complains */
#endif
    }       

    void rCFD_user_access_data_before_shift(const short i_phase, const short i_layer)
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
            
            loop_cells{
                
                i_frame = Rec.global_frame[_C.island_id[i_cell]];
                
                if(_C.x[i_cell][2] < 0.02){
                    
                    _C.data[i_phase][i_cell][c_gas_A] = 0.0;
                    
                    _C.data[i_phase][i_cell][c_gas_B] = 1.0;
                    
                    volume_in += _C.volume[i_cell] * _C.vof[i_frame][i_cell][i_phase];
                }

                if(_C.x[i_cell][1] < -0.075){
                    
                    _C.data[i_phase][i_cell][c_gas_A] = 1.0;
                    
                    _C.data[i_phase][i_cell][c_gas_B] = 0.0;
                    
                    volume_in2 += _C.volume[i_cell] * _C.vof[i_frame][i_cell][i_phase];
                }
                
                if(_C.x[i_cell][2] > 0.4){
                    
                    mean_value[c_gas_A] += _C.data[i_phase][i_cell][c_gas_A] * _C.volume[i_cell] * _C.vof[i_frame][i_cell][i_phase];
                    
                    mean_value[c_gas_B] += _C.data[i_phase][i_cell][c_gas_B] * _C.volume[i_cell] * _C.vof[i_frame][i_cell][i_phase];
                    
                    volume_out += _C.volume[i_cell] * _C.vof[i_frame][i_cell][i_phase];
                }
            }
            
            volume_in_global =  PRF_GRSUM1(volume_in);
            volume_in2_global = PRF_GRSUM1(volume_in2);
            volume_out_global = PRF_GRSUM1(volume_out);
            
            if(volume_in > 0.0){
                
                Balance[gas][c_gas_B].mass_source += primary_inflow * Solver.timestep_width[i_layer] * volume_in / volume_in_global;
            }
            
            if(volume_in2 > 0.0){
                
                Balance[gas][c_gas_A].mass_source += secondary_inflow * Solver.timestep_width[i_layer] * volume_in2 / volume_in2_global;
            }

            if(volume_out > 0.0){
                
                mean_value[c_gas_A] /= volume_out;
                mean_value[c_gas_B] /= volume_out;
                            
                Balance[gas][c_gas_A].mass_source -= mean_value[c_gas_A] * (primary_inflow + secondary_inflow) * Solver.timestep_width[i_layer] * volume_out / volume_out_global;                
                Balance[gas][c_gas_B].mass_source -= mean_value[c_gas_B] * (primary_inflow + secondary_inflow) * Solver.timestep_width[i_layer] * volume_out / volume_out_global;
            }
        }

        if(Solver_Dict.verbose){
            
            Message0("\n\nUser before Shift: i_phase %d, i_layer %d", i_phase, i_layer);
        }
#endif
    }

    void rCFD_user_access_data_after_swap(const short i_phase, const short i_layer)
    {

    }

    void rCFD_user_post()
    {
#if RP_NODE

        int i_layer, i_phase, i_cell, i_data, i_UDMI, i_frame;

        Domain  *d=Get_Domain(1);
        Thread  *t;
        
        double mass_c_gas_A, mixing_index, mixing_nom, mixing_denom;
        
        i_layer = 0;
        
        if(Solver.current_layer > 0){
            
            rCFD_map_from_to_layer(Solver.current_layer, i_layer);
        }      
        
        /* P.1. eval mass_c_gas_A and mixing_index */
        {   
            mass_c_gas_A = 0.0;
            
            i_phase = gas;
            
            i_data = c_gas_A;
            
            loop_int_cells{
                        
                i_frame = Rec.global_frame[_C.island_id[i_cell]];
                
                mass_c_gas_A +=  _C.data[_i_data] * _C.volume[i_cell] * _C.vof[_i_vof];                                              
            }
            
            mass_c_gas_A = PRF_GRSUM1(mass_c_gas_A);                                                                                                                                                                                    /* sum of mass of gas B */
            
            mixing_nom = 0.0;
            
            mixing_denom = 0.0;
            
            i_phase = solid;
            
            thread_loop_c(t, d){
                
                if(FLUID_CELL_THREAD_P(t)){
                    
                    begin_c_loop_int(i_cell,t){
                        
                        i_frame = Rec.global_frame[_C.island_id[i_cell]];
                        
                        mixing_nom += fabs( _C.data[i_phase][i_cell][c_drift_0] - mean_c_drift_0_conc) * _C.volume[i_cell] * _C.vof[i_frame][i_cell][i_phase];                                 
                        
                        /* species mass fraction - mean solid concentration --------- that should tend to value 0 */
                        
                        mixing_denom += _C.volume[i_cell]  * _C.vof[i_frame][i_cell][i_phase];
                        
                    }end_c_loop_int(i_cell,t)
                }
            }
            
            mixing_nom = PRF_GRSUM1(mixing_nom);                                                                                                                                                                                                  /* sum of species mixing */
            
            mixing_denom = PRF_GRSUM1(mixing_denom);
            
            if(mixing_denom > 0.0){
                
                mixing_index = 1.0 - 2.0 * mixing_nom / mixing_denom;
                
                /* mixing index should be within the value of 0 and 1 */
            }
            else{
                
                mixing_index = 0.0;
            }
        }       
               
        /* P.2. node-0 writes values to ref file */
        if(myid == 0){

            FILE    *f_out = NULL;
            
            if(ref_species_initiated == 0){
                
                f_out = fopen("./monitor_rCFD.out","w");
            }
            else{
                f_out = fopen("./monitor_rCFD.out","a");
            }
            
            if(f_out == NULL){
                
                Message0("\nERROR: Could not open ref monitor file\n");
                
                return;
            }
            
            fprintf(f_out, "%e %e %e\n", Solver.global_time, mass_c_gas_A, mixing_index);
            
            fclose(f_out);
        }
        
        /* P.3. nodes write species conc to UDMI */
        {

            thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){begin_c_loop_int(i_cell, t){
                             
                i_frame = Rec.global_frame[_C.island_id[i_cell]];
                     
                i_UDMI = 1;     /* start index */
                     
                i_phase = gas;
                     
                i_data = c_gas_A;
                            
                C_UDMI(i_cell, t, i_UDMI) = _C.data[_i_data];
                
                i_UDMI = 2;
                        
                i_phase = solid;
                     
                i_data = c_drift_0;
                
                C_UDMI(i_cell, t, i_UDMI) = _C.data[_i_data] * _C.vof[_i_vof];
                                    
            }end_c_loop_int(i_cell, t)}}
        }
        
        ref_species_initiated = 1;
        
        Message0("\n\n...CFD_ref_monitors\n");
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
