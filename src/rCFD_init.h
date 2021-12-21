#ifndef RCFD_INIT
#define RCFD_INIT

#include "rCFD_types.h"
#include "rCFD_globals.h"
#include "rCFD_macros.h"
#include "rCFD_parallel.h"
#include "rCFD_user.h"

/* (C)  2021 
    Stefan Pirker
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/  

void init_all(void)
{

    /* D1. Solver_Dict & Solver */
    {
        rCFD_default_Solver_Dict(&Solver_Dict);
        
        rCFD_user_set_Solver_Dict(&Solver_Dict);
        
#if RP_NODE     
        rCFD_default_Solver(&Solver);
#endif      
    }

    /* D2. File_Dict */
    {
        rCFD_default_File_Dict(&Solver_Dict, &File_Dict);
        
        rCFD_user_set_File_Dict(&Solver_Dict, &File_Dict);
    }

    /* D3. Phase_Dict */
    {
#if RP_NODE     
        Phase_Dict = (Phase_Dict_type*)malloc(Solver_Dict.number_of_phases * sizeof(Phase_Dict_type));
        
        rCFD_default_Phase_Dict(&Solver_Dict, Phase_Dict);
        
        rCFD_user_set_Phase_Dict(&Solver_Dict, Phase_Dict); 
#endif      
    }       
    
    /* D4. Tracer_Dict */
    {
#if RP_NODE
        Tracer_Dict.random_walk = (short*)malloc(Solver_Dict.number_of_phases * sizeof(short));
        
        rCFD_default_Tracer_Dict(&Solver_Dict, &Tracer_Dict);

        rCFD_user_set_Tracer_Dict(&Solver_Dict, &Tracer_Dict);
#endif      
    }
    
    /* D5. Norm_Dict */
    {
#if RP_NODE
        rCFD_default_Norm_Dict(&Solver_Dict, &Norm_Dict);

        rCFD_user_set_Norm_Dict(&Solver_Dict, &Norm_Dict);
#endif      
    }   

    /* D6. Rec_Dict */
    {
        rCFD_default_Rec_Dict(&Solver_Dict, &Rec_Dict);

        rCFD_user_set_Rec_Dict(&Solver_Dict, &Rec_Dict);    
    }   

    /* D7. Data_Dict */ 
    {
#if RP_NODE     
        int i_phase;
        
        Data_Dict = (Data_Dict_type**)malloc(Solver_Dict.number_of_phases * sizeof(Data_Dict_type*));
        
        loop_phases{
            
            Data_Dict[i_phase] = (Data_Dict_type*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(Data_Dict_type));
        }   
        
        rCFD_default_Data_Dict(&Solver_Dict, Phase_Dict, Data_Dict);
        
        rCFD_user_set_Data_Dict(&Solver_Dict, Phase_Dict, Data_Dict);
#endif      
    }
    
    /* D8. Balance_Dict */
    {
#if RP_NODE
        int i_phase;
        
        Balance_Dict = (Balance_Dict_type**)malloc(Solver_Dict.number_of_phases * sizeof(Balance_Dict_type*));
        
        loop_phases{
            
            Balance_Dict[i_phase] = (Balance_Dict_type*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(Balance_Dict_type));
        }
        
        rCFD_default_Balance_Dict(&Solver_Dict, Phase_Dict, Balance_Dict);
        
        rCFD_user_set_Balance_Dict(&Solver_Dict, Phase_Dict, Balance_Dict);
#endif      
    }

    /* D9. Topo_Dict */
    {
        rCFD_default_Topo_Dict(&Solver_Dict, &Topo_Dict);
        
        rCFD_user_set_Topo_Dict(&Solver_Dict, &Topo_Dict);
        
#if RP_NODE     
        int i_layer;
        
        Topo_Dict.Cell_Dict = (Cell_Dict_type*)malloc(Topo_Dict.number_of_layers * sizeof(Cell_Dict_type));
        
        Topo_Dict.Face_Dict = (Face_Dict_type*)malloc(Topo_Dict.number_of_layers * sizeof(Face_Dict_type));

        loop_layers{
                
            rCFD_default_Cell_Dict(&Solver_Dict, &Topo_Dict.Cell_Dict[i_layer], i_layer);

            rCFD_user_set_Cell_Dict(&Solver_Dict, &Topo_Dict.Cell_Dict[i_layer], i_layer);

            rCFD_default_Face_Dict(&Solver_Dict, &Topo_Dict.Face_Dict[i_layer], i_layer);

            rCFD_user_set_Face_Dict(&Solver_Dict, &Topo_Dict.Face_Dict[i_layer], i_layer);          
        }   
#endif  
    
    } 
    
        
    /* G1. Cells (first initialization) */
    {
#if RP_NODE     
        int     i_phase, i_frame, i_cell, i_layer;
        
        loop_layers{
            
            if(i_layer == 0){
            
                C.x = (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
                
                loop_cells{
                    
                    C.x[i_cell] = (double*)malloc( 3 * sizeof(double));
                }
                
                C.volume = (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));
                
                C.average_velocity =    (double**)malloc(Solver_Dict.number_of_phases * sizeof(double*));
                C.crossing_time =       (double**)malloc(Solver_Dict.number_of_phases * sizeof(double*));
                
                loop_phases{
                    
                    C.average_velocity[i_phase] =   (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));            
                    C.crossing_time[i_phase] =      (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));
                }

                C.hit_by_other_cell =   (short*)malloc(_Cell_Dict.number_of_cells * sizeof(short));

                C.island_id =           (short*)malloc(_Cell_Dict.number_of_cells * sizeof(short));
                
                C.weight_after_shift =  (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));
                C.weight_after_swap =   (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));

                C.vof = (double***)malloc(Solver_Dict.number_of_frames * sizeof(double**));
                
                loop_frames{
                    
                    C.vof[i_frame] = (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
                    
                    loop_cells{
                        
                        C.vof[i_frame][i_cell] = (double*)malloc(Solver_Dict.number_of_phases * sizeof(double));
                    }
                }

                C.data =        (double***)malloc(Solver_Dict.number_of_phases * sizeof(double**));
                C.data_shift =  (double***)malloc(Solver_Dict.number_of_phases * sizeof(double**));
                C.data_swap =   (double***)malloc(Solver_Dict.number_of_phases * sizeof(double**));
                
                loop_phases{
                    
                    C.data[i_phase] =       (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
                    C.data_shift[i_phase] = (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
                    C.data_swap[i_phase] =  (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
                    
                    loop_cells{
                        
                        C.data[i_phase][i_cell] =       (double*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(double));
                        C.data_shift[i_phase][i_cell] = (double*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(double));
                        C.data_swap[i_phase][i_cell] =  (double*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(double));
                    }
                }       
                
                if(Solver_Dict.data_drifting_on){
                    
                    C.drift_exchange = (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));
                }
                else{
                    
                    C.drift_exchange = NULL;
                }       
                
                if(_Cell_Dict.number_of_user_vars > 0){
                    
                    C.user = (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
                    
                    loop_cells{
                        
                        C.user[i_cell] = (double*)malloc(_Cell_Dict.number_of_user_vars * sizeof(double));
                    }
                }
                else{
                    
                    C.user = NULL;
                }
                        
                rCFD_default_Cell(&Solver_Dict, &File_Dict, Phase_Dict, &Topo_Dict, &C, i_layer);
            }
        }
#endif  
    }   

    /* G2. Faces */
    {
#if RP_NODE
        int i_face, i_layer;
        
        loop_layers{
            
            if(i_layer == 0){
                
                F.c0 = (int*)malloc(_Face_Dict.number_of_faces * sizeof(int));
        
                F.c1 = (int*)malloc(_Face_Dict.number_of_faces * sizeof(int));
                
                F.area = (double**)malloc(_Face_Dict.number_of_faces * sizeof(double*));
                
                loop_faces{
                    
                    F.area[i_face] = (double*)malloc( 3 * sizeof(double));
                }

                rCFD_default_Face(&Solver_Dict, &Topo_Dict, &F, i_layer);
            }
        }
#endif      
    }
    
    /* G3. Tracer */
    {
#if RP_NODE

        Tracer.monitor_counter = (int*)malloc(Solver_Dict.number_of_phases * sizeof(int));

        Tracer.number_of_shifts = (int*)malloc(Solver_Dict.number_of_phases * sizeof(int));
        
        Tracer.shifts = (C2C_shift_type**)malloc(Solver_Dict.number_of_phases * sizeof(C2C_shift_type*));
        
        rCFD_default_Tracer(&Solver_Dict, &Tracer);

#endif      
    }

    /* G4. Norms */
    {
#if RP_NODE
        Domain  *d=Get_Domain(1);
        Thread  *t;
        cell_t  i_cell;
        
        int     number_of_norms = 0;
        
        thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){begin_c_loop_int(i_cell,t){

            if((i_cell % Norm_Dict.coarse_graining) == 0){
                
                number_of_norms++;              
            }
            
        }end_c_loop_int(i_cell,t)}}
        
        Norms.number_of_norms = number_of_norms;
        
        Norms.norm = (double*)malloc(number_of_norms * sizeof(double));
        
        rCFD_default_Norms(&Solver_Dict, &Norms);
        
#endif      
    }

    /* G5. Rec */
    {
        int     i_state, i_state2, i_island;
        
        Rec.global_frame = (int*)malloc(Solver_Dict.number_of_islands * sizeof(int));
        
        loop_islands{
            
            Rec.global_frame[i_island] = 0;
        }
        
        Rec.jumps = (int****)malloc(Solver_Dict.number_of_states * sizeof(int***));
        
        loop_states{
            
            Rec.jumps[i_state] = (int***)malloc(Solver_Dict.number_of_states * sizeof(int**));
            
            loop_states2{
                
                Rec.jumps[i_state][i_state2] = (int**)malloc(Solver_Dict.number_of_islands * sizeof(int*));
                
                loop_islands{
                    
                    Rec.jumps[i_state][i_state2][i_island] = (int*)malloc(Solver_Dict.number_of_frames * sizeof(int));
                }
            }
        }
        
        rCFD_default_Rec(&Solver_Dict, &File_Dict, &Rec_Dict, &Rec);
    }   

    /* G6. Balance (first initialization) */    
    {
#if RP_NODE
        int i_phase, i_data, i_node, i_node2;
        
        Balance = (Balance_type**)malloc(Solver_Dict.number_of_phases * sizeof(Balance_type*));
        
        loop_phases{
            
            Balance[i_phase] = (Balance_type*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(Balance_type));
            
            loop_data{ 
            
                /* global_balancing */
                if(Balance_Dict[i_phase][i_data].type == global_balancing){
                                        
                    Balance[i_phase][i_data].node2node_flux = NULL;
                    Balance[i_phase][i_data].node2node_data_flux = NULL;
                }
                
                /* per node balancing */
                else if(Balance_Dict[i_phase][i_data].type == per_node_balancing){
                    
                    Balance[i_phase][i_data].node2node_flux = (double**)malloc((node_last + 1) * sizeof(double*));
                    Balance[i_phase][i_data].node2node_data_flux = (double**)malloc((node_last + 1) * sizeof(double*));
                    
                    for(i_node = 0; i_node < (node_last + 1); i_node++){
                        
                        Balance[i_phase][i_data].node2node_flux[i_node] = (double*)malloc((node_last + 1) * sizeof(double));
                        Balance[i_phase][i_data].node2node_data_flux[i_node] = (double*)malloc((node_last + 1) * sizeof(double));
                        
                        for(i_node2 = 0; i_node2 < (node_last + 1); i_node2++){
                            
                            Balance[i_phase][i_data].node2node_flux[i_node][i_node2] = 0.0;
                            Balance[i_phase][i_data].node2node_data_flux[i_node][i_node2] = 0.0;
                        }
                    }
                }
                
                Balance[i_phase][i_data].mass_integral = 0.0;
                
                Balance[i_phase][i_data].mass_integral_target = Balance[i_phase][i_data].mass_integral;
                
                Balance[i_phase][i_data].mass_integral_global = PRF_GRSUM1(Balance[i_phase][i_data].mass_integral);

                Balance[i_phase][i_data].mass_integral_target_global = Balance[i_phase][i_data].mass_integral_global;
                
                Balance[i_phase][i_data].mass_source = 0.0;

                Balance[i_phase][i_data].mass_source_global = 0.0;
                
                Balance[i_phase][i_data].mass_error = 0.0;

                Balance[i_phase][i_data].mass_error_global = 0.0;

                Balance[i_phase][i_data].mass_error_prev = 0.0;

                Balance[i_phase][i_data].face_swap = Solver_Dict.face_swap_min;

                Balance[i_phase][i_data].face_swap_loops = 1;               
            }           
        }
#endif      
    }

    /* G1 + G6: Init data fields and update balances */
    {
#if RP_NODE     
        int i_layer, i_frame, i_phase, i_data, i_cell;
        
        i_layer = 0;    /* user initialization just for base layer of grid, might be changed to loop_layers in future */
                
        rCFD_user_init_Data(&Solver_Dict, Balance, Phase_Dict, Data_Dict, &Topo_Dict, &C, i_layer);
        
        i_frame = 0;
        
        loop_int_cells{
            
            loop_phases{
            
                loop_data{
                    
                    switch (Data_Dict[i_phase][i_data].type){
                        
                        case generic_data:
                            
                            Balance[_i_balance].mass_integral += C.data[_i_data];
                            
                            break;
                            
                        case concentration_data:
                        
                            Balance[_i_balance].mass_integral += C.data[_i_data] * 
                            
                                Phase_Dict[i_phase].density * C.volume[i_cell] * C.vof[_i_vof];
                            
                            break;
                            
                        case temperature_data:
                        
                            Balance[_i_balance].mass_integral += C.data[_i_data] * 
                    
                                Phase_Dict[i_phase].heat_capacity * Phase_Dict[i_phase].density * C.volume[i_cell] * C.vof[_i_vof];
                                
                            break;
                    
                        default: break;
                    }
                }
            }       
        }

        loop_phases{
            
            loop_data{
                            
                Balance[i_phase][i_data].mass_integral_target = Balance[i_phase][i_data].mass_integral;
                
                Balance[i_phase][i_data].mass_integral_global = PRF_GRSUM1(Balance[i_phase][i_data].mass_integral);

                Balance[i_phase][i_data].mass_integral_target_global = Balance[i_phase][i_data].mass_integral_global;               
            }
        }
        
        rCFD_user_post(&Solver_Dict, Phase_Dict, &Topo_Dict, &C, &Rec);   /* post-process initialization */
#endif
    }       

    /* P: Parallel grid communication */
    {
#if RP_NODE
        /* for the moment, just consider MPI communication for i_layer = 0 */
        int i_layer = 0;
        
        init_parallel_grid(&Solver_Dict, &Topo_Dict, i_layer);
#endif
    }       
}

#endif