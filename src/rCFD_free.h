#ifndef RCFD_FREE
#define RCFD_FREE

#include "rCFD_types.h"
#include "rCFD_globals.h"
#include "rCFD_macros.h"

/* (C)  2021 
    Stefan Pirker
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/  

    void free_all(void)
    {       
        int     i_state, i_state2, i_island;
        
        if(Rec.global_frame != NULL){
            
            free(Rec.global_frame);
        }
        
        if(Rec.jumps != NULL){
            
            loop_states{
            
                if(Rec.jumps[i_state] != NULL){
                    
                    loop_states2{
                    
                        if(Rec.jumps[i_state][i_state2] != NULL){
                            
                            loop_islands{
                                
                                if(Rec.jumps[i_state][i_state2][i_island] != NULL){
                                    
                                    free(Rec.jumps[i_state][i_state2][i_island]);
                                }
                            }
                            
                            free(Rec.jumps[i_state][i_state2]);
                        }
                    }
                    
                    free(Rec.jumps[i_state]);
                }
            }
            
            free(Rec.jumps);
        }
        
#if RP_NODE
        int     i_phase, i_frame, i_data, i_shift, i_node, i_cell, i_face, i_layer;
        
        int     number_of_shifts;
        
        if(Tracer_Dict.random_walk != NULL){ 
        
            free(Tracer_Dict.random_walk);  
        }
        
        i_layer = 0;    /* for the moment, we only consider the base grid */

        if(C.x != NULL){
            
            loop_cells{
                
                if(C.x[i_cell] != NULL){
                    
                    free(C.x[i_cell]);
                }
            }
                
            free(C.x);
        }
        
        if(C.volume != NULL){
            
            free(C.volume);
        }
        
        if(C.average_velocity != NULL){
            
            loop_phases{
                
                if(C.average_velocity[i_phase] != NULL){
                    
                    free(C.average_velocity[i_phase]);
                }
            }
            
            free(C.average_velocity);
        }   

        if(C.crossing_time != NULL){
            
            loop_phases{
                
                if(C.crossing_time[i_phase] != NULL){
                    
                    free(C.crossing_time[i_phase]);
                }
            }
            
            free(C.crossing_time);
        }   

        if(C.hit_by_other_cell != NULL){ 
        
            free(C.hit_by_other_cell);  
        }

        if(C.island_id != NULL){ 
        
            free(C.island_id);  
        }
        
        if(C.weight_after_shift != NULL){
            
            free(C.weight_after_shift);
        }
        
        if(C.weight_after_swap != NULL){
            
            free(C.weight_after_swap);
        }
        
        if(C.vof != NULL){
            
            loop_frames{
                
                if(C.vof[i_frame] != NULL){
                    
                    loop_cells{
                        
                        if(C.vof[i_frame][i_cell] != NULL){
                            
                            free(C.vof[i_frame][i_cell]);
                        }
                    }
                    
                    free(C.vof[i_frame]);
                }
            }
            
            free(C.vof);
        }
        
        if(C.data != NULL){
            
            loop_phases{
                
                if(C.data[i_phase] != NULL){
                
                    loop_cells{
                        
                        if(C.data[i_phase][i_cell] != NULL){
                            
                            free(C.data[i_phase][i_cell]);
                        }
                    }
                    
                    free(C.data[i_phase]);
                }
            }
            
            free(C.data);
        }

        if(C.data_shift != NULL){
            
            loop_phases{
                
                if(C.data_shift[i_phase] != NULL){
                
                    loop_cells{
                        
                        if(C.data_shift[i_phase][i_cell] != NULL){
                            
                            free(C.data_shift[i_phase][i_cell]);
                        }
                    }
                    
                    free(C.data_shift[i_phase]);
                }
            }
            
            free(C.data_shift);
        }
        
        if(C.data_swap != NULL){
            
            loop_phases{
                
                if(C.data_swap[i_phase] != NULL){
                
                    loop_cells{
                        
                        if(C.data_swap[i_phase][i_cell] != NULL){
                            
                            free(C.data_swap[i_phase][i_cell]);
                        }
                    }
                    
                    free(C.data_swap[i_phase]);
                }
            }
            
            free(C.data_swap);
        }

        if(C.drift_exchange != NULL){
            
            free(C.drift_exchange);
        }
        
        if(C.user != NULL){
            
            loop_cells{
                
                if(C.user[i_cell] !=  NULL){
                    
                    free(C.user[i_cell]);
                }
            }
            
            free(C.user);
        }

        if(F.c0 != NULL){
            
            free(F.c0);
        }
        
        if(F.c1 != NULL){
            
            free(F.c1);
        }
        
        if(F.area != NULL){
            
            loop_faces{
                
                if(F.area[i_face] != NULL){
                    
                    free(F.area[i_face]);
                }
            }
            
            free(F.area);
        }

        /* free Topo_Dict after C and F */
        
        if(Topo_Dict.Cell_Dict != NULL){
            
            free(Topo_Dict.Cell_Dict);
        }
        
        if(Topo_Dict.Face_Dict != NULL){
            
            free(Topo_Dict.Face_Dict);
        }
                        
        if(Tracer.monitor_counter != NULL){ 
        
            free(Tracer.monitor_counter);   
        }

        if(Tracer.number_of_shifts != NULL){ 
        
            free(Tracer.number_of_shifts);  
        }

        if(Tracer.shifts != NULL){

            loop_phases{
                
                if(Tracer.shifts[i_phase] != NULL){
                    
                    free(Tracer.shifts[i_phase]);
                }
            }
        
            free(Tracer.shifts);    
        }   

        if(Norms.norm != NULL){
            
            free(Norms.norm);
        }
        
        if(C2Cs != NULL){
            
            loop_states{
                
                if(C2Cs[i_state] != NULL){
                    
                    loop_phases{
                        
                        if(C2Cs[i_state][i_phase] != NULL){
                            
                            loop_frames{
                                    
                                if(C2Cs[_i_C2C].shifts != NULL){
                                    
                                    free(C2Cs[_i_C2C].shifts);
                                }
                                
                                if(C2Cs[_i_C2C].shifts_in != NULL){
                                    
                                    free(C2Cs[_i_C2C].shifts_in);
                                }
                                
                                if(C2Cs[_i_C2C].shifts_out != NULL){
                                    
                                    free(C2Cs[_i_C2C].shifts_out);
                                }
                                
                                if(C2Cs[_i_C2C].island_offsets != NULL){
                                    
                                    free(C2Cs[_i_C2C].island_offsets);
                                }

                                if(C2Cs[_i_C2C].island_offsets_in != NULL){
                                    
                                    free(C2Cs[_i_C2C].island_offsets_in);
                                }
                                
                                if(myid == 0){
                                    
                                    if(C2Cs[_i_C2C].number_of_shifts_to_node_zero != NULL){
                                        
                                        free(C2Cs[_i_C2C].number_of_shifts_to_node_zero);
                                    }

                                    if(C2Cs[_i_C2C].number_of_shifts_from_node_zero != NULL){
                                        
                                        free(C2Cs[_i_C2C].number_of_shifts_from_node_zero);
                                    }
                                    
                                    if(C2Cs[_i_C2C].in2out != NULL){
                                        
                                        for(i_node = 0; i_node < (node_last + 1); i_node++){
                                            
                                            if(C2Cs[_i_C2C].in2out[i_node] != NULL){
                                                
                                                free(C2Cs[_i_C2C].in2out[i_node]);
                                            }
                                        }   
                                        
                                        free(C2Cs[_i_C2C].in2out);
                                    }
                                }
                            }
                                
                            free(C2Cs[i_state][i_phase]);                   
                        }
                    }
                    
                    free(C2Cs[i_state]);
                }
            }
            
            free(C2Cs);
        }

        if(C2Cs_MPI.shifts_in != NULL){
            
            number_of_shifts = C2Cs_MPI.max_number_of_MPI_shifts;
            
            loop_shifts{
                
                free(C2Cs_MPI.shifts_in[_i_shift].data);
            }
            
            free(C2Cs_MPI.shifts_in);
        }

        if(Balance != NULL){
            
            loop_phases{
                
                if(Balance[i_phase] != NULL){
                
                    loop_data{
                                            
                        if(Balance[i_phase][i_data].node2node_flux != NULL){
                        
                            for(i_node = 0; i_node < (node_last + 1); i_node++){
                            
                                if(Balance[i_phase][i_data].node2node_flux[i_node] != NULL){
                                    
                                    free(Balance[i_phase][i_data].node2node_flux[i_node]);
                                }
                            }
                            
                            free(Balance[i_phase][i_data].node2node_flux);
                        }

                        if(Balance[i_phase][i_data].node2node_data_flux != NULL){
                        
                            for(i_node = 0; i_node < (node_last + 1); i_node++){
                            
                                if(Balance[i_phase][i_data].node2node_data_flux[i_node] != NULL){
                                    
                                    free(Balance[i_phase][i_data].node2node_data_flux[i_node]);
                                }
                            }
                            
                            free(Balance[i_phase][i_data].node2node_data_flux);
                        }
                    }
                
                    free(Balance[i_phase]);
                }
            }

            free(Balance);
        }

        if(Balance_Dict != NULL){
            
            loop_phases{
                
                if(Balance_Dict[i_phase] != NULL){
                    
                    free(Balance_Dict[i_phase]);
                }
            }
            
            free(Balance_Dict);
        }
        
        if(Data_Dict != NULL){
        
            loop_phases{
                
                if(Data_Dict[i_phase] != NULL){
                    
                    loop_data{
                        
                        if(Data_Dict[i_phase][i_data].drift_velocity != NULL){
                            
                            free(Data_Dict[i_phase][i_data].drift_velocity);
                        }
                        
                        if(Data_Dict[i_phase][i_data].user != NULL){
                            
                            free(Data_Dict[i_phase][i_data].user);
                        }
                    }
                    
                    free(Data_Dict[i_phase]);
                }
            }
            
            free(Data_Dict);
        }
        
        /* free last, because previously loop_data was needed */
        if(Phase_Dict != NULL){
            
            loop_phases{
                
                if(Phase_Dict[i_phase].user != NULL){
                    
                    free(Phase_Dict[i_phase].user);
                }
            }
            
            free(Phase_Dict);
        }
#endif      
    }
#endif