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
        int i_state;

        if(Rec.global_frame != NULL){

            free(Rec.global_frame);
        }

        free_i_4d(Rec.jumps);

        if(Solver.timestep_width != NULL){

            free(Solver.timestep_width);
        }

#if RP_NODE
        int     i_phase, i_frame, i_data, i_shift, i_node, i_cell, i_layer;

        int     number_of_shifts;

        if(Tracer_Dict.random_walk != NULL){

            free(Tracer_Dict.random_walk);
        }

        if(Topo.Cell != NULL){

            loop_layers{

                free_r_2d(_C.x);

                if(_C.volume != NULL){

                    free(_C.volume);
                }

                free_r_2d(_C.average_velocity);

                free_r_2d(_C.crossing_time);

                if(_C.hit_by_other_cell != NULL){

                    free(_C.hit_by_other_cell);
                }

                if(_C.island_id != NULL){

                    free(_C.island_id);
                }

                if(_C.weight_after_shift != NULL){

                    free(_C.weight_after_shift);
                }

                if(_C.weight_after_swap != NULL){

                    free(_C.weight_after_swap);
                }

                free_r_3d(_C.vof);

                free_r_3d(_C.data);

                free_r_3d(_C.data_shift);

                free_r_3d(_C.data_swap);

                if(_C.drift_exchange != NULL){

                    free(_C.drift_exchange);
                }

                free_r_2d(_C.user);

                if(_C.marked != NULL){

                    free(_C.marked);
                }

                if(_C.parent_cell != NULL){

                    free(_C.parent_cell);
                }

                if(_C.number_of_children != NULL){

                    free(_C.number_of_children);
                }

                if(_C.child_index != NULL){

                    loop_cells{

                        if(_C.child_index[i_cell] !=  NULL){

                            free(_C.child_index[i_cell]);
                        }
                    }

                    free(_C.child_index);
                }
            }

            free(Topo.Cell);
        }

        if(Topo.Face != NULL){

            loop_layers{

                if(_F.c0 != NULL){

                    free(_F.c0);
                }

                if(_F.c1 != NULL){

                    free(_F.c1);
                }

                free_r_2d(_F.area);
            }

            free(Topo.Face);
        }

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

                                if(C2Cs[_i_C2C].number_of_shifts != NULL){

                                    free(C2Cs[_i_C2C].number_of_shifts);
                                }

                                if(C2Cs[_i_C2C].number_of_shifts_in != NULL){

                                    free(C2Cs[_i_C2C].number_of_shifts_in);
                                }

                                if(C2Cs[_i_C2C].number_of_shifts_out != NULL){

                                    free(C2Cs[_i_C2C].number_of_shifts_out);
                                }

                                if(C2Cs[_i_C2C].shifts != NULL){

                                    loop_layers{

                                        if(C2Cs[_i_C2C].shifts[i_layer] != NULL){

                                            free(C2Cs[_i_C2C].shifts[i_layer]);
                                        }
                                    }

                                    free(C2Cs[_i_C2C].shifts);
                                }

                                if(C2Cs[_i_C2C].shifts_in != NULL){

                                    loop_layers{

                                        if(C2Cs[_i_C2C].shifts_in[i_layer] != NULL){

                                            free(C2Cs[_i_C2C].shifts_in[i_layer]);
                                        }
                                    }

                                    free(C2Cs[_i_C2C].shifts_in);
                                }

                                if(C2Cs[_i_C2C].shifts_out != NULL){

                                    loop_layers{

                                        if(C2Cs[_i_C2C].shifts_out[i_layer] != NULL){

                                            free(C2Cs[_i_C2C].shifts_out[i_layer]);
                                        }
                                    }

                                    free(C2Cs[_i_C2C].shifts_out);
                                }

                                free_i_2d(C2Cs[_i_C2C].island_offsets);

                                free_i_2d(C2Cs[_i_C2C].island_offsets_in);

                                if(myid == 0){

                                    if(C2Cs[_i_C2C].number_of_shifts_to_node_zero != NULL){

                                        loop_layers{

                                            if(C2Cs[_i_C2C].number_of_shifts_to_node_zero[i_layer] != NULL){

                                                free(C2Cs[_i_C2C].number_of_shifts_to_node_zero[i_layer]);
                                            }
                                        }

                                        free(C2Cs[_i_C2C].number_of_shifts_to_node_zero);
                                    }

                                    if(C2Cs[_i_C2C].number_of_shifts_from_node_zero != NULL){

                                        loop_layers{

                                            if(C2Cs[_i_C2C].number_of_shifts_from_node_zero[i_layer] != NULL){

                                                free(C2Cs[_i_C2C].number_of_shifts_from_node_zero[i_layer]);
                                            }
                                        }

                                        free(C2Cs[_i_C2C].number_of_shifts_from_node_zero);
                                    }

                                    if(C2Cs[_i_C2C].in2out != NULL){

                                        loop_layers{

                                            if(C2Cs[_i_C2C].in2out[i_layer] != NULL){

                                                for(i_node = 0; i_node < (node_last + 1); i_node++){

                                                    if(C2Cs[_i_C2C].in2out[i_layer][i_node] != NULL){

                                                        free(C2Cs[_i_C2C].in2out[i_layer][i_node]);
                                                    }
                                                }

                                                free(C2Cs[_i_C2C].in2out[i_layer]);
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

                free(C2Cs_MPI.shifts_in[i_shift].data);
            }

            free(C2Cs_MPI.shifts_in);
        }

        if(Balance != NULL){

            loop_phases{

                if(Balance[i_phase] != NULL){

                    loop_data{

                        free_r_2d(Balance[i_phase][i_data].node2node_flux);

                        free_r_2d(Balance[i_phase][i_data].node2node_data_flux);
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
