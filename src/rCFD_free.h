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
        free(Rec.global_frame);

        free_i_4d(Rec.jumps);

        free(Solver.timestep_width);

#if RP_NODE
        int     i_phase, i_frame, i_data, i_shift, i_node, i_cell, i_layer;

        int     number_of_shifts;

        free(Tracer_Dict.random_walk);

        if(Topo.Cell != NULL){

            loop_layers{

                free_r_2d(_C.x);

                free(_C.volume);

                free_r_2d(_C.average_velocity);

                free_r_2d(_C.crossing_time);

                free(_C.hit_by_other_cell);

                free(_C.island_id);

                free(_C.weight_after_shift);

                free(_C.weight_after_swap);

                free_r_3d(_C.vof);

                if(_C.data != NULL){

                    loop_phases{

                        if(_C.data[i_phase] != NULL){

                            loop_cells{

                                free(_C.data[i_phase][i_cell]);
                            }

                            free(_C.data[i_phase]);
                        }
                    }

                    free(_C.data);
                }

                if(_C.data_shift != NULL){

                    loop_phases{

                        if(_C.data_shift[i_phase] != NULL){

                            loop_cells{

                                free(_C.data_shift[i_phase][i_cell]);
                            }

                            free(_C.data_shift[i_phase]);
                        }
                    }

                    free(_C.data_shift);
                }

                if(_C.data_swap != NULL){

                    loop_phases{

                        if(_C.data_swap[i_phase] != NULL){

                            loop_cells{

                                free(_C.data_swap[i_phase][i_cell]);
                            }

                            free(_C.data_swap[i_phase]);
                        }
                    }

                    free(_C.data_swap);
                }

                free(_C.drift_exchange);

                free_r_2d(_C.user);

                free(_C.marked);

                free(_C.parent_cell);

                free(_C.number_of_children);

                if(_C.child_index != NULL){

                    loop_cells{

                        free(_C.child_index[i_cell]);
                    }

                    free(_C.child_index);
                }
            }

            free(Topo.Cell);
        }

        if(Topo.Face != NULL){

            loop_layers{

                free(_F.c0);

                free(_F.c1);

                free_r_2d(_F.area);
            }

            free(Topo.Face);
        }

        free(Topo_Dict.Cell_Dict);

        free(Topo_Dict.Face_Dict);


        free(Tracer.monitor_counter);

        free(Tracer.number_of_shifts);

        if(Tracer.shifts != NULL){

            loop_phases{

                free(Tracer.shifts[i_phase]);
            }

            free(Tracer.shifts);
        }

        free(Norms.norm);

        if(C2Cs != NULL){

            int i_state;

            loop_states{

                if(C2Cs[i_state] != NULL){

                    loop_phases{

                        if(C2Cs[i_state][i_phase] != NULL){

                            loop_frames{

                                free(C2Cs[_i_C2C].number_of_shifts);

                                free(C2Cs[_i_C2C].number_of_shifts_in);

                                free(C2Cs[_i_C2C].number_of_shifts_out);

                                if(C2Cs[_i_C2C].shifts != NULL){

                                    loop_layers{

                                        free(C2Cs[_i_C2C].shifts[i_layer]);
                                    }

                                    free(C2Cs[_i_C2C].shifts);
                                }

                                if(C2Cs[_i_C2C].shifts_in != NULL){

                                    loop_layers{

                                        free(C2Cs[_i_C2C].shifts_in[i_layer]);
                                    }

                                    free(C2Cs[_i_C2C].shifts_in);
                                }

                                if(C2Cs[_i_C2C].shifts_out != NULL){

                                    loop_layers{

                                        free(C2Cs[_i_C2C].shifts_out[i_layer]);
                                    }

                                    free(C2Cs[_i_C2C].shifts_out);
                                }

                                free_i_2d(C2Cs[_i_C2C].island_offsets);

                                free_i_2d(C2Cs[_i_C2C].island_offsets_in);

                                if(myid == 0){

                                    if(C2Cs[_i_C2C].number_of_shifts_to_node_zero != NULL){

                                        loop_layers{

                                            free(C2Cs[_i_C2C].number_of_shifts_to_node_zero[i_layer]);
                                        }

                                        free(C2Cs[_i_C2C].number_of_shifts_to_node_zero);
                                    }

                                    if(C2Cs[_i_C2C].number_of_shifts_from_node_zero != NULL){

                                        loop_layers{

                                            free(C2Cs[_i_C2C].number_of_shifts_from_node_zero[i_layer]);
                                        }

                                        free(C2Cs[_i_C2C].number_of_shifts_from_node_zero);
                                    }

                                    if(C2Cs[_i_C2C].in2out != NULL){

                                        loop_layers{

                                            if(C2Cs[_i_C2C].in2out[i_layer] != NULL){

                                                for(i_node = 0; i_node < (node_last + 1); i_node++){

                                                    free(C2Cs[_i_C2C].in2out[i_layer][i_node]);
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

                free(Balance_Dict[i_phase]);
            }

            free(Balance_Dict);
        }

        if(Data_Dict != NULL){

            loop_phases{

                if(Data_Dict[i_phase] != NULL){

                    loop_data{

                        free(Data_Dict[i_phase][i_data].drift_velocity);

                        free(Data_Dict[i_phase][i_data].user);
                    }

                    free(Data_Dict[i_phase]);
                }
            }

            free(Data_Dict);
        }

        /* free last, because previously loop_data was needed */
        if(Phase_Dict != NULL){

            loop_phases{

                free(Phase_Dict[i_phase].user);
            }

            free(Phase_Dict);
        }
#endif
    }
#endif
