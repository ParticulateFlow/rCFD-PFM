#ifndef RCFD_LAYERS
#define RCFD_LAYERS

#include "rCFD_macros.h"
#include "rCFD_globals.h"


/* (C)  2022
    Stefan Pirker
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/


    void rCFD_allocate_layer(const short i_layer)
    {
    #if RP_NODE

        int i_cell, i_phase;

        _C.x = malloc_r_2d(_Cell_Dict.number_of_cells, 3);

        _C.volume = (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));
        _C.grid_spacing = (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));

        _C.average_velocity = malloc_r_2d(Solver_Dict.number_of_phases, _Cell_Dict.number_of_cells);
        _C.crossing_time    = malloc_r_2d(Solver_Dict.number_of_phases, _Cell_Dict.number_of_cells);

        _C.hit_by_other_cell =   (short*)malloc(_Cell_Dict.number_of_cells * sizeof(short));

        _C.island_id =           (short*)malloc(_Cell_Dict.number_of_cells * sizeof(short));

        _C.weight_after_shift =  (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));
        _C.weight_after_swap =   (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));

        _C.vof = malloc_r_3d(Solver_Dict.number_of_frames, _Cell_Dict.number_of_cells, Solver_Dict.number_of_phases);

        _C.data =       (double***)malloc(Solver_Dict.number_of_phases * sizeof(double**));
        _C.data_shift = (double***)malloc(Solver_Dict.number_of_phases * sizeof(double**));
        _C.data_swap =  (double***)malloc(Solver_Dict.number_of_phases * sizeof(double**));

        loop_phases{

            _C.data[i_phase] =       (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
            _C.data_shift[i_phase] = (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
            _C.data_swap[i_phase] =  (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));

            loop_cells{

                _C.data[i_phase][i_cell] =       (double*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(double));
                _C.data_shift[i_phase][i_cell] = (double*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(double));
                _C.data_swap[i_phase][i_cell] =  (double*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(double));
            }
        }

        if(Solver_Dict.data_drifting_on){

            _C.drift_exchange = (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));
        }
        else{

            _C.drift_exchange = NULL;
        }

        if(_Cell_Dict.number_of_user_vars > 0){

            _C.user = malloc_r_2d(_Cell_Dict.number_of_cells, _Cell_Dict.number_of_user_vars);
        }
        else{

            _C.user = NULL;
        }

        /* multi-layer vars */

        _C.marked =                 (short*)malloc(_Cell_Dict.number_of_cells * sizeof(short));

        if((i_layer + 1) < Solver_Dict.number_of_layers){

            _C.parent_cell =        (int*)malloc(_Cell_Dict.number_of_cells * sizeof(int));
        }
        else{

            _C.parent_cell =        NULL;
        }

        if(i_layer > 0){

            _C.number_of_children = (short*)malloc(_Cell_Dict.number_of_cells * sizeof(short));
            _C.child_index =           (int**) malloc(_Cell_Dict.number_of_cells * sizeof(int*));      /* child_index[i_cell] will be allocated later */
        }
        else{

            _C.number_of_children = NULL;
            _C.child_index =           NULL;
        }
    #endif
    }

    void rCFD_map_parent(const short i_layer)
    {
#if RP_NODE
        short debug_this_code = 0;

        if(upper_layer >= Solver_Dict.number_of_layers){

            Message("\nWARNING myid %d: Tried to map to non-existing parent of layer %d", myid, i_layer);

            return;
        }

        int i_cell, i_data, i_phase;

        int upper_cell;

        /* initialize upper-layer data */
        loop_phases{

            loop_cells_of_upper_layer{

                loop_data{

                    Topo.Cell[upper_layer].data[_i_data] = 0.0;
                }
            }
        }

        /* accumulate upper-layer data */
        loop_phases{

            loop_cells{

                upper_cell = _C.parent_cell[i_cell];

                loop_data{

                    Topo.Cell[upper_layer].data[i_phase][upper_cell][i_data] += _C.data[_i_data];
                }
            }
        }

        /* average upper-layer data */
        loop_phases{

            loop_cells_of_upper_layer{

                loop_data{

                    Topo.Cell[upper_layer].data[_i_data] /= (double) Topo.Cell[upper_layer].number_of_children[i_cell];
                }
            }
        }

        if(debug_this_code){

            i_data = 0;

            i_phase = 0;

            loop_cells_of_upper_layer{

                Topo.Cell[upper_layer].data[_i_data] = 0.0;
            }

            loop_cells{

                upper_cell = _C.parent_cell[i_cell];

                Topo.Cell[upper_layer].data[i_phase][upper_cell][i_data] += 1.0;
            }
        }
#endif
    }

    void rCFD_map_children(const short i_layer)
    {
#if RP_NODE
        short debug_this_code = 0;

        if(i_layer < 1){

            Message("\nWARNING myid %d: Tried to map to non-existing children of layer %d", myid, i_layer);

            return;
        }

        int i_cell, i_data, i_phase, i_child;

        int lower_cell;

        loop_phases{

            loop_cells{

                loop_data{

                    loop_children{

                        lower_cell = _C.child_index[i_cell][i_child];

                        Topo.Cell[lower_layer].data[i_phase][lower_cell][i_data] = _C.data[_i_data];
                    }
                }
            }
        }

        if(debug_this_code){

            i_data = 0;

            i_phase = 0;

            loop_cells_of_lower_layer{

                Topo.Cell[lower_layer].data[_i_data] = 0.0;
            }

            loop_cells{

                loop_children{

                    lower_cell = _C.child_index[i_cell][i_child];

                    Topo.Cell[lower_layer].data[i_phase][lower_cell][i_data] = _C.data[_i_data];
                }
            }
        }
#endif
    }

    void rCFD_map_from_to_layer(const short from_layer, const short to_layer)
    {
#if RP_NODE

        if((to_layer >= Solver_Dict.number_of_layers) || (to_layer < 0)){

            Message("\nWARNING myid %d: Tried to map to non-existing layer %d", myid, to_layer);

            return;
        }

        if((from_layer >= Solver_Dict.number_of_layers) || (from_layer < 0)){

            Message("\nWARNING myid %d: Unreasonable mapping layer %d", myid, from_layer);

            return;
        }

        int i_layer;

        if(from_layer > to_layer){

            i_layer = from_layer;

            while(i_layer > to_layer){

                rCFD_map_children(i_layer);

                i_layer--;
            }
        }
        else if(from_layer < to_layer){

            i_layer = from_layer;

            while(i_layer < to_layer){

                rCFD_map_parent(i_layer);

                i_layer++;
            }
        }
#endif
    }

#endif
