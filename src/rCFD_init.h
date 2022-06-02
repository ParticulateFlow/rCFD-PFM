#ifndef RCFD_INIT
#define RCFD_INIT

#include "rCFD_types.h"
#include "rCFD_globals.h"
#include "rCFD_macros.h"
#include "rCFD_parallel.h"
#include "rCFD_user.h"
#include "rCFD_layer.h"

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
        rCFD_default_Solver_Dict();

        rCFD_user_set_Solver_Dict();

        rCFD_default_Solver();
    }

    /* D2. File_Dict */
    {
        rCFD_default_File_Dict();

        rCFD_user_set_File_Dict();
    }

    /* D3. Phase_Dict */
    {
#if RP_NODE
        Phase_Dict = (Phase_Dict_type*)malloc(Solver_Dict.number_of_phases * sizeof(Phase_Dict_type));

        rCFD_default_Phase_Dict();

        rCFD_user_set_Phase_Dict();
#endif
    }

    /* D4. Tracer_Dict */
    {
#if RP_NODE
        Tracer_Dict.random_walk = (short*)malloc(Solver_Dict.number_of_phases * sizeof(short));

        rCFD_default_Tracer_Dict();

        rCFD_user_set_Tracer_Dict();
#endif
    }

    /* D5. Norm_Dict */
    {
#if RP_NODE
        rCFD_default_Norm_Dict();

        rCFD_user_set_Norm_Dict();
#endif
    }

    /* D6. Rec_Dict */
    {
        rCFD_default_Rec_Dict();

        rCFD_user_set_Rec_Dict();
    }

    /* D7. Data_Dict */
    {
#if RP_NODE
        int i_phase;

        Data_Dict = (Data_Dict_type**)malloc(Solver_Dict.number_of_phases * sizeof(Data_Dict_type*));

        loop_phases{

            Data_Dict[i_phase] = (Data_Dict_type*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(Data_Dict_type));
        }

        rCFD_default_Data_Dict();

        rCFD_user_set_Data_Dict();
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

        rCFD_default_Balance_Dict();

        rCFD_user_set_Balance_Dict();
#endif
    }

    /* D9. Topo_Dict */
    {
        rCFD_default_Topo_Dict();

        rCFD_user_set_Topo_Dict();

#if RP_NODE

        Topo_Dict.Cell_Dict = (Cell_Dict_type*)malloc(Solver_Dict.number_of_layers * sizeof(Cell_Dict_type));

        rCFD_default_Cell_Dict_L0();

        int i_layer;
        
        loop_layers{
            
            rCFD_user_set_Cell_Dict(i_layer);
        }
        
        Topo_Dict.Face_Dict = (Face_Dict_type*)malloc(Solver_Dict.number_of_layers * sizeof(Face_Dict_type));

        rCFD_default_Face_Dict_L0();
        
        loop_layers{
            
            rCFD_user_set_Face_Dict(i_layer);
        }       
#endif
    }

    /* G1. Topo */
    {
        rCFD_default_Topo();

#if RP_NODE
        int i_layer, i_face;

        Topo.Cell = (Cell_type*)malloc(Solver_Dict.number_of_layers * sizeof(Cell_type));

        Topo.Face = (Face_type*)malloc(Solver_Dict.number_of_layers * sizeof(Face_type));

        i_layer = 0;

        /* T.1. allocate cells for L0 */
        {
            rCFD_allocate_layer(i_layer);
        }

        /* T.2. set cell default values for L0 */
        {
            rCFD_default_Cell_L0();
        }

        /* T.3. allocate faces for L0 */
        {
            _F.c0 = (int*)malloc(_Face_Dict.number_of_faces * sizeof(int));

            _F.c1 = (int*)malloc(_Face_Dict.number_of_faces * sizeof(int));

            _F.area = (double**)malloc(_Face_Dict.number_of_faces * sizeof(double*));

            loop_faces{

                _F.area[i_face] = (double*)malloc( 3 * sizeof(double));
            }
        }

        /* T.4. set face default values for L0 */
        {
            rCFD_default_Face_L0();
        }

        /* T.5. set pointers of upper grid layers to NULL */
        if(Solver_Dict.number_of_layers > 1){

            loop_layers_but_L0{

                _C.x                    = NULL;
                _C.volume               = NULL;

                _C.average_velocity     = NULL;
                _C.crossing_time        = NULL;

                _C.hit_by_other_cell    = NULL;
                _C.island_id            = NULL;

                _C.weight_after_shift   = NULL;
                _C.weight_after_swap    = NULL;

                _C.vof                  = NULL;

                _C.data                 = NULL;
                _C.data_shift           = NULL;
                _C.data_swap            = NULL;

                _C.drift_exchange       = NULL;

                _C.user                 = NULL;
                _C.rec_user             = NULL;

                _C.marked               = NULL;
                _C.parent_cell          = NULL;
                _C.number_of_children   = NULL;
                _C.child_index             = NULL;

                _F.c0                   = NULL;
                _F.c1                   = NULL;

                _F.area                 = NULL;
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

        rCFD_default_Tracer();

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

        rCFD_default_Norms();

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

        rCFD_default_Rec();
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
            }
        }
#endif
    }

    /* G1 + G6: Init data fields and update balances */
    {
#if RP_NODE
        int i_layer, i_frame, i_phase, i_data, i_cell;

        i_layer = 0;    /* user initialization just for base layer of grid, might be changed to loop_layers in future */

        rCFD_user_init_Data(i_layer);

        i_frame = 0;

        loop_int_cells{

            loop_phases{

                loop_data{

                    switch (Data_Dict[i_phase][i_data].type){

                        case generic_data:

                            Balance[_i_balance].mass_integral += _C.data[_i_data];

                            break;

                        case concentration_data:

                            Balance[_i_balance].mass_integral += _C.data[_i_data] *

                                Phase_Dict[i_phase].density * _C.volume[i_cell] * _C.vof[_i_vof];

                            break;

                        case temperature_data:

                            Balance[_i_balance].mass_integral += _C.data[_i_data] *

                                Phase_Dict[i_phase].heat_capacity * Phase_Dict[i_phase].density * _C.volume[i_cell] * _C.vof[_i_vof];

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

        rCFD_user_post();    /* post-process initialization */
#endif
    }

    /* P: Parallel grid communication */
    {
#if RP_NODE

        init_parallel_grid_L0();
#endif
    }

    /* L: Initialize layers */
    if(Solver_Dict.number_of_layers > 1){
#if RP_NODE

#if 1   /* local vars */

        short   debug_this_code = 0;

        int     i_layer, i_cell, i_face, i_cell_face, i_face2, i_cell_face2, i_child, i_dim, i_tmp, i_tmp2, i_parent_face, i_redist, i_phase, i_frame;

        int     cell_index_in_upper_layer, number_of_cells_in_upper_layer, upper_cell, c0, c1, c1_c0, c1_c1;

        int     min_dist_parent_cell, max_number_of_faces_per_cell, max_number_of_children, min_number_of_children;

        int     number_of_unassigned_parents, number_of_unassigned_children, number_of_negative_child_indices;

        double  dist, min_dist, x_cluster[3];

        int     **tmp_faces_of_cell = NULL;     /* [i_cell (i_layer)][i_cell_face] */

        int     parent_face_count, upper_c0, upper_c1;

        int     *tmp_parent_face_index = NULL;  /* [i_face] */

        int     **tmp_parent_cell_cell_index = NULL;    /* [i_cell (upper_layer)][i_cell_face] */

        int     **tmp_parent_cell_face_index = NULL;

        short   *tmp_cluster_already_visited = NULL;

        enum{
            not_yet_visited,
            center_of_cluster,
            member_of_cluster,
            neighbor_of_cluster
        };
#endif

        i_layer = 0;

        while(upper_layer < Solver_Dict.number_of_layers)
        {

            /* L.1 create temporary faces_of_cell structure (tmp_faces_of_cell) */
            {

                tmp_faces_of_cell = (int**)malloc(_Cell_Dict.number_of_cells * sizeof(int*));

                loop_cells{

                    tmp_faces_of_cell[i_cell] = (int*)malloc(Solver_Dict.max_number_of_faces_per_cell * sizeof(int));

                    for(i_cell_face = 0; i_cell_face < Solver_Dict.max_number_of_faces_per_cell; i_cell_face++){

                        tmp_faces_of_cell[i_cell][i_cell_face] = -1;
                    }
                }

                max_number_of_faces_per_cell = 0;

                loop_faces{

                    c0 = _F.c0[i_face];
                    c1 = _F.c1[i_face];

                    i_cell_face = 0;

                    while(i_cell_face < Solver_Dict.max_number_of_faces_per_cell){

                        if(tmp_faces_of_cell[c0][i_cell_face] < 0){

                            tmp_faces_of_cell[c0][i_cell_face] = i_face;

                            if((i_cell_face + 1) > max_number_of_faces_per_cell){

                                max_number_of_faces_per_cell = (i_cell_face + 1);
                            }

                            i_cell_face = 2 * Solver_Dict.max_number_of_faces_per_cell;
                        }
                        else{

                            i_cell_face++;
                        }
                    }

                    if(i_cell_face == Solver_Dict.max_number_of_faces_per_cell){

                        Message("\nWARNING: init layers, L1: myid %d  i_layer %d c0 %d i_cell_face %d\n", myid, i_layer, c0, i_cell_face);
                    }


                    i_cell_face = 0;

                    while(i_cell_face < Solver_Dict.max_number_of_faces_per_cell){

                        if(tmp_faces_of_cell[c1][i_cell_face] < 0){

                            tmp_faces_of_cell[c1][i_cell_face] = i_face;

                            if((i_cell_face + 1) > max_number_of_faces_per_cell){

                                max_number_of_faces_per_cell = (i_cell_face + 1);
                            }

                            i_cell_face = 2 * Solver_Dict.max_number_of_faces_per_cell;
                        }
                        else{

                            i_cell_face++;
                        }
                    }

                    if(i_cell_face == Solver_Dict.max_number_of_faces_per_cell){

                        Message("\nWARNING: init layers, L1: myid %d  i_layer %d c1 %d i_cell_face %d\n", myid, i_layer, c1, i_cell_face);
                    }

                }

                if(debug_this_code){

                    Message("\nDEBUG L.1 myid %d i_layer %d max_number_of_faces_per_cell %d", myid, i_layer, max_number_of_faces_per_cell);
                }
            }

            /* L.2 create cell clusters, set number of upper layer cells */
            {
                /* initialize _C.marked and _C.parent_cell */
                loop_cells{

                    _C.marked[i_cell] = not_yet_visited;
                    _C.parent_cell[i_cell] = -1;
                }

                cell_index_in_upper_layer = 0;

                /* mark cell clusters, get number of upper layer cells */
                loop_cells{

                    if(_C.marked[i_cell] == not_yet_visited){

                        _C.marked[i_cell] = center_of_cluster;

                        _C.parent_cell[i_cell] = cell_index_in_upper_layer;

                        /* mark neighbors as members_of_cluster */
                        {
                            i_cell_face = 0;

                            while ((tmp_faces_of_cell[i_cell][i_cell_face] > -1) && (i_cell_face < Solver_Dict.max_number_of_faces_per_cell)){

                                i_face = tmp_faces_of_cell[i_cell][i_cell_face];

                                c0 = _F.c0[i_face];
                                c1 = _F.c1[i_face];

                                if(_C.marked[c1] == center_of_cluster){

                                    i_tmp = c1; c1 = c0; c0 = i_tmp;
                                }

                                _C.marked[c1] = member_of_cluster;

                                _C.parent_cell[c1] = cell_index_in_upper_layer;

                                /* mark neighbors of neighbor as neighbors_of_cluster */
                                {
                                    i_cell_face2 = 0;

                                    while ((tmp_faces_of_cell[c1][i_cell_face2] > -1) && (i_cell_face2 < Solver_Dict.max_number_of_faces_per_cell)){

                                        i_face2 = tmp_faces_of_cell[c1][i_cell_face2];

                                        c1_c0 = _F.c0[i_face2];
                                        c1_c1 = _F.c1[i_face2];

                                        if(_C.marked[c1_c0] == not_yet_visited){

                                            _C.marked[c1_c0] = neighbor_of_cluster;
                                        }
                                        else if(_C.marked[c1_c1] == not_yet_visited){

                                            _C.marked[c1_c1] = neighbor_of_cluster;
                                        }

                                        i_cell_face2++;
                                    }
                                }

                                i_cell_face++;
                            }
                        }

                        cell_index_in_upper_layer++;
                    }
                }

                number_of_cells_in_upper_layer = cell_index_in_upper_layer;

                Topo_Dict.Cell_Dict[upper_layer].number_of_cells = number_of_cells_in_upper_layer;

                Topo_Dict.Cell_Dict[upper_layer].number_of_int_cells = number_of_cells_in_upper_layer;
                Topo_Dict.Cell_Dict[upper_layer].number_of_ext_cells = 0;

                if(debug_this_code){

                    i_tmp = 0;

                    loop_cells{

                        if(_C.marked[i_cell] == neighbor_of_cluster){

                            i_tmp++;
                        }
                    }

                    Message("\nDEBUG L.2 myid %d Topo_Dict.Cell_Dict[%d].number_of_cells %d; unused neighbors %d",

                        myid, upper_layer, Topo_Dict.Cell_Dict[upper_layer].number_of_cells, i_tmp);
                }
            }

            /* L.3 allocate and initialize upper layer grid vars */
            {
                rCFD_allocate_layer(upper_layer);

                /* initialize multi-layer vars */
                loop_cells_of_upper_layer{

                    loop_dim{

                        Topo.Cell[upper_layer].x[i_cell][i_dim] = 0.0;
                    }

                    Topo.Cell[upper_layer].volume[i_cell] = 0.0;

                    Topo.Cell[upper_layer].number_of_children[i_cell] = 0;

                    /* TODO account for multiple islands */
                    Topo.Cell[upper_layer].island_id[i_cell] = 0;
                }

                loop_frames{

                    loop_cells_of_upper_layer{

                        loop_phases{

                            Topo.Cell[upper_layer].vof[_i_vof] = 1.0;
                        }
                    }
                }

                if(debug_this_code){

                    Message("\nDEBUG L.3 myid %d i_layer %d", myid, i_layer);
                }
            }

            /* L.4 fill upper layer grid vars (by first set of clusters) */
            {
                max_number_of_children = 0;
                min_number_of_children = 9999;

                loop_cells{

                    if((_C.marked[i_cell] == center_of_cluster) || (_C.marked[i_cell] == member_of_cluster)){

                        upper_cell = _C.parent_cell[i_cell];

                        loop_dim{

                            Topo.Cell[upper_layer].x[upper_cell][i_dim] += _C.x[i_cell][i_dim];
                        }

                        Topo.Cell[upper_layer].volume[upper_cell] += _C.volume[i_cell];

                        Topo.Cell[upper_layer].number_of_children[upper_cell] += 1;

                        if(Topo.Cell[upper_layer].number_of_children[upper_cell] > max_number_of_children){

                            max_number_of_children = Topo.Cell[upper_layer].number_of_children[upper_cell];
                        }

                        if(Topo.Cell[upper_layer].number_of_children[upper_cell] < min_number_of_children){

                            min_number_of_children = Topo.Cell[upper_layer].number_of_children[upper_cell];
                        }
                    }
                }

                /* average coordinates */
                loop_cells_of_upper_layer{

                    if(Topo.Cell[upper_layer].number_of_children[i_cell] > 0){

                        loop_dim{

                            Topo.Cell[upper_layer].x[i_cell][i_dim] /= (double)Topo.Cell[upper_layer].number_of_children[i_cell];
                        }
                    }
                }

                if(debug_this_code){

                    Message("\nDEBUG L.4 myid %d i_layer %d, max/min_number_of_children %d/%d", myid, i_layer, max_number_of_children, min_number_of_children);
                }

            }

            /* L.5 add still unassigned cells */
            {
                i_tmp2 = 0;

                loop_cells{

                    if(_C.marked[i_cell] == neighbor_of_cluster){

                        i_cell_face = 0;

                        min_dist = 1.0e10;

                        min_dist_parent_cell = 0;

                        while ((tmp_faces_of_cell[i_cell][i_cell_face] > -1) && (i_cell_face < Solver_Dict.max_number_of_faces_per_cell)){

                            i_face = tmp_faces_of_cell[i_cell][i_cell_face];

                            c0 = _F.c0[i_face];
                            c1 = _F.c1[i_face];

                            if(c1 == i_cell){

                                i_tmp = c1; c1 = c0; c0 = i_tmp;
                            }

                            if((_C.marked[c1] == center_of_cluster) || (_C.marked[c1] == member_of_cluster)){

                                /* calc min distance to nearest cluster midpoint */

                                upper_c1 = _C.parent_cell[c1];

                                loop_dim{

                                    x_cluster[i_dim] = Topo.Cell[upper_layer].x[upper_c1][i_dim];
                                }

                                dist = (x_cluster[0] - _C.x[c1][0])*(x_cluster[0] - _C.x[c1][0]) +
                                       (x_cluster[1] - _C.x[c1][1])*(x_cluster[1] - _C.x[c1][1]) +
                                       (x_cluster[2] - _C.x[c1][2])*(x_cluster[2] - _C.x[c1][2]);

                                if(dist < min_dist){

                                    min_dist = dist;

                                    min_dist_parent_cell = _C.parent_cell[c1];
                                }
                            }

                            i_cell_face++;
                        }

                        upper_cell = min_dist_parent_cell;

                        _C.marked[i_cell] = member_of_cluster;

                        _C.parent_cell[i_cell] = upper_cell;

                        loop_dim{

                            Topo.Cell[upper_layer].x[upper_cell][i_dim] = (_C.x[i_cell][0] +

                                (double)Topo.Cell[upper_layer].number_of_children[upper_cell] * Topo.Cell[upper_layer].x[upper_cell][i_dim]) /

                                ((double)Topo.Cell[upper_layer].number_of_children[upper_cell] + 1.0);
                        }

                        Topo.Cell[upper_layer].volume[upper_cell] += _C.volume[i_cell];

                        Topo.Cell[upper_layer].number_of_children[upper_cell] += 1;

                        if(Topo.Cell[upper_layer].number_of_children[upper_cell] > max_number_of_children){

                            max_number_of_children = Topo.Cell[upper_layer].number_of_children[upper_cell];
                        }

                        i_tmp2++;
                    }
                }

                /* Warning, if unassigned cells exist */
                {
                    i_tmp = 0;

                    loop_cells{

                        if(_C.marked[i_cell] == neighbor_of_cluster){

                            i_tmp++;
                        }
                    }

                    if(i_tmp > 0){

                        Message("\nWARNING Initialize layers, myid %d: %d unassigned cells exist in layer %d\n", myid, i_tmp, i_layer);
                    }
                }

                if(debug_this_code){

                    Message("\nDEBUG L.5 myid %d i_layer %d, added cells %d unassigned cells %d max_number_of_children %d", myid, i_layer, i_tmp2, i_tmp, max_number_of_children);
                }

            }

            /* L.5.a fill upper layer vof (TODO island id)*/
            {
                /* init vof(upper_layer) */
                loop_frames{

                    loop_cells{

                        if(_C.parent_cell[i_cell] > -1){

                            upper_cell = _C.parent_cell[i_cell];

                            loop_phases{

                                Topo.Cell[upper_layer].vof[i_frame][upper_cell][i_phase] = 0.0;
                            }
                        }
                    }
                }

                /* accu. map vof(upper_layer) */
                loop_frames{

                    loop_cells{

                        if(_C.parent_cell[i_cell] > -1){

                            upper_cell = _C.parent_cell[i_cell];

                            loop_phases{

                                Topo.Cell[upper_layer].vof[i_frame][upper_cell][i_phase] += _C.vof[_i_vof] * _C.volume[i_cell];
                            }
                        }
                    }
                }

                /* average vof(upper_layer) */
                loop_frames{

                    loop_cells_of_upper_layer{

                        if(Topo.Cell[upper_layer].volume[i_cell] > 0){

                            loop_phases{

                                Topo.Cell[upper_layer].vof[i_frame][i_cell][i_phase] /= Topo.Cell[upper_layer].volume[i_cell];
                            }
                        }
                    }
                }

            }

            /* L.6 fill upper layer list of child_index */
            {
                /* allocate upper layer list of child_index */
                loop_cells_of_upper_layer{

                    Topo.Cell[upper_layer].child_index[i_cell] =

                        (int*)malloc(Topo.Cell[upper_layer].number_of_children[i_cell] * sizeof(int));

                    for(i_child = 0; i_child < Topo.Cell[upper_layer].number_of_children[i_cell]; i_child++){

                        Topo.Cell[upper_layer].child_index[i_cell][i_child] = -1;
                    }
                }

                /* fill upper layer list of child_index */

                number_of_unassigned_parents =      0;
                number_of_unassigned_children =     0;
                number_of_negative_child_indices =  0;

                loop_cells{

                    if(_C.parent_cell[i_cell] > -1){

                        upper_cell = _C.parent_cell[i_cell];

                        i_child = 0;

                        while(i_child < Topo.Cell[upper_layer].number_of_children[upper_cell]){

                            if (Topo.Cell[upper_layer].child_index[upper_cell][i_child] < 0){

                                Topo.Cell[upper_layer].child_index[upper_cell][i_child] = i_cell;

                                i_child =  2 * Topo.Cell[upper_layer].number_of_children[upper_cell];
                            }
                            else{

                                i_child++;
                            }
                        }

                        if(i_child == Topo.Cell[upper_layer].number_of_children[upper_cell]){

                            number_of_unassigned_children++;
                        }
                    }
                    else{

                        number_of_unassigned_parents++;
                    }
                }

                loop_cells_of_upper_layer{

                    for(i_child = 0; i_child < Topo.Cell[upper_layer].number_of_children[i_cell]; i_child++){

                        if(Topo.Cell[upper_layer].child_index[i_cell][i_child] == -1){

                            number_of_negative_child_indices++;
                        }
                    }
                }

                if(debug_this_code){

                    Message("\nDEBUG L.6 myid %d i_layer %d number_of_unassigned_parents %d, number_of_unassigned_children %d neg_child_index %d",

                        myid, i_layer, number_of_unassigned_parents, number_of_unassigned_children, number_of_negative_child_indices);

                    if(myid == 0){

                        FILE *f_out = fopen("tmp", "w");

                        loop_cells_of_upper_layer{

                            fprintf(f_out, "%d :", Topo.Cell[upper_layer].number_of_children[i_cell]);

                            for(i_child = 0; i_child < Topo.Cell[upper_layer].number_of_children[i_cell]; i_child++){

                                fprintf(f_out, " %d", Topo.Cell[upper_layer].child_index[i_cell][i_child]);
                            }

                            fprintf(f_out, "\n");
                        }

                        fclose(f_out);
                    }
                }
            }

            /* L.6.a re-distribute upper-layer's children */
            {
                i_redist = 0;

                while(((max_number_of_children - min_number_of_children) > 2) && (i_redist < Solver_Dict.number_redist_loops_per_layer))
                {
                    /* adapt min clusters */
                    loop_faces{

                        c0 = _F.c0[i_face];
                        c1 = _F.c1[i_face];

                        upper_c0 = _C.parent_cell[c0];
                        upper_c1 = _C.parent_cell[c1];

                        if((Topo.Cell[upper_layer].number_of_children[upper_c1] == min_number_of_children) &&
                           (Topo.Cell[upper_layer].number_of_children[upper_c0] > (min_number_of_children + 1))){

                            i_tmp = c0; c0 = c1; c1 = i_tmp; i_tmp = upper_c0; upper_c0 = upper_c1; upper_c1 = i_tmp;
                        }

                        if((Topo.Cell[upper_layer].number_of_children[upper_c0] == min_number_of_children) &&
                           (Topo.Cell[upper_layer].number_of_children[upper_c1] > (min_number_of_children + 1))){

                            _C.parent_cell[c1] = _C.parent_cell[c0];

                            Topo.Cell[upper_layer].number_of_children[upper_c0] += 1;
                            Topo.Cell[upper_layer].number_of_children[upper_c1] -= 1;
                        }
                    }

                    /* adapt max clusters */
                    loop_faces{

                        c0 = _F.c0[i_face];
                        c1 = _F.c1[i_face];

                        upper_c0 = _C.parent_cell[c0];
                        upper_c1 = _C.parent_cell[c1];

                        if((Topo.Cell[upper_layer].number_of_children[upper_c1] == max_number_of_children) &&
                           (Topo.Cell[upper_layer].number_of_children[upper_c0] < (max_number_of_children - 1))){

                            i_tmp = c0; c0 = c1; c1 = i_tmp; i_tmp = upper_c0; upper_c0 = upper_c1; upper_c1 = i_tmp;
                        }

                        if((Topo.Cell[upper_layer].number_of_children[upper_c0] == max_number_of_children) &&
                           (Topo.Cell[upper_layer].number_of_children[upper_c1] < (max_number_of_children - 1))){

                            _C.parent_cell[c0] = _C.parent_cell[c1];

                            Topo.Cell[upper_layer].number_of_children[upper_c0] -= 1;
                            Topo.Cell[upper_layer].number_of_children[upper_c1] += 1;
                        }
                    }

                    /* adapt clusters in between min and max */
                    {
                        tmp_cluster_already_visited = (short*)malloc(Topo_Dict.Cell_Dict[upper_layer].number_of_cells * sizeof(short));

                        loop_cells_of_upper_layer{

                            tmp_cluster_already_visited[i_cell] = 0;
                        }

                        loop_faces{

                            c0 = _F.c0[i_face];
                            c1 = _F.c1[i_face];

                            upper_c0 = _C.parent_cell[c0];
                            upper_c1 = _C.parent_cell[c1];

                            if((tmp_cluster_already_visited[upper_c0] + tmp_cluster_already_visited[upper_c1]) == 0){

                                if(Topo.Cell[upper_layer].number_of_children[upper_c0] > (Topo.Cell[upper_layer].number_of_children[upper_c1] + 1)){

                                    i_tmp = c0; c0 = c1; c1 = i_tmp; i_tmp = upper_c0; upper_c0 = upper_c1; upper_c1 = i_tmp;
                                }

                                if(Topo.Cell[upper_layer].number_of_children[upper_c0] < (Topo.Cell[upper_layer].number_of_children[upper_c1] - 1)){

                                    _C.parent_cell[c1] = _C.parent_cell[c0];

                                    Topo.Cell[upper_layer].number_of_children[upper_c0] += 1;
                                    Topo.Cell[upper_layer].number_of_children[upper_c1] -= 1;

                                    tmp_cluster_already_visited[upper_c0] = 1;
                                    tmp_cluster_already_visited[upper_c1] = 1;
                                }
                            }
                        }

                        free(tmp_cluster_already_visited);
                    }

                    /* re-allocate upper_layer_children, re-initialize upper layer vars */
                    loop_cells_of_upper_layer{

                        free(Topo.Cell[upper_layer].child_index[i_cell]);

                        Topo.Cell[upper_layer].child_index[i_cell] = (int*)malloc(Topo.Cell[upper_layer].number_of_children[i_cell] * sizeof(int));

                        for(i_child = 0; i_child < Topo.Cell[upper_layer].number_of_children[i_cell]; i_child++){

                            Topo.Cell[upper_layer].child_index[i_cell][i_child] = -1;
                        }

                        loop_dim{

                            Topo.Cell[upper_layer].x[i_cell][i_dim] = 0.0;
                        }

                        Topo.Cell[upper_layer].volume[i_cell] = 0.0;
                    }

                    /* set upper_layer grid vars */
                    loop_cells{

                        upper_cell = _C.parent_cell[i_cell];

                        loop_dim{

                            Topo.Cell[upper_layer].x[upper_cell][i_dim] += _C.x[i_cell][i_dim];
                        }

                        Topo.Cell[upper_layer].volume[upper_cell] += _C.volume[i_cell];

                        /* add child_index */
                        if(upper_cell > -1){

                            i_child = 0;

                            while(i_child < Topo.Cell[upper_layer].number_of_children[upper_cell]){

                                if (Topo.Cell[upper_layer].child_index[upper_cell][i_child] < 0){

                                    Topo.Cell[upper_layer].child_index[upper_cell][i_child] = i_cell;

                                    i_child =  2 * Topo.Cell[upper_layer].number_of_children[upper_cell];
                                }
                                else{

                                    i_child++;
                                }
                            }
                        }
                    }

                    /* average upper-layer coords, get new min/max number_of_children */
                    {
                        max_number_of_children = 0;
                        min_number_of_children = 9999;

                        loop_cells_of_upper_layer{

                            loop_dim{

                                Topo.Cell[upper_layer].x[i_cell][i_dim] /= (double)Topo.Cell[upper_layer].number_of_children[i_cell];
                            }

                            if(Topo.Cell[upper_layer].number_of_children[i_cell] > max_number_of_children){

                                max_number_of_children = Topo.Cell[upper_layer].number_of_children[i_cell];
                            }

                            if(Topo.Cell[upper_layer].number_of_children[i_cell] < min_number_of_children){

                                min_number_of_children = Topo.Cell[upper_layer].number_of_children[i_cell];
                            }
                        }
                    }

                    if(debug_this_code){

                        Message("\nDEBUG myid %d: L.6.a: i_redist %d, min/max number_of_children %d/%d", myid, i_redist, min_number_of_children, max_number_of_children);
                    }

                    i_redist++;
                }
            }

            /* L.7 create tmp_parent_xy vars for upper_layer faces */
            {
                tmp_parent_face_index = (int*)malloc(_Face_Dict.number_of_faces * sizeof(int));

                loop_faces{

                    tmp_parent_face_index[i_face] = -1;
                }

                tmp_parent_cell_cell_index = (int**)malloc(Topo_Dict.Cell_Dict[upper_layer].number_of_cells * sizeof(int*));

                loop_cells_of_upper_layer{

                    tmp_parent_cell_cell_index[i_cell] = (int*)malloc(Solver_Dict.max_number_of_faces_per_cell * sizeof(int));

                    for(i_cell_face = 0; i_cell_face < Solver_Dict.max_number_of_faces_per_cell; i_cell_face++){

                        tmp_parent_cell_cell_index[i_cell][i_cell_face] = -1;
                    }
                }

                tmp_parent_cell_face_index = (int**)malloc(Topo_Dict.Cell_Dict[upper_layer].number_of_cells * sizeof(int*));

                loop_cells_of_upper_layer{

                    tmp_parent_cell_face_index[i_cell] = (int*)malloc(Solver_Dict.max_number_of_faces_per_cell * sizeof(int));

                    for(i_cell_face = 0; i_cell_face < Solver_Dict.max_number_of_faces_per_cell; i_cell_face++){

                        tmp_parent_cell_face_index[i_cell][i_cell_face] = -1;
                    }
                }

                if(debug_this_code){

                    Message("\nDEBUG L.7 myid %d i_layer %d", myid, i_layer);
                }

            }

            /* L.8 fill tmp_parent_xy vars */
            {
                parent_face_count = 0;

                loop_faces{

                    c0 = _F.c0[i_face];
                    c1 = _F.c1[i_face];

                    upper_c0 = _C.parent_cell[c0];
                    upper_c1 = _C.parent_cell[c1];

                    if(upper_c0 == upper_c1){

                        tmp_parent_face_index[i_face] = -1;
                    }
                    else{

                        i_cell_face = 0;

                        while(i_cell_face < Solver_Dict.max_number_of_faces_per_cell){

                            /* parent face already exists */
                            if(tmp_parent_cell_cell_index[upper_c0][i_cell_face] == upper_c1){

                                tmp_parent_face_index[i_face] = tmp_parent_cell_face_index[upper_c0][i_cell_face];

                                i_cell_face = Solver_Dict.max_number_of_faces_per_cell;
                            }

                            /* add new parent face to upper_c0, upper_c1 */
                            if(tmp_parent_cell_cell_index[upper_c0][i_cell_face] < 0){

                                tmp_parent_cell_cell_index[upper_c0][i_cell_face] = upper_c1;

                                tmp_parent_cell_face_index[upper_c0][i_cell_face] = parent_face_count;

                                /* also add this parent face to the list of upper_c1 */
                                {
                                    i_cell_face2 = 0;

                                    while((tmp_parent_cell_cell_index[upper_c1][i_cell_face2] >= 0) &&

                                        (i_cell_face < Solver_Dict.max_number_of_faces_per_cell)){

                                        i_cell_face2++;
                                    }

                                    if(i_cell_face2 < Solver_Dict.max_number_of_faces_per_cell){

                                        tmp_parent_cell_cell_index[upper_c1][i_cell_face2] = upper_c0;

                                        tmp_parent_cell_face_index[upper_c1][i_cell_face2] = parent_face_count;
                                    }
                                }

                                tmp_parent_face_index[i_face] = parent_face_count;

                                parent_face_count++;

                                i_cell_face = Solver_Dict.max_number_of_faces_per_cell;
                            }

                            i_cell_face++;
                        }
                    }
                }

                Topo_Dict.Face_Dict[upper_layer].number_of_faces = parent_face_count;

                /* TODO account for upper-layer parallel grid communication */
                Topo_Dict.Face_Dict[upper_layer].number_of_int_faces = parent_face_count;
                Topo_Dict.Face_Dict[upper_layer].number_of_ext_faces = 0;

                if(debug_this_code){

                    Message("\nDEBUG L.8 myid %d, upper_layer %d, number of faces in upper layer %d", myid, upper_layer, parent_face_count);
                }
            }

            /* L.9 create and fill upper faces */
            {
                Topo.Face[upper_layer].c0 = (int*)malloc(Topo_Dict.Face_Dict[upper_layer].number_of_faces * sizeof(int));
                Topo.Face[upper_layer].c1 = (int*)malloc(Topo_Dict.Face_Dict[upper_layer].number_of_faces * sizeof(int));

                Topo.Face[upper_layer].area = (double**)malloc(Topo_Dict.Face_Dict[upper_layer].number_of_faces * sizeof(double*));

                loop_faces_of_upper_layer{

                    Topo.Face[upper_layer].area[i_face] = (double*)malloc (3 * sizeof(double));
                }

                loop_faces{

                    if(tmp_parent_face_index[i_face] >= 0){

                        c0 = _F.c0[i_face];
                        c1 = _F.c1[i_face];

                        upper_c0 = _C.parent_cell[c0];
                        upper_c1 = _C.parent_cell[c1];

                        i_parent_face = tmp_parent_face_index[i_face];

                        Topo.Face[upper_layer].c0[i_parent_face] = upper_c0;
                        Topo.Face[upper_layer].c1[i_parent_face] = upper_c1;

                        Topo.Face[upper_layer].area[i_parent_face][0] += _F.area[i_face][0];
                        Topo.Face[upper_layer].area[i_parent_face][1] += _F.area[i_face][1];
                        Topo.Face[upper_layer].area[i_parent_face][2] += _F.area[i_face][2];
                    }
                }

                if(debug_this_code){

                    Message("\nDEBUG L.9 myid %d i_layer %d", myid, i_layer);
                }

            }

            /* L.10 consistency checks */
            {

                if(Topo_Dict.Cell_Dict[upper_layer].number_of_cells < Solver_Dict.min_number_of_cells_per_layer){

                    Message("\nWARNING - myid %d: layer-%d number_of_cells (%d) is smaller than Solver_Dict.min_number_of_cells_per_layer (%d)\n",

                        myid, upper_layer, Topo_Dict.Cell_Dict[upper_layer].number_of_cells, Solver_Dict.min_number_of_cells_per_layer);
                }

                if(max_number_of_children > Solver_Dict.max_number_of_children){

                    Message0("\nWARNING - myid %d: layer-%d number_of_children (%d) is larger than Solver_Dict.max_number_of_children (%d)\n",

                        myid, upper_layer, max_number_of_children, Solver_Dict.max_number_of_children);
                }

                if((number_of_unassigned_children + number_of_unassigned_parents + number_of_negative_child_indices) > 0){

                    Message("\nWARNING - myid %d: Inconsistency between layer %d and %d:  unassigned_children (%d) unassigned_parents (%d) negative_child_indices (%d)\n",

                        myid, i_layer, upper_layer, number_of_unassigned_children, number_of_unassigned_parents, number_of_negative_child_indices);
                }

                if(debug_this_code){

                    Message("\nDEBUG L.10 myid %d i_layer %d", myid, i_layer);
                }
            }

            /* L.11 free local vars */
            {
                if(tmp_faces_of_cell != NULL){

                    loop_cells{

                        if(tmp_faces_of_cell[i_cell] != NULL){

                            free(tmp_faces_of_cell[i_cell]);
                        }
                    }

                    free(tmp_faces_of_cell);
                }

                if(tmp_parent_face_index != NULL){

                    free(tmp_parent_face_index);
                }

                if(tmp_parent_cell_cell_index != NULL){

                    loop_cells_of_upper_layer{

                        if(tmp_parent_cell_cell_index[i_cell] != NULL){

                            free(tmp_parent_cell_cell_index[i_cell]);
                        }
                    }

                    free(tmp_parent_cell_cell_index);
                }

                if(tmp_parent_cell_face_index != NULL){

                    loop_cells_of_upper_layer{

                        if(tmp_parent_cell_face_index[i_cell] != NULL){

                            free(tmp_parent_cell_face_index[i_cell]);
                        }
                    }

                    free(tmp_parent_cell_face_index);
                }

                if(debug_this_code){

                    Message("\nDEBUG L.11 myid %d i_layer %d ", myid, i_layer);
                }

            }

            i_layer++;
        }

        if(Solver_Dict.verbose){

            Message("\nmulti-layer grid: myid %d\n", myid);

            loop_layers{

                Message("   Layer %d: Cells %d Faces %d\n", i_layer, _Cell_Dict.number_of_cells, _Face_Dict.number_of_faces);
            }
        }

        /* Test mapping */
        if(0==1){
            i_layer = 0;

            rCFD_map_parent(i_layer);
            /* layer-1 should have number of layer-0 children in data[0] */

            i_layer = 1;

            rCFD_map_parent(i_layer);
            /* layer-2 should have number of their layer-0 children in data[0] */
#if 0
            i_layer = 2;

            rCFD_map_parent(i_layer);
            /* layer-3 cells should have their layer-2 siblings */

            i_layer = 3;

            rCFD_map_children(i_layer);
            /* layer-1 cells should have their layer-2 siblings */
#endif
            i_layer = 2;

            rCFD_map_children(i_layer);
            /* layer-0 cells should have their layer-2 siblings */

            i_layer = 1;

            rCFD_map_children(i_layer);
            /* layer-0 cells should have their layer-2 siblings */

            i_layer = 0;

            rCFD_user_post();
            /* data[0] should visualize clusters */
        }
#endif
    }

}
#endif
