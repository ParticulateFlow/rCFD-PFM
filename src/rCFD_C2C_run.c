#include <stdio.h>
#include <string.h> //memset

#include "rCFD_types.h"
#include "rCFD_globals.h"
#include "rCFD_memory.h"
#include "rCFD_parallel.h"
#include "rCFD_defaults.h"
#include "rCFD_macros.h"
#include "rCFD_init.h"
#include "rCFD_layer.h"
#include "rCFD_free.h"

#include "rCFD_user.h"

/* (C)  2021-22
    Stefan Pirker
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/

/* Contents:

    rCFD_init_all

    rCFD_read_C2Cs

    rCFD_run

    rCFD_free_all
*/

/*************************************************************************************/
DEFINE_ON_DEMAND(rCFD_init_all)
/*************************************************************************************/
{
    init_all();

    Message0("\n\n...rCFD_init_all\n");
}

/*************************************************************************************/
DEFINE_ON_DEMAND(rCFD_read_C2Cs)
/*************************************************************************************/
{
#if RP_NODE
    int     i_state, i_phase, i_frame, i_shift, i_read, i_read_int, i_read_incoming, i_island, i_layer;

    int     C2C_index;

    int     *i_island_int, *i_island_incoming;

    int     total_number_of_C2Cs_read = 0;

    tmp_C2C_type    tmp_C2Cs;
#endif

    Solver.clock = clock();

    /* A. Allocate C2C's and local indices */
    {
#if RP_NODE

        /* global C2Cs [state, phase, frame] */
        {
            C2Cs = (C2C_type***)malloc(Solver_Dict.number_of_states * sizeof(C2C_type**));

            loop_states{

                C2Cs[i_state] = (C2C_type**)malloc(Solver_Dict.number_of_phases * sizeof(C2C_type*));

                loop_phases{

                    C2Cs[i_state][i_phase] = (C2C_type*)malloc(Solver_Dict.number_of_frames * sizeof(C2C_type));

                    loop_frames{

                        /* init C2Cs structs */
                        C2Cs[_i_C2C].format = c0_n0_c1_n1_w0;

                        C2Cs[_i_C2C].number_of_shifts =        (int*)malloc(Solver_Dict.number_of_layers * sizeof(int));
                        C2Cs[_i_C2C].number_of_shifts_in =     (int*)malloc(Solver_Dict.number_of_layers * sizeof(int));
                        C2Cs[_i_C2C].number_of_shifts_out =    (int*)malloc(Solver_Dict.number_of_layers * sizeof(int));

                        loop_layers{

                            C2Cs[_i_C2C].number_of_shifts[i_layer] =        0;
                            C2Cs[_i_C2C].number_of_shifts_in[i_layer] =     0;
                            C2Cs[_i_C2C].number_of_shifts_out[i_layer] =    0;
                        }

                        C2Cs[_i_C2C].shifts =       (C2C_shift_type**)malloc(Solver_Dict.number_of_layers * sizeof(C2C_shift_type*));
                        C2Cs[_i_C2C].shifts_in =    (C2C_shift_type**)malloc(Solver_Dict.number_of_layers * sizeof(C2C_shift_type*));
                        C2Cs[_i_C2C].shifts_out =   (C2C_shift_type**)malloc(Solver_Dict.number_of_layers * sizeof(C2C_shift_type*));

                        loop_layers{

                            C2Cs[_i_C2C].shifts[i_layer] =      NULL;
                            C2Cs[_i_C2C].shifts_in[i_layer] =   NULL;
                            C2Cs[_i_C2C].shifts_out[i_layer] =  NULL;
                        }

                        C2Cs[_i_C2C].island_offsets    = malloc_i_2d(Solver_Dict.number_of_layers, Solver_Dict.number_of_islands+1);
                        C2Cs[_i_C2C].island_offsets_in = malloc_i_2d(Solver_Dict.number_of_layers, Solver_Dict.number_of_islands+1);

                        memset(C2Cs[_i_C2C].island_offsets[0],    0, Solver_Dict.number_of_layers * (Solver_Dict.number_of_islands+1) * sizeof(int));
                        memset(C2Cs[_i_C2C].island_offsets_in[0], 0, Solver_Dict.number_of_layers * (Solver_Dict.number_of_islands+1) * sizeof(int));


                        /* following vars will be allocated in rCFD_parallel.h */

                        C2Cs[_i_C2C].number_of_shifts_to_node_zero = NULL;
                        C2Cs[_i_C2C].number_of_shifts_from_node_zero = NULL;
                        C2Cs[_i_C2C].in2out = NULL;
                    }
                }
            }
        }

        /* local island counters */
        {
            i_island_int = (int*)malloc(Solver_Dict.number_of_islands*sizeof(int));
            i_island_incoming = (int*)malloc(Solver_Dict.number_of_islands*sizeof(int));
        }

#endif
    }

    /* B. Read files into tmp_shifts; on the fly analysis of number_of_shifts and island_offsets, write tmp_C2Cs into C2Cs */
    {
#if RP_NODE
        char    filename[40];
        FILE    *fi = NULL;
        int     i_tmp, c0, node0, c1, node1, my_island_id, i_layer;
        double  tmp, w0;

        i_layer = 0;

        loop_states{

            loop_phases{

                /* try to open file[state][phase] */
                {
                    sprintf(filename,"%s_%d_%d_%d", File_Dict.C2C_filename, i_state, i_phase, myid);

                    fi = fopen(filename,"r");
                }

                /* read C2Cs per [state, phase, frame] */
                if(fi){

                    loop_frames{

                        /* B.1. initializations per frame */
                        {
                            i_read = 0; i_read_int = 0; i_read_incoming = 0;

                            loop_islands{

                                i_island_int[i_island] = 0;
                                i_island_incoming[i_island] = 0;
                            }
                        }

                        /* B.2. allocate tmp_C2Cs, read tmp_C2Cs, count indices */
                        {
                            fscanf(fi,"%d \n",&tmp_C2Cs.number_of_shifts);

                            /* in future, we might read this from file */
                            tmp_C2Cs.format = c0_n0_c1_n1_w0;

                            /* allocate tmp_C2Cs.shifts */
                            {
                                tmp_C2Cs.shifts = (C2C_shift_type*)malloc(tmp_C2Cs.number_of_shifts * sizeof(C2C_shift_type));

                                for(i_shift = 0; i_shift < tmp_C2Cs.number_of_shifts; i_shift++){

                                    tmp_C2Cs.shifts[i_shift].c0 =     -1;
                                    tmp_C2Cs.shifts[i_shift].node0 =  -1;
                                    tmp_C2Cs.shifts[i_shift].c1 =     -1;
                                    tmp_C2Cs.shifts[i_shift].node1 =  -1;
                                    tmp_C2Cs.shifts[i_shift].w0 =     0.0;
                                }

                                tmp_C2Cs.island_offsets =       (int*)malloc((Solver_Dict.number_of_islands+1) * sizeof(int));
                                tmp_C2Cs.island_offsets_in =    (int*)malloc((Solver_Dict.number_of_islands+1) * sizeof(int));

                                loop_islands{

                                    tmp_C2Cs.island_offsets[i_island] =     0;
                                    tmp_C2Cs.island_offsets_in[i_island] =  0;
                                }

                            }

                            /* fill tmp_C2Cs and monitor indices */
                            for(i_shift = 0; i_shift < tmp_C2Cs.number_of_shifts; i_shift++){

                                if((i_shift % Solver_Dict.C2C_loading_reduction) == 0){

                                    if(tmp_C2Cs.format == c0_n0_c1_n1_w0){

                                        fscanf(fi,"%d %d %d %d %le\n", &c0, &node0, &c1, &node1, &w0);
                                    }

                                    tmp_C2Cs.shifts[i_read].c0 = c0;
                                    tmp_C2Cs.shifts[i_read].node0 = node0;
                                    tmp_C2Cs.shifts[i_read].c1 = c1;
                                    tmp_C2Cs.shifts[i_read].node1 = node1;
                                    tmp_C2Cs.shifts[i_read].w0 = w0;

                                    i_read++;

                                    my_island_id = _C.island_id[c1];

                                    if(node0 == node1){

                                        i_read_int++;
                                        i_island_int[my_island_id]++;
                                    }
                                    else{

                                        i_read_incoming++;
                                        i_island_incoming[my_island_id]++;
                                    }
                                }
                                else{

                                    if(tmp_C2Cs.format == c0_n0_c1_n1_w0){

                                        fscanf(fi,"%d %d %d %d %le\n", &i_tmp, &i_tmp, &i_tmp, &i_tmp, &tmp);
                                    }
                                }
                            }
                        }

                        /* B.3 set i_read_outgoing, set island indices */
                        {
                            tmp_C2Cs.island_offsets[0] =        0;
                            tmp_C2Cs.island_offsets_in[0] =     0;

                            loop_islands{

                                tmp_C2Cs.island_offsets[(i_island + 1)] =
                                    tmp_C2Cs.island_offsets[i_island] + i_island_int[i_island];

                                tmp_C2Cs.island_offsets_in[(i_island + 1)] =
                                    tmp_C2Cs.island_offsets_in[i_island] + i_island_incoming[i_island];
                            }
                        }

                        /* B.4 fill tmp_C2Cs into C2Cs.shifts ... */
                        {
                            /* B.4.1. initializations per frame */
                            {
                                loop_islands{

                                    i_island_int[i_island] = 0;
                                    i_island_incoming[i_island] = 0;
                                }
                            }

                            /* B.4.2. C2C.numbers, allocate and init C2Cs.shifts/shifts_in */
                            {
                                C2Cs[_i_C2C].number_of_shifts[0] = i_read_int;

                                C2Cs[_i_C2C].shifts[0] = (C2C_shift_type*)malloc(i_read_int * sizeof(C2C_shift_type));

                                for(i_shift = 0; i_shift < C2Cs[_i_C2C].number_of_shifts[0]; i_shift++){

                                    C2Cs[_i_C2C].shifts[0][i_shift].c0 =        -1;
                                    C2Cs[_i_C2C].shifts[0][i_shift].node0 =     -1;
                                    C2Cs[_i_C2C].shifts[0][i_shift].c1 =        -1;
                                    C2Cs[_i_C2C].shifts[0][i_shift].node1 =     -1;
                                    C2Cs[_i_C2C].shifts[0][i_shift].w0 =        0.0;
                                }

                                C2Cs[_i_C2C].number_of_shifts_in[0] = i_read_incoming;

                                C2Cs[_i_C2C].shifts_in[0] = (C2C_shift_type*)malloc(i_read_incoming * sizeof(C2C_shift_type));

                                for(i_shift = 0; i_shift < C2Cs[_i_C2C].number_of_shifts_in[0]; i_shift++){

                                    C2Cs[_i_C2C].shifts_in[0][i_shift].c0 =     -1;
                                    C2Cs[_i_C2C].shifts_in[0][i_shift].node0 =  -1;
                                    C2Cs[_i_C2C].shifts_in[0][i_shift].c1 =     -1;
                                    C2Cs[_i_C2C].shifts_in[0][i_shift].node1 =  -1;
                                    C2Cs[_i_C2C].shifts_in[0][i_shift].w0 =     0.0;
                                }

                                total_number_of_C2Cs_read +=

                                    C2Cs[_i_C2C].number_of_shifts[0] + C2Cs[_i_C2C].number_of_shifts_in[0];
                            }

                            /* B.4.3. fill C2Cs.shifts/shifts_in by tmp_C2Cs*/
                            {
                                for(i_shift = 0; i_shift < tmp_C2Cs.number_of_shifts; i_shift++){

                                    c1 = tmp_C2Cs.shifts[i_shift].c1;

                                    my_island_id = _C.island_id[c1];

                                    node0 = tmp_C2Cs.shifts[i_shift].node0;
                                    node1 = tmp_C2Cs.shifts[i_shift].node1;


                                    if(node0 == node1){

                                        C2C_index = tmp_C2Cs.island_offsets[my_island_id] + i_island_int[my_island_id];

                                        C2Cs[_i_C2C].shifts[0][C2C_index].c0 =    tmp_C2Cs.shifts[i_shift].c0;
                                        C2Cs[_i_C2C].shifts[0][C2C_index].node0 = tmp_C2Cs.shifts[i_shift].node0;
                                        C2Cs[_i_C2C].shifts[0][C2C_index].c1 =    tmp_C2Cs.shifts[i_shift].c1;
                                        C2Cs[_i_C2C].shifts[0][C2C_index].node1 = tmp_C2Cs.shifts[i_shift].node1;
                                        C2Cs[_i_C2C].shifts[0][C2C_index].w0 =    tmp_C2Cs.shifts[i_shift].w0;

                                        i_island_int[my_island_id]++;
                                    }
                                    else{

                                        C2C_index = tmp_C2Cs.island_offsets_in[my_island_id] + i_island_incoming[my_island_id];

                                        C2Cs[_i_C2C].shifts_in[0][C2C_index].c0 =     tmp_C2Cs.shifts[i_shift].c0;
                                        C2Cs[_i_C2C].shifts_in[0][C2C_index].node0 =  tmp_C2Cs.shifts[i_shift].node0;
                                        C2Cs[_i_C2C].shifts_in[0][C2C_index].c1 =     tmp_C2Cs.shifts[i_shift].c1;
                                        C2Cs[_i_C2C].shifts_in[0][C2C_index].node1 =  tmp_C2Cs.shifts[i_shift].node1;
                                        C2Cs[_i_C2C].shifts_in[0][C2C_index].w0 =     tmp_C2Cs.shifts[i_shift].w0;

                                        i_island_incoming[my_island_id]++;
                                    }
                                }

                                C2Cs[_i_C2C].island_offsets[0][0] = 0;
                                C2Cs[_i_C2C].island_offsets_in[0][0] = 0;

                                loop_islands{

                                    C2Cs[_i_C2C].island_offsets[0][(i_island + 1)] = tmp_C2Cs.island_offsets[(i_island + 1)];
                                    C2Cs[_i_C2C].island_offsets_in[0][(i_island + 1)] = tmp_C2Cs.island_offsets_in[(i_island + 1)];
                                }
                            }
                        }

                        /* B.5 free tmp_C2Cs stuff */
                        {
                            free(tmp_C2Cs.shifts);

                            free(tmp_C2Cs.island_offsets);
                            free(tmp_C2Cs.island_offsets_in);
                        }
                    }

                    fclose(fi);
                }
            }
        }

#endif
    }

    /* C. Free locally allocated memory */
    {
#if RP_NODE
        free(i_island_int);
        free(i_island_incoming);
#endif
    }

    /* D. Set-up C2Cs for upper layers */
    {
#if RP_NODE

        short   debug_this_code = 0;

        if(Solver_Dict.number_of_layers > 1){

#if 1       /* local vars */

            int     i_layer, i_upper_shift;

            int     number_of_shifts, number_of_shifts_in_upper_layer, c0, c1, upper_c0, upper_c1;

            double  w0;
#endif

            loop_states{

                loop_phases{

                    loop_frames{

                        i_layer = 0;

                        while(upper_layer < Solver_Dict.number_of_layers){

                            /* D.1 fill shifts of upper_layer */
                            {
                                /* TODO remove doublets */

                                number_of_shifts = C2Cs[_i_C2C].number_of_shifts[i_layer];

                                number_of_shifts_in_upper_layer = 0;

                                loop_shifts{

                                    c0 = C2Cs[_i_C2C].shifts[_i_shift].c0;
                                    c1 = C2Cs[_i_C2C].shifts[_i_shift].c1;

                                    upper_c0 = _C.parent_cell[c0];
                                    upper_c1 = _C.parent_cell[c1];

                                    if(upper_c0 != upper_c1){

                                        number_of_shifts_in_upper_layer++;
                                    }
                                }

                                C2Cs[_i_C2C].number_of_shifts[upper_layer] = number_of_shifts_in_upper_layer;

                                /* TODO account for multiple islands */
                                C2Cs[_i_C2C].island_offsets[upper_layer][0] = 0;
                                C2Cs[_i_C2C].island_offsets[upper_layer][1] = number_of_shifts_in_upper_layer;

                                C2Cs[_i_C2C].shifts[upper_layer] = (C2C_shift_type*)malloc(number_of_shifts_in_upper_layer * sizeof(C2C_shift_type));

                                i_upper_shift = 0;

                                loop_shifts{

                                    c0 = C2Cs[_i_C2C].shifts[_i_shift].c0;
                                    c1 = C2Cs[_i_C2C].shifts[_i_shift].c1;

                                    upper_c0 = _C.parent_cell[c0];
                                    upper_c1 = _C.parent_cell[c1];

                                    w0 = C2Cs[_i_C2C].shifts[_i_shift].w0;

                                    if(upper_c0 != upper_c1){

                                        C2Cs[_i_C2C].shifts[upper_layer][i_upper_shift].node0 = myid;
                                        C2Cs[_i_C2C].shifts[upper_layer][i_upper_shift].c0 =    upper_c0;
                                        C2Cs[_i_C2C].shifts[upper_layer][i_upper_shift].node1 = myid;
                                        C2Cs[_i_C2C].shifts[upper_layer][i_upper_shift].c1 =    upper_c1;
                                        C2Cs[_i_C2C].shifts[upper_layer][i_upper_shift].w0 =    w0;

                                        i_upper_shift++;
                                    }
                                }

                                if((debug_this_code)&&(i_frame == 0)){

                                    Message("\nD.2 myid %d: i_layer %d, number_of_shifts %d/%d", myid, i_layer, number_of_shifts, number_of_shifts_in_upper_layer);
                                }
                            }

                            i_layer++;
                        }
                    }
                }
            }
        }
#endif
    }

    /* E. set up parallel C2C communication */
    {
#if RP_NODE
        if(Solver_Dict.number_of_layers > 1){

            init_upper_layer_parallel_C2Cs();
        }

        init_parallel_C2Cs();
#endif
    }

    /* F. Message & Transcript */
    {
#if RP_NODE

        total_number_of_C2Cs_read = PRF_GISUM1(total_number_of_C2Cs_read);

        if(myid == 0){

            FILE    *f_trn = NULL;
            
            char    file_name[80];
            
            sprintf(file_name,"%s", File_Dict.Run_Transscript_filename);

            f_trn = fopen(file_name, "a");

            if(f_trn){

                fprintf(f_trn,"\n\nrCFD_read_C2Cs");

                fprintf(f_trn,"\n\n   Read %d C2Cs from %s in %f seconds",

                    total_number_of_C2Cs_read, File_Dict.C2C_filename,

                    (double)(clock() - Solver.clock)/(double)CLOCKS_PER_SEC);

                fclose(f_trn);
            }

            Message("\n...rCFD_read_C2Cs ->  Read %d C2Cs from %s in %f seconds",

                total_number_of_C2Cs_read, File_Dict.C2C_filename,

                (double)(clock() - Solver.clock)/(double)CLOCKS_PER_SEC);
        }
#endif
    }
}

/*************************************************************************************/
DEFINE_ON_DEMAND(rCFD_run)
/*************************************************************************************/
{
    /* Contents:

        ST:     Next state

        N+1:    Next frame

        I:      Initializations

        S1:     Source term 1

        BC:     Update BC

        C:      Convection (cell shift, fill holes)

            C1      internal shifts
            C2      MPI shifts
            C3      fill holes

        D:      Diffusion

            D1      Face-swaps for balance correction
            D2      Binarization

        S2:     Source term 2

        B:      Update balances

            B1      Update
            B2      Correction

        PP:     On the fly post-processing

    */

#if 1   /* global in function definitions */

    int         i_run, i_island, i_phase, i_layer, i_state, i_state2, i_rec;

    int         prev_layer, i_rec_max;

#if RP_NODE
    int         i_frame, i_data, i_shift, i_node, i_node2, i_cell, i_face, i_drift, i_dim, i_while;

    int         loop_offset0, loop_offset1;
    int         number_of_fill_loops, number_of_unhit_cells;
    int         c, c0, c1;
    int         i_frame_c0, i_frame_c1;
    int         i_tmp;

    int         i_warning_fill_1, i_warning_drift_1;

    short       balance_error_exists;

    double      w0, data0, vol_flip;
    double      drift_volume, local_drift_exchange, local_mass_c0, local_mass_c1, local_mass, hindering_factor, available_c1_mass;
    double      sum_of_conc, flux_in, flux_out, data_in_mean, data_out_mean, flux_mean;
    double      available_exchange, exchange_ratio;

    FILE        *f_out = NULL;
#else

    double      rand_real;
#endif

#endif

    Solver.clock = clock();

    /* W: initialize i_warning's */
    {
#if RP_NODE

        i_warning_fill_1    = 0;
        i_warning_drift_1   = 0;
#endif
    }

    loop_runs{

        /* S: select state and layer for all phases */
        {
            i_state = 0, i_state2 = 0;

            /* set layer */
            {
                prev_layer = Solver.current_layer;

                i_layer = rCFD_user_set_layer(Solver.current_layer);

                if(i_layer != prev_layer){

                    rCFD_map_from_to_layer(prev_layer, i_layer);

                    Solver.current_layer = i_layer;
                }
            }
        }

        /* N+1: Get next frames[islands] for all phases */
        {

            if(Solver_Dict.recurrence_process_on){

#if RP_HOST     /* init rand. gen */

                if(Solver.global_run_counter == 0){

                    srand(time(0));

                    Message("\n\nInitialized random generator");
                }
#endif

                i_rec_max = 1;

                for(i_rec = 0; i_rec < i_layer; i_rec++){

                    i_rec_max *= 2;
                }

                for(i_rec = 0; i_rec < i_rec_max; i_rec++){

                    if(Rec.frame_in_sequence < Rec.sequence_length){

                        Rec.frame_in_sequence++;

                        /* set next frame (within sequence) */
                        loop_islands{

                            if(Rec.global_frame[i_island] < (Solver_Dict.number_of_frames - 1)){

                                Rec.global_frame[i_island]++;
                            }
                            else{

                                Rec.global_frame[i_island] = Rec.jumps[i_state][i_state2][i_island][(Solver_Dict.number_of_frames-1)];
                            }
                        }
                    }
                    else{ /* new rCFD_seq */

                        Rec.frame_in_sequence = 0;
#if RP_HOST
                        rand_real = (double)rand()/(double)RAND_MAX;

                        Rec.sequence_length = Rec_Dict.min_seq_length + (int)(rand_real*(double)(Rec_Dict.max_seq_length - Rec_Dict.min_seq_length));
#endif
                        host_to_node_int_1(Rec.sequence_length);

                        /* set next frame (new sequence) */
                        loop_islands{

                            Rec.global_frame[i_island] = Rec.jumps[i_state][i_state2][i_island][Rec.global_frame[i_island]];
                        }
                    }
                }
            }
            else{

                loop_islands{

                    Rec.global_frame[i_island] = 0;
                }
            }

            if(Solver_Dict.verbose){

                Message0("\n\nNext Frame: i_run %d, i_frame[0] %d", i_run, Rec.global_frame[0]);
            }
        }

        loop_phases{

            /* I: Initialization (weights, data_shift, data_swap, mass_drift) */
            {
#if RP_NODE
                i_phase = rCFD_user_phase_switch(i_phase);

                loop_cells{

                    _C.weight_after_shift[i_cell] =  0.0;
                    _C.weight_after_swap[i_cell] =   0.0;


                    loop_data{

                        _C.data_shift[_i_data] = 0.0;
                        _C.data_swap[_i_data] =  0.0;
                    }

                }

                /* init Balance.node2node_flux */
                loop_data{

                    if(Balance_Dict[i_phase][i_data].type == per_node_balancing){

                        for(i_node = 0; i_node < (node_last + 1); i_node++){

                            for(i_node2 = 0; i_node2 < (node_last + 1); i_node2++){

                                Balance[i_phase][i_data].node2node_flux[i_node][i_node2] = 0.0;
                                Balance[i_phase][i_data].node2node_data_flux[i_node][i_node2] = 0.0;
                            }
                        }
                    }
                }

                if(Solver_Dict.verbose){

                    Message0("\n\nInit: i_run %d, i_phase %d, i_layer %d", i_run, i_phase, i_layer);
                }
#endif
            }

            /* AD1: Access data before shift */
            {
#if RP_NODE
                rCFD_user_access_data_before_shift(i_phase, i_layer);
#endif
            }

            /* C: Convection, cell-to-cell shifts */
            {
#if RP_NODE
                if(Solver_Dict.data_convection_on){

                    /* C1,2 C2C shifts (local, cross-partitions) */
                    loop_islands{

                        i_frame = Rec.global_frame[i_island];

                        /* C1: local shifts */
                        {
                            loop_offset0    = C2Cs[_i_C2C].island_offsets[i_layer][i_island];
                            loop_offset1    = C2Cs[_i_C2C].island_offsets[i_layer][(i_island + 1)];

                            /*Message0("\n\n loop_offset0 %d, loop_offset1 %d\\", loop_offset0, loop_offset1);*/

                            i_tmp = 0;

                            for(i_shift = loop_offset0; i_shift < loop_offset1; i_shift++){

                                c0 = C2Cs[_i_C2C].shifts[_i_shift].c0;
                                c1 = C2Cs[_i_C2C].shifts[_i_shift].c1;
                                w0 = C2Cs[_i_C2C].shifts[_i_shift].w0;

                                if(w0 > 0.0){

                                    i_tmp++;

                                    loop_data{

                                        data0 = _C.data[i_phase][c0][i_data];

                                        _C.data_shift[i_phase][c1][i_data] =

                                            (w0 * data0 + _C.weight_after_shift[c1] * _C.data_shift[i_phase][c1][i_data])/

                                            (w0 +  _C.weight_after_shift[c1]);
                                    }

                                    _C.weight_after_shift[c1] += w0;
                                }
                            }

                            /*Message0("\n\n i_tmp = %d\n\n", i_tmp);*/
                        }

                        /* C2: cross-partition shifts */
                        {
                            shift_parallel_C2C_data(i_state, i_phase, i_frame, i_island, i_layer);
                        }

                    } /* loop_islands */

                    /* C3: fill holes */
                    {
                        number_of_unhit_cells = 1;

                        number_of_fill_loops = 0;

                        while((number_of_unhit_cells > 0) && (number_of_fill_loops < Solver_Dict.max_fill_loops)){

                            number_of_fill_loops++;

                            /* fill holes by face swaps */
                            loop_faces{

                                c0 = _F.c0[i_face];
                                c1 = _F.c1[i_face];

                                if((_C.weight_after_shift[c0] > 0.)&&(_C.weight_after_shift[c1] == 0.0)){

                                    c = c0; c0 = c1; c1 = c;
                                }

                                if((_C.weight_after_shift[c1] > 0.)&&(_C.weight_after_shift[c0] == 0.0)){

                                    loop_data{

                                        _C.data_swap[i_phase][c0][i_data] =

                                            (_C.weight_after_shift[c1] * _C.data_shift[i_phase][c1][i_data] +

                                             _C.weight_after_swap[c0] * _C.data_swap[i_phase][c0][i_data]) /

                                            (_C.weight_after_shift[c1] + _C.weight_after_swap[c0]);
                                    }

                                    _C.weight_after_swap[c0] += _C.weight_after_shift[c1];
                                }
                            }

                            number_of_unhit_cells = 0;

                            loop_cells{

                                if(_C.weight_after_swap[i_cell] > 0.0){

                                    loop_data{

                                        _C.data_shift[i_phase][i_cell][i_data] = _C.data_swap[i_phase][i_cell][i_data];
                                    }

                                    _C.weight_after_shift[i_cell] = _C.weight_after_swap[i_cell];

                                    _C.weight_after_swap[i_cell] = 0.0;
                                }

                                i_frame = Rec.global_frame[_C.island_id[i_cell]];

                                if((_C.weight_after_shift[i_cell] == 0.0) && (_C.vof[i_frame][i_cell][i_phase] > 0.0)){

                                    number_of_unhit_cells++;
                                }

                            }
                        }

                        /* data = data_update */
                        loop_cells{

                            loop_data{

                                _C.data[_i_data] = _C.data_shift[_i_data];
                            }
                        }

                        i_warning_fill_1 = PRF_GISUM1(number_of_unhit_cells);
                    }
                }

                if(Solver_Dict.verbose){

                    Message0("\n\nC2C shifts: i_run %d, i_phase %d, i_layer %d", i_run, i_phase, i_layer);

                    if(i_warning_fill_1){

                        Message0("\n\nWARNING rCFD_run, c2c shifts, fill holes: %d unvisited cells exist", i_warning_fill_1);
                    }
                }
#endif
            }

            /* D: Diffusion, face swaps */
            {
                /* D.1 Drifting */
                if((Solver_Dict.data_drifting_on) && (Solver.current_layer == 0)){

                    /* Existing Limitations (11/21)

                        L1: drifting pattern might depend on face order

                        L2: drifting of tmperature_data not yet implemented

                        L3: drifting only implemented for i_layer == 0
                    */
#if RP_NODE
                    loop_data{

                        if(Data_Dict[i_phase][i_data].drifting_data){

                            for(i_drift = 0; i_drift < Solver_Dict.number_of_drift_loops; i_drift++){

                                /* init _C.drift_exchange */
                                loop_cells{

                                    _C.drift_exchange[i_cell] = 0.0;
                                }

                                switch (Data_Dict[i_phase][i_data].type){

                                    case generic_data:
                                    {

                                        /* set _C.drift_exchange */
                                        loop_faces{

                                            if(valid_parallel_face){

                                                c0 = _F.c0[i_face];
                                                c1 = _F.c1[i_face];

                                                i_frame_c0 = Rec.global_frame[_C.island_id[c0]];
                                                i_frame_c1 = Rec.global_frame[_C.island_id[c1]];

                                                local_mass_c0 = _C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase] * Phase_Dict[i_phase].density;
                                                local_mass_c1 = _C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase] * Phase_Dict[i_phase].density;

                                                /* don't drift across interfaces */
                                                if(local_mass_c0 * local_mass_c1 > 0.0){

                                                    drift_volume = 0.0;

                                                    loop_dim{

                                                        drift_volume += Data_Dict[i_phase][i_data].drift_velocity[i_dim] * _F.area[i_face][i_dim] * Solver.timestep_width[i_layer] /

                                                            (double)Solver_Dict.number_of_drift_loops;
                                                    }

                                                    /* flip cells, such that flux is from c0 to c1 */
                                                    if(drift_volume < 0.0){

                                                        c = c0; c0 = c1; c1 = c;

                                                        i_frame = i_frame_c0; i_frame_c0 = i_frame_c1; i_frame_c1 = i_frame;

                                                        local_mass = local_mass_c0; local_mass_c0 = local_mass_c1; local_mass_c1 = local_mass;

                                                        drift_volume *= -1.0;
                                                    }

                                                    if(drift_volume > _C.volume[c0]){

                                                        drift_volume = _C.volume[c0];
                                                    }

                                                    local_drift_exchange = 0.0;

                                                    if(_C.volume[c0] > 0.0){

                                                        local_drift_exchange = _C.data[i_phase][c0][i_data] * drift_volume /_C.volume[c0];
                                                    }

                                                    _C.drift_exchange[c0] -= local_drift_exchange;
                                                    _C.drift_exchange[c1] += local_drift_exchange;
                                                }
                                            }
                                        }

                                        /* parallel exchange of local_drift_exchange */
                                        {
                                            sum_up_parallel_corona_cells(_C.drift_exchange, i_layer);
                                        }

                                        /* adjust data by _C.drift_exchange */
                                        loop_cells{

                                            _C.data[i_phase][i_cell][i_data] += _C.drift_exchange[i_cell];
                                        }

                                        break;
                                    }

                                    case concentration_data:
                                    {
                                        /* set _C.drift_exchange */
                                        loop_faces{

                                            if(valid_parallel_face){

                                                c0 = _F.c0[i_face];
                                                c1 = _F.c1[i_face];

                                                i_frame_c0 = Rec.global_frame[_C.island_id[c0]];
                                                i_frame_c1 = Rec.global_frame[_C.island_id[c1]];

                                                local_mass_c0 = _C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase] * Phase_Dict[i_phase].density;
                                                local_mass_c1 = _C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase] * Phase_Dict[i_phase].density;

                                                /* don't drift across interfaces */
                                                if(local_mass_c0 * local_mass_c1 > 0.0){

                                                    drift_volume = 0.0;

                                                    loop_dim{

                                                        drift_volume += Data_Dict[i_phase][i_data].drift_velocity[i_dim] * _F.area[i_face][i_dim] *

                                                            Solver.timestep_width[i_layer] / (double)Solver_Dict.number_of_drift_loops;

                                                        /* drift_volume is defined per face, w/o information of vof's of c0/c1 */
                                                    }

                                                    /* flip cells, such that flux is from c0 to c1 */
                                                    if(drift_volume < 0.0){

                                                        c = c0; c0 = c1; c1 = c;

                                                        i_frame = i_frame_c0; i_frame_c0 = i_frame_c1; i_frame_c1 = i_frame;

                                                        local_mass = local_mass_c0; local_mass_c0 = local_mass_c1; local_mass_c1 = local_mass;

                                                        drift_volume *= -1.0;
                                                    }

                                                    /* in granular flow, drifting ability depend on solid volume fraction */
                                                    if(Phase_Dict[i_phase].hindered_drift_on){

                                                        hindering_factor = 1.0 - _C.vof[i_frame_c0][c0][i_phase] / Phase_Dict[i_phase].vof_max;

                                                        if(hindering_factor < 0.0){

                                                            hindering_factor = 0.0;
                                                        }

                                                        drift_volume *= hindering_factor;
                                                    }

#if 0   /* code sketch changing_c2c_bulk_diameter_on */

                                                    if(Phase_Dict[i_phase].changing_c2c_bulk_diameter_on){

                                                        loop_data_2{

                                                            if((_Data_Dict[i_phase][i_data].fractional_diameter > 0.0)&&(_C.data[i_phase][c0][i_data] > 0.0)){

                                                                c0_mean_diam += _Data_Dict[i_phase][i_data].fractional_diameter * _C.data[i_phase][c0][i_data];

                                                                c0_mass_fraction += _C.data[i_phase][c0][i_data];
                                                            }
                                                        }

                                                        if(c0_mass_fraction > 0.0){

                                                            c0_mean_diam /= c0_mass_fraction;
                                                        }
                                                        else{
                                                            c0_mean_diam = Phase_Dict[i_phase].c2c_bulk_diameter;
                                                        }

                                                        fractional_diameter = _Data_Dict[i_phase][i_data].fractional_diameter;

                                                        if(fractional_diameter != Phase_Dict[i_phase].c2c_bulk_diameter){

                                                            drift_scale = (fractional_diameter - c0_mw_mean_diam)/(fractional_diameter - Phase_Dict[i_phase].c2c_bulk_diameter);
                                                        }
                                                        else{

                                                            drift_scale = 1.0;
                                                        }

                                                        drift_volume *= drift_scale;
                                                    }

                                                    /* end of code sketch */
#endif

                                                    if(drift_volume > _C.volume[c0]){

                                                        drift_volume = _C.volume[c0];

                                                        i_warning_drift_1 ++;
                                                    }

                                                    local_drift_exchange = _C.data[i_phase][c0][i_data] * drift_volume * _C.vof[i_frame_c0][c0][i_phase] *

                                                        Phase_Dict[i_phase].density;     /* [kg] */

                                                    available_c1_mass = (1.0 - _C.data[i_phase][c1][i_data]) * _C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase] *

                                                        Phase_Dict[i_phase].density;

                                                    if(local_drift_exchange > available_c1_mass){

                                                        local_drift_exchange = available_c1_mass;
                                                    }

                                                    _C.drift_exchange[c0] -= local_drift_exchange;
                                                    _C.drift_exchange[c1] += local_drift_exchange;
                                                }
                                            }
                                        }

                                        /* parallel exchange of _C.drift_exchange */
                                        {
                                            sum_up_parallel_corona_cells(_C.drift_exchange, i_layer);
                                        }

                                        /* adapt _C.data by _C.drift_exchange */
                                        loop_cells{

                                            i_frame = Rec.global_frame[_C.island_id[i_cell]];

                                            local_mass = _C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density;

                                            if(local_mass > 1.0e-10){

                                                _C.data[_i_data] = ((_C.data[_i_data] * local_mass) + _C.drift_exchange[i_cell]) / local_mass;
                                            }
                                        }

                                        break;
                                    }

                                    case temperature_data:
                                    {

                                        /* to do: replace local_drift_exchange by drift_exchange */

                                        break;
                                    }

                                    default: break;
                                }
                            }
                        }
                    }
#endif
                }

                /* D.2 Physical diffusion */
                if(Solver_Dict.face_diffusion_on){
#if RP_NODE
                    /* TODO: allow for more diff. loops, account for heterogeneous diffusion */

                    loop_data{

                        switch (Data_Dict[i_phase][i_data].type){

                            case concentration_data:
                            {
                                /* Init data_swap */
                                loop_int_cells{

                                    _C.data_swap[_i_data] = 0.0;
                                }

                                /* define swap masses */
                                loop_faces{

                                    c0 = _F.c0[i_face];
                                    c1 = _F.c1[i_face];

                                    if(_C.data[_c0_data] > _C.data[_c1_data]){

                                        c = c0; c0 = c1; c1 = c;
                                    }

                                    if(_C.data[_c1_data] > _C.data[_c0_data]){

                                        i_frame_c0 = Rec.global_frame[_C.island_id[c0]];
                                        i_frame_c1 = Rec.global_frame[_C.island_id[c1]];

                                        if((_C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase]) < (_C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase])){

                                            vol_flip = _C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase];
                                        }
                                        else{
                                            vol_flip = _C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase];
                                        }

                                        _C.data_swap[i_phase][c1][i_data] -= vol_flip *

                                            (_C.data[i_phase][c1][i_data] - _C.data[i_phase][c0][i_data]) / 2. * Data_Dict[i_phase][i_data].physical_diff;

                                        _C.data_swap[i_phase][c0][i_data] += vol_flip *

                                            (_C.data[i_phase][c1][i_data] - _C.data[i_phase][c0][i_data]) / 2. * Data_Dict[i_phase][i_data].physical_diff;
                                    }
                                }

                                /* update data by data_swap */
                                loop_int_cells{

                                    i_frame = Rec.global_frame[_C.island_id[i_cell]];

                                    if((_C.volume[i_cell] * _C.vof[_i_vof]) > 0.0){

                                        _C.data[_i_data] += _C.data_swap[_i_data] / (_C.volume[i_cell] * _C.vof[_i_vof]);
                                    }

                                    _C.data_swap[_i_data] = 0.0;
                                }

                                break;
                            }

                        }

                    }
#endif
                }

                /* D.3 Binarization */
                if(Solver_Dict.data_binarization_on){
#if RP_NODE
                    /* data as outer loop, because it is unlikely that many data get binarized */
                    loop_data{

                        if(Data_Dict[i_phase][i_data].type == binary_data){

                            /* Init data_swap */
                            loop_cells{

                                _C.data_swap[i_phase][i_cell][i_data] = 0.0;
                            }

                            /* Artificial Diffusion (mass neutral) */
                            loop_faces{

                                c0 = _F.c0[i_face];
                                c1 = _F.c1[i_face];

                                if(_C.data[i_phase][c0][i_data] > _C.data[i_phase][c1][i_data]){

                                    c = c0; c0 = c1; c1 = c;
                                }

                                if(_C.data[i_phase][c1][i_data] > _C.data[i_phase][c0][i_data]){

                                    if(_C.volume[c0] < _C.volume[c1]){

                                        vol_flip = _C.volume[c0] * Data_Dict[i_phase][i_data].binarization_art_diff;
                                    }
                                    else{
                                        vol_flip = _C.volume[c1] * Data_Dict[i_phase][i_data].binarization_art_diff;
                                    }

                                    _C.data_swap[i_phase][c1][i_data] -= vol_flip *

                                        (_C.data[i_phase][c1][i_data] - _C.data[i_phase][c0][i_data]) / 2.;

                                    _C.data_swap[i_phase][c0][i_data] += vol_flip *

                                        (_C.data[i_phase][c1][i_data] - _C.data[i_phase][c0][i_data]) / 2.;
                                }
                            }

                            /* Binarization (mass acting) */
                            loop_cells{

                                if(_C.data[i_phase][i_cell][i_data] > 0.5){

                                    _C.data[i_phase][i_cell][i_data] = 1.0;
                                }
                                else{

                                    _C.data[i_phase][i_cell][i_data] = 0.0;
                                }
                            }
                        }
                    }
#endif
                }

                if(Solver_Dict.verbose){
#if RP_NODE
                    Message0("\n\nFace swaps: i_run %d, i_phase %d, i_layer %d", i_run, i_phase, i_layer);

                    i_warning_drift_1 = PRF_GISUM1(i_warning_drift_1);

                    if(i_warning_drift_1){

                        Message0("\n\nWARNING rCFD_run, face swap, data drifting: limited drift by c0-volume at %d faces", i_warning_drift_1);
                    }
#endif
                }
            }

            /* AD2: Access data after swap */
            {
#if RP_NODE
                rCFD_user_access_data_after_swap(i_phase, i_layer);
#endif
            }

            /* Balance correction */
            {

                /* Update balances */
                {
                    /* B1: mass_integral, mass_integral_global */
                    {
#if RP_NODE
                        loop_data{

                            Balance[_i_balance].mass_integral = 0.0;
                        }

                        loop_int_cells{

                            i_frame = Rec.global_frame[_C.island_id[i_cell]];

                            loop_data{

                                switch (Data_Dict[i_phase][i_data].type){

                                    case temperature_data:

                                        Balance[_i_balance].mass_integral += _C.data[_i_data] *

                                            _C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density * Phase_Dict[i_phase].heat_capacity;

                                        break;

                                    case concentration_data:

                                        Balance[_i_balance].mass_integral += _C.data[_i_data] *

                                            _C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density;

                                        break;

                                    case generic_data:

                                        Balance[i_phase][i_data].mass_integral += _C.data[i_phase][i_cell][i_data];

                                        break;

                                    default: break;
                                }
                            }
                        }

                        loop_data{

                            Balance[i_phase][i_data].mass_integral_global = PRF_GRSUM1(Balance[i_phase][i_data].mass_integral);
                        }
#endif
                    }

                    /* B2: data_int_target += sources */
                    {
#if RP_NODE
                        loop_data{

                            Balance[_i_balance].mass_integral_target += Balance[_i_balance].mass_source;

                            Balance[_i_balance].mass_source_global = PRF_GRSUM1(Balance[_i_balance].mass_source);

                            Balance[_i_balance].mass_integral_target_global += Balance[_i_balance].mass_source_global;

                            Balance[_i_balance].mass_source = 0.0;

                            Balance[_i_balance].mass_source_global = 0.0;
                        }
#endif
                    }

                    /* B3: mass_integral_target += node2node fluxes */
                    if(0==1){
#if RP_NODE
                        loop_data{

                            if(Balance_Dict[i_phase][i_data].type == per_node_balancing){

                                /* ToDo efficiency of multiple MPI comm ? */
                                for(i_node = 0; i_node < (node_last + 1); i_node++){

                                    for(i_node2 = 0; i_node2 < (node_last + 1); i_node2++){

                                        Balance[i_phase][i_data].node2node_flux[i_node][i_node2] =

                                            PRF_GRSUM1(Balance[i_phase][i_data].node2node_flux[i_node][i_node2]);

                                        Balance[i_phase][i_data].node2node_data_flux[i_node][i_node2] =

                                            PRF_GRSUM1(Balance[i_phase][i_data].node2node_data_flux[i_node][i_node2]);
                                    }
                                }

                                flux_in = 0.0; flux_out = 0.0;

                                data_in_mean = 0.0; data_out_mean = 0.0;

                                for(i_node = 0; i_node < (node_last + 1); i_node++){

                                    flux_in += Balance[i_phase][i_data].node2node_flux[i_node][myid];
                                    data_in_mean += Balance[i_phase][i_data].node2node_data_flux[i_node][myid];

                                    flux_out += Balance[i_phase][i_data].node2node_flux[myid][i_node];
                                    data_out_mean += Balance[i_phase][i_data].node2node_data_flux[myid][i_node];
                                }

                                if(flux_in > 0.0){

                                    data_in_mean /= flux_in;
                                }

                                if(flux_out > 0.0){

                                    data_out_mean /= flux_out;
                                }


                                flux_mean = (flux_in + flux_out) / 2.;

                                Balance[i_phase][i_data].mass_integral_target += flux_mean * (data_in_mean - data_out_mean);

                                /* ToDo: check that sum of local balances fulfill global balance */
                            }
                        }
#endif
                    }

                    /* B4: data_int_error  */
                    {
#if RP_NODE
                        loop_data{

                            Balance[_i_balance].mass_error =

                                Balance[_i_balance].mass_integral_target - Balance[_i_balance].mass_integral;

                            Balance[_i_balance].mass_error_global = PRF_GRSUM1(Balance[_i_balance].mass_error);
                        }
#endif
                    }

                }

                /* Balance correction by face swaps */
                if(Solver_Dict.balance_correction_on){

                    /* existing limitations:

                        L1: only implemented for generic_data, concentration_data and temperature_data

                        L2: only implemented/tested for global balancing

                        L3: balancing of concentration_data only works for constant phase properties (density, heat capacity)
                    */
#if RP_NODE
                    if((i_run % Solver_Dict.balance_correction_update) == 0){

                        loop_data{

                            switch (Data_Dict[i_phase][i_data].type){

                                case generic_data:
                                {
                                    i_while = 0;

                                    balance_error_exists = 1;

                                    while ((balance_error_exists) && (i_while < Balance_Dict[i_phase][i_data].max_correction_loops)){

                                        loop_int_cells{

                                            _C.data_swap[i_phase][i_cell][i_data] = 0.0;
                                        }

                                        available_exchange = 0.0;

                                        loop_int_faces{

                                            c0 = _F.c0[i_face];
                                            c1 = _F.c1[i_face];

                                            available_exchange += 0.5 * fabs(_C.data[i_phase][c1][i_data] - _C.data[i_phase][c0][i_data]) * Solver_Dict.face_swap_max_per_loop;
                                        }

                                        if(Balance_Dict[i_phase][i_data].type == global_balancing){

                                            available_exchange = PRF_GRSUM1(available_exchange);
                                        }

                                        if(available_exchange > fabs(Balance[i_phase][i_data].mass_error_global)){

                                            if(available_exchange > 0.0){

                                                exchange_ratio = fabs(Balance[i_phase][i_data].mass_error_global) / available_exchange;
                                            }
                                            else{
                                                exchange_ratio = 0.0;
                                            }
                                        }
                                        else{

                                            exchange_ratio = 1.0;
                                        }

                                        loop_int_faces{

                                            c0 = _F.c0[i_face];
                                            c1 = _F.c1[i_face];

                                            if(_C.data[i_phase][c0][i_data] > _C.data[i_phase][c1][i_data]){

                                                c = c0, c1 = c0, c0 = c;
                                            }

                                            if(_C.data[i_phase][c0][i_data] < _C.data[i_phase][c1][i_data]){

                                                if(Balance[i_phase][i_data].mass_error_global > 0.0){

                                                    _C.data_swap[i_phase][c0][i_data] += 0.5 * exchange_ratio *

                                                        (_C.data[i_phase][c1][i_data] - _C.data[i_phase][c0][i_data]) * Solver_Dict.face_swap_max_per_loop;
                                                }

                                                if(Balance[i_phase][i_data].mass_error_global < 0.0){

                                                    _C.data_swap[i_phase][c1][i_data] -= 0.5 * exchange_ratio *

                                                        (_C.data[i_phase][c1][i_data] - _C.data[i_phase][c0][i_data]) * Solver_Dict.face_swap_max_per_loop;
                                                }
                                            }
                                        }

                                        Balance[i_phase][i_data].mass_integral = 0.0;

                                        loop_int_cells{

                                            _C.data[i_phase][i_cell][i_data] += _C.data_swap[i_phase][i_cell][i_data];

                                            Balance[i_phase][i_data].mass_integral += _C.data[i_phase][i_cell][i_data];
                                        }

                                        Balance[i_phase][i_data].mass_integral_global = PRF_GRSUM1(Balance[i_phase][i_data].mass_integral);

                                        Balance[i_phase][i_data].mass_error_global = Balance[i_phase][i_data].mass_integral_target_global - Balance[i_phase][i_data].mass_integral_global;

                                        if(fabs(Balance[i_phase][i_data].mass_error_global) < (Balance_Dict[i_phase][i_data].accuracy_level * Balance[i_phase][i_data].mass_integral_target_global)){

                                            balance_error_exists = 0;
                                        }

                                        i_while++;
                                    }

                                    if(i_while == Balance_Dict[i_phase][i_data].max_correction_loops){

                                        Message0("\nWARNING: Balance for i_phase %d i_data %d beyond accuracy level %f", i_phase, i_data, Balance_Dict[i_phase][i_data].accuracy_level);
                                    }

                                    break;
                                }

                                case concentration_data:
                                {
                                    i_while = 0;

                                    balance_error_exists = 1;

                                    while ((balance_error_exists) && (i_while < Balance_Dict[i_phase][i_data].max_correction_loops)){

                                        loop_int_cells{

                                            _C.data_swap[_i_data] = 0.0;
                                        }

                                        available_exchange = 0.0;

                                        loop_int_faces{

                                            c0 = _F.c0[i_face];
                                            c1 = _F.c1[i_face];

                                            if(_C.data[_c0_data] > _C.data[_c1_data]){

                                                c = c0; c0 = c1; c1 = c;
                                            }

                                            if(_C.data[_c1_data] > _C.data[_c0_data]){

                                                i_frame_c0 = Rec.global_frame[_C.island_id[c0]];
                                                i_frame_c1 = Rec.global_frame[_C.island_id[c1]];

                                                if((_C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase]) < (_C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase])){

                                                    vol_flip = _C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase];
                                                }
                                                else{
                                                    vol_flip = _C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase];
                                                }

                                                available_exchange += 0.5 * (_C.data[_c1_data] - _C.data[_c0_data]) *

                                                    vol_flip * Phase_Dict[i_phase].density * Solver_Dict.face_swap_max_per_loop;
                                            }
                                        }

                                        if(Balance_Dict[i_phase][i_data].type == global_balancing){

                                            available_exchange = PRF_GRSUM1(available_exchange);

                                            if(available_exchange > fabs(Balance[i_phase][i_data].mass_error_global)){

                                                if(available_exchange > 0.0){

                                                    exchange_ratio = fabs(Balance[i_phase][i_data].mass_error_global) /

                                                        available_exchange;
                                                }
                                                else{
                                                    exchange_ratio = 0.0;
                                                }
                                            }
                                            else{

                                                exchange_ratio = 1.0;
                                            }
                                        }
                                        else{

                                            exchange_ratio = 0.0;
                                        }

                                        loop_int_faces{

                                            c0 = _F.c0[i_face];
                                            c1 = _F.c1[i_face];

                                            if(_C.data[_c0_data] > _C.data[_c1_data]){

                                                c = c0, c1 = c0, c0 = c;
                                            }

                                            if(_C.data[_c0_data] < _C.data[_c1_data]){

                                                i_frame_c0 = Rec.global_frame[_C.island_id[c0]];
                                                i_frame_c1 = Rec.global_frame[_C.island_id[c1]];

                                                if((_C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase]) < (_C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase])){

                                                    vol_flip = _C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase];
                                                }
                                                else{
                                                    vol_flip = _C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase];
                                                }

                                                if(Balance[_i_balance].mass_error_global > 0.0){

                                                    _C.data_swap[_c0_data] += exchange_ratio * (_C.data[_c1_data] - _C.data[_c0_data]) *

                                                        vol_flip * Phase_Dict[i_phase].density * Solver_Dict.face_swap_max_per_loop;
                                                }

                                                if(Balance[_i_balance].mass_error_global < 0.0){

                                                    _C.data_swap[_c1_data] -= exchange_ratio * (_C.data[_c1_data] - _C.data[_c0_data]) *

                                                        vol_flip * Phase_Dict[i_phase].density * Solver_Dict.face_swap_max_per_loop;
                                                }
                                            }
                                        }

                                        Balance[_i_balance].mass_integral = 0.0;

                                        loop_int_cells{

                                            i_frame = Rec.global_frame[_C.island_id[i_cell]];

                                            if((_C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density) > 0.0){

                                                _C.data[_i_data] = (_C.data[_i_data] * _C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density + _C.data_swap[_i_data]) /

                                                    (_C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density);
                                            }

                                            Balance[_i_balance].mass_integral += _C.data[_i_data] *

                                                _C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density;
                                        }

                                        Balance[i_phase][i_data].mass_integral_global = PRF_GRSUM1(Balance[i_phase][i_data].mass_integral);

                                        Balance[i_phase][i_data].mass_error_global = Balance[i_phase][i_data].mass_integral_target_global - Balance[i_phase][i_data].mass_integral_global;

                                        if(fabs(Balance[i_phase][i_data].mass_error_global) <

                                            (Balance_Dict[i_phase][i_data].accuracy_level * Balance[i_phase][i_data].mass_integral_target_global)){

                                            balance_error_exists = 0;
                                        }

                                        i_while++;
                                    }

                                    if(i_while == Balance_Dict[i_phase][i_data].max_correction_loops){

                                        Message0("\nWARNING: Balance for i_phase %d i_data %d beyond accuracy level %f", i_phase, i_data, Balance_Dict[i_phase][i_data].accuracy_level);
                                    }

                                    break;
                                }

                                case temperature_data:
                                {
                                    i_while = 0;

                                    balance_error_exists = 1;

                                    while ((balance_error_exists) && (i_while < Balance_Dict[i_phase][i_data].max_correction_loops)){

                                        loop_int_cells{

                                            _C.data_swap[_i_data] = 0.0;
                                        }

                                        available_exchange = 0.0;

                                        loop_int_faces{

                                            c0 = _F.c0[i_face];
                                            c1 = _F.c1[i_face];

                                            if(_C.data[_c0_data] > _C.data[_c1_data]){

                                                c = c0; c0 = c1; c1 = c;
                                            }

                                            if(_C.data[_c1_data] > _C.data[_c0_data]){

                                                i_frame_c0 = Rec.global_frame[_C.island_id[c0]];
                                                i_frame_c1 = Rec.global_frame[_C.island_id[c1]];

                                                if((_C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase]) < (_C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase])){

                                                    vol_flip = _C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase];
                                                }
                                                else{
                                                    vol_flip = _C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase];
                                                }

                                                available_exchange += 0.5 * (_C.data[_c1_data] - _C.data[_c0_data]) *

                                                    vol_flip * Phase_Dict[i_phase].density * Phase_Dict[i_phase].heat_capacity *

                                                    Solver_Dict.face_swap_max_per_loop;
                                            }
                                        }

                                        if(Balance_Dict[i_phase][i_data].type == global_balancing){

                                            available_exchange = PRF_GRSUM1(available_exchange);

                                            if(available_exchange > fabs(Balance[i_phase][i_data].mass_error_global)){

                                                if(available_exchange > 0.0){

                                                    exchange_ratio = fabs(Balance[i_phase][i_data].mass_error_global) /

                                                        available_exchange;
                                                }
                                                else{
                                                    exchange_ratio = 0.0;
                                                }
                                            }
                                            else{

                                                exchange_ratio = 1.0;
                                            }
                                        }
                                        else{

                                            exchange_ratio = 0.0;
                                        }

                                        loop_int_faces{

                                            c0 = _F.c0[i_face];
                                            c1 = _F.c1[i_face];

                                            if(_C.data[_c0_data] > _C.data[_c1_data]){

                                                c = c0, c1 = c0, c0 = c;
                                            }

                                            if(_C.data[_c0_data] < _C.data[_c1_data]){

                                                i_frame_c0 = Rec.global_frame[_C.island_id[c0]];
                                                i_frame_c1 = Rec.global_frame[_C.island_id[c1]];

                                                if((_C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase]) < (_C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase])){

                                                    vol_flip = _C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase];
                                                }
                                                else{
                                                    vol_flip = _C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase];
                                                }

                                                if(Balance[i_phase][i_data].mass_error_global > 0.0){

                                                    _C.data_swap[_c0_data] += exchange_ratio * (_C.data[_c1_data] - _C.data[_c0_data]) *

                                                        vol_flip * Phase_Dict[i_phase].density * Phase_Dict[i_phase].heat_capacity *

                                                        Solver_Dict.face_swap_max_per_loop;
                                                }

                                                if(Balance[i_phase][i_data].mass_error_global < 0.0){

                                                    _C.data_swap[_c1_data] -= exchange_ratio * (_C.data[_c1_data] - _C.data[_c0_data]) *

                                                        vol_flip * Phase_Dict[i_phase].density * Phase_Dict[i_phase].heat_capacity *

                                                        Solver_Dict.face_swap_max_per_loop;
                                                }
                                            }
                                        }

                                        Balance[i_phase][i_data].mass_integral = 0.0;

                                        loop_int_cells{

                                            i_frame = Rec.global_frame[_C.island_id[i_cell]];

                                            if((_C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density * Phase_Dict[i_phase].heat_capacity) > 0.0){

                                                _C.data[_i_data] = (_C.data[_i_data] * _C.volume[i_cell] * _C.vof[i_frame][i_cell][i_phase] *

                                                    Phase_Dict[i_phase].density * Phase_Dict[i_phase].heat_capacity + _C.data_swap[_i_data]) /

                                                    (_C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density * Phase_Dict[i_phase].heat_capacity);
                                            }

                                            Balance[_i_balance].mass_integral += _C.data[_i_data] *

                                                _C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density * Phase_Dict[i_phase].heat_capacity;
                                        }

                                        Balance[_i_balance].mass_integral_global = PRF_GRSUM1(Balance[_i_balance].mass_integral);

                                        Balance[_i_balance].mass_error_global = Balance[_i_balance].mass_integral_target_global -

                                            Balance[_i_balance].mass_integral_global;

                                        if(fabs(Balance[_i_balance].mass_error_global) <

                                            (Balance_Dict[_i_balance].accuracy_level * Balance[_i_balance].mass_integral_target_global)){

                                            balance_error_exists = 0;
                                        }

                                        i_while++;
                                    }

                                    if(i_while == Balance_Dict[_i_balance].max_correction_loops){

                                        Message0("\nWARNING: Balance for i_phase %d i_data %d beyond accuracy level %f",

                                            i_phase, i_data, Balance_Dict[_i_balance].accuracy_level);
                                    }

                                    break;
                                }
                            }
                        }
                    }
#endif
                }

                /* Write global Balance (node-0) */
                {
#if RP_NODE
                    if(myid == 0){

                        loop_data{

                            if( (Balance_Dict[i_phase][i_data].write_balance_to_file) &&
                                ((Solver.global_run_counter % Balance_Dict[i_phase][i_data].write_balance_to_file_interval) == 0)){

                                if(Solver.balance_file_opened == 0){

                                    f_out = fopen(File_Dict.Balance_filename, "w");

                                    Solver.balance_file_opened = 1;
                                }
                                else{
                                    f_out = fopen(File_Dict.Balance_filename, "a");
                                }

                                if(f_out){

                                    fprintf(f_out, "%d, %d, %d, %e, %e\n", Solver.global_run_counter, i_phase, i_data,

                                        Balance[i_phase][i_data].mass_integral_global, Balance[i_phase][i_data].mass_integral_target_global);

                                    fclose(f_out);
                                }
                            }
                        }
                    }
#endif
                }

                /* Adjust conc. data, such that sum(conc) = 1 */
                if(Solver_Dict.control_conc_sum_on){

#if RP_NODE
                    loop_cells{

                        sum_of_conc = 0.0;

                        loop_data{

                            if(Data_Dict[i_phase][i_data].type == concentration_data){

                                sum_of_conc += _C.data[_i_data];
                            }
                        }

                        loop_data{

                            if(Data_Dict[i_phase][i_data].type == concentration_data){

                                if(sum_of_conc > 0.0){

                                    _C.data[_i_data] /= sum_of_conc;
                                }
                            }
                        }
                    }
#endif
                }

                if(Solver_Dict.verbose){

                    Message0("\n\nBalance: i_run %d, i_phase %d, i_layer %d", i_run, i_phase, i_layer);
                }
            }

        }   /* loop_phases */

        Solver.global_run_counter++;

        Solver.global_time += Solver.timestep_width[i_layer];

        /* Post processing */
        {
#if RP_NODE

            rCFD_user_post();
#endif
        }

    }   /* loop_runs */

    /* D. Message & Transcript */
    {
#if RP_NODE
        if(myid == 0){

            FILE    *f_trn = NULL;
            
            char    file_name[80];
            
            sprintf(file_name,"%s", File_Dict.Run_Transscript_filename);

            f_trn = fopen(file_name, "a");

            if(f_trn){

                fprintf(f_trn,"\n\nrCFD_run");

                fprintf(f_trn,"\n\n   %d Runs @ global run counter %d",

                    Solver_Dict.number_of_runs, Solver.global_run_counter);

                if(i_warning_fill_1){

                    fprintf(f_trn,"\n\n   WARNING: %d fill cell warnings exist, consider increasing fill loops", i_warning_fill_1);
                }

                if(i_warning_drift_1){

                    fprintf(f_trn,"\n\n   WARNING: %d drift warnings exist, consider increasing number of drift loops", i_warning_drift_1);
                }

                fprintf(f_trn,"\n\n   %f seconds real world time took %f seconds compute time",

                    Solver_Dict.global_time_step * (double)Solver_Dict.number_of_runs,

                    (double)(clock() - Solver.clock)/(double)CLOCKS_PER_SEC);

                fclose(f_trn);
            }
            
            Message("\n...rCFD_run -> %d Runs in %f seconds @ global run counter %d and global time %f\n",

                Solver_Dict.number_of_runs, (double)(clock() - Solver.clock)/(double)CLOCKS_PER_SEC,

                Solver.global_run_counter, Solver.global_time);
        }
#endif
    }
}

/*************************************************************************************/
DEFINE_ON_DEMAND(rCFD_free_all)
/*************************************************************************************/
{
    free_all();

#if RP_NODE
    free_all_parallel();
#endif

    /* Transcript and Message */
    if(myid == 0){
#if RP_NODE 
        FILE    *f_trn = NULL;
        
        char    file_name[80];
        
        if(Solver_Dict.mode == preparation_mode){
            
            sprintf(file_name,"%s", File_Dict.Prep_Transscript_filename);
        }
        else{
            
            sprintf(file_name,"%s", File_Dict.Run_Transscript_filename);
        }
                    
        f_trn = fopen(file_name, "a" );

        if(f_trn){

            time_t      current_time = time(NULL);

            char        *c_time_string = ctime(&current_time);
            
            if(Solver_Dict.mode == preparation_mode){
                
                fprintf(f_trn,"\n\nrCFD_prep ends @ %s", c_time_string);
            }
            else{
                
                fprintf(f_trn,"\n\nrCFD_run ends @ %s", c_time_string);
            }

            fclose(f_trn);
        }
        
        Message0("\n\n...rCFD_free_all");
#endif      
    }
}
