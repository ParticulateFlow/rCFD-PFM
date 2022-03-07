#ifndef RCFD_MACROS
#define RCFD_MACROS

#include "udf.h"

/* (C)  2021-22
    Stefan Pirker
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/

#if 1   /* global definitions */

    #define     Transcript      ((myid == 0) && (Solver_Dict.verbose == 1))

#endif

#if 1   /* loops */

    #define     loop_cells                      for(i_cell = 0; i_cell < _Cell_Dict.number_of_cells; i_cell++)
    #define     loop_int_cells                  for(i_cell = 0; i_cell < _Cell_Dict.number_of_int_cells; i_cell++)
    #define     loop_ext_cells                  for(i_cell = _Cell_Dict.number_of_int_cells; i_cell < _Cell_Dict.number_of_cells; i_cell++)
    #define     loop_cells_of_upper_layer       for(i_cell = 0; i_cell < Topo_Dict.Cell_Dict[(i_layer + 1)].number_of_cells; i_cell++)
    #define     loop_cells_of_lower_layer       for(i_cell = 0; i_cell < Topo_Dict.Cell_Dict[(i_layer - 1)].number_of_cells; i_cell++)

    #define     loop_children                   for(i_child = 0; i_child < Topo.Cell[i_layer].number_of_children[i_cell]; i_child++)

    #define     loop_data                       for(i_data = 0; i_data < Phase_Dict[i_phase].number_of_data; i_data++)
    #define     loop_data_user                  for(i_user = 0; i_user < Data_Dict[i_phase][i_data].number_of_user_vars; i_user++)

    #define     loop_dim                        for(i_dim = 0; i_dim < 3; i_dim++)

    #define     loop_faces                      for(i_face = 0; i_face < _Face_Dict.number_of_faces; i_face++)
    #define     loop_int_faces                  for(i_face = 0; i_face < _Face_Dict.number_of_int_faces; i_face++)
    #define     loop_faces_of_upper_layer       for(i_face = 0; i_face < Topo_Dict.Face_Dict[(i_layer + 1)].number_of_faces; i_face++)

    #define     loop_frames                     for(i_frame = 0; i_frame < Solver_Dict.number_of_frames; i_frame++)
    #define     loop_frames2                    for(i_frame2 = 0; i_frame2 < Solver_Dict.number_of_frames; i_frame2++)

    #define     loop_islands                    for(i_island = 0; i_island < Solver_Dict.number_of_islands; i_island++)

    #define     loop_layers                     for(i_layer = 0; i_layer < Solver_Dict.number_of_layers; i_layer++)
    #define     loop_layers_but_L0              for(i_layer = 1; i_layer < Solver_Dict.number_of_layers; i_layer++)

    #define     loop_max_swap_loops             for(i_swap = 0; i_swap < max_swap_loops; i_swap ++)

    #define     loop_norms                      for(i_norm = 0; i_norm < Norms.number_of_norms; i_norm++)

    #define     loop_phases                     for(i_phase = 0; i_phase < Solver_Dict.number_of_phases; i_phase++)
    #define     loop_phase_user                 for(i_user = 0; i_user < Phase_Dict[i_phase].number_of_user_vars; i_user++)

    #define     loop_runs                       for(i_run = 0; i_run < Solver_Dict.number_of_runs; i_run++)

    #define     loop_shifts                     for(i_shift = 0; i_shift < number_of_shifts; i_shift++)

    #define     loop_states                     for(i_state = 0; i_state < Solver_Dict.number_of_states; i_state++)
    #define     loop_states2                    for(i_state2 = 0; i_state2 < Solver_Dict.number_of_states; i_state2++)

#endif

#if 1   /* indices */

    #define     _Cell_Dict                      Topo_Dict.Cell_Dict[i_layer]
    #define     _Face_Dict                      Topo_Dict.Face_Dict[i_layer]

    #define     _C                              Topo.Cell[i_layer]
    #define     _F                              Topo.Face[i_layer]

    #define     _i_balance                      i_phase][i_data

    #define     _i_C2C                          i_state][i_phase][i_frame

    #define     _i_data                         i_phase][i_cell][i_data
    #define     _c0_data                        i_phase][c0][i_data
    #define     _c1_data                        i_phase][c1][i_data

    #define     _i_shift                        i_layer][i_shift

    #define     _i_vof                          i_frame][i_cell][i_phase

    #define     upper_layer                     (i_layer + 1)
    #define     lower_layer                     (i_layer - 1)

#endif

#if 1   /* macros for rCFD_prep (mostly Fluent - specific) */

    /* rCFD_write_Tracer_Positions */

    #define     no_ROI_defined                  (Tracer_Dict.region_of_interest_exists == 0)

    #define     x_in_ROI                        ((x[0] > Tracer_Dict.ROI_x_min) && (x[0] < Tracer_Dict.ROI_x_max))
    #define     y_in_ROI                        ((x[1] > Tracer_Dict.ROI_y_min) && (x[1] < Tracer_Dict.ROI_y_max))
    #define     z_in_ROI                        ((x[2] > Tracer_Dict.ROI_z_min) && (x[2] < Tracer_Dict.ROI_z_max))

    #define     cell_in_ROI                     ((x_in_ROI) && (y_in_ROI) && (z_in_ROI))

    /* rCFD_update_Tracers */

    #define     Tracer_just_started                 (p->user[p_just_started] > 0.0)
    #define     not_Tracer_start_interval           (((N_TIME + (int)(Solver_Dict.time_steps_per_monitoring_interval/2)) % Solver_Dict.time_steps_per_monitoring_interval) != 0)
    #define     too_early_to_start_Tracers          (CURRENT_TIME < Solver_Dict.start_time_for_monitoring)
    #define     Tracer_Database_full                (Tracer.frame_counter >= Solver_Dict.number_of_frames)
    #define     Tracer_Database_not_full            !(Tracer_Database_full)
    #define     kill_Tracer_by_coarse_graining      (rand_real > (1./(double)Tracer_Dict.coarse_graining))
    #define     excess_Tracer_cell_crossing_time    (_C.crossing_time[i_phase][i_cell] > (2. * Solver_Dict.global_time_step))
    #define     Tracer_has_crossed_cell_border      ((int)p->user[p_c_old] != i_cell)
    #define     Tracer_from_slow_cell               (p->user[p_time_ratio] > 1.0)
    #define     Tracer_not_stored_yet               (p->user[p_just_killed] < 1.0)
    #define     Tracer_lifetime_expired             ((CURRENT_TIME - p->user[p_start_time]) >= Solver_Dict.global_time_step)

    /* rCFD_write_Norms */

    #define     Norms_write_interval                ((N_TIME % Solver_Dict.time_steps_per_monitoring_interval) == 0)
    #define     Norm_Database_full                  (Norms.frame_counter >= Solver_Dict.number_of_frames)
    #define     Norm_Database_not_full              !(Norm_Database_full)


    /* Fluent Particle Vars */

    #define     p_vars_used_by_others           0
    #define     p_just_started                  (0 + p_vars_used_by_others)
    #define     p_start_time                    (1 + p_vars_used_by_others)
    #define     p_just_killed                   (2 + p_vars_used_by_others)
    #define     p_phase_id                      (3 + p_vars_used_by_others)
    #define     p_phase_fraction                (4 + p_vars_used_by_others)
    #define     p_c0                            (5 + p_vars_used_by_others)
    #define     p_node0                         (6 + p_vars_used_by_others)
    #define     p_w0                            (7 + p_vars_used_by_others)
    #define     p_c_old                         (8 + p_vars_used_by_others)
    #define     p_time_ratio                    (9 + p_vars_used_by_others)
    #define     p_u_rwm                         (10 + p_vars_used_by_others)
    #define     p_v_rwm                         (11 + p_vars_used_by_others)
    #define     p_w_rwm                         (12 + p_vars_used_by_others)
    #define     p_vel_rwm_old                   (13 + p_vars_used_by_others)
#endif

#endif
