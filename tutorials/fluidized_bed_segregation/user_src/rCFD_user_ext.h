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

        static short    rCFD_monitors_initiated = 0;
#endif

        /* names of phases */
        enum{
            ph_gas,
            ph_solid,
            number_of_phase_names
        };

        /* names of data - phase_0 */
        enum{
            number_of_phase0_data_names
        };

        /* names of data - phase_1 */
        enum{
            c_drift_y_pos,  /* small particles */
            c_drift_y_neg,
            c_drift_0,
            number_of_phase1_data_names
        };

        /* names of ANSYS/fluent UDMIs */
        enum{
            c_vof_small,
            c_vof_large,
            c_vof_mix,
            c_vof,
            c_d23,
            c_vof_small_ref,            /* 5    */
            c_vof_large_ref,
            c_d23_ref,
            c_vof_average,              /* 8    */
            c_vof_small_average,
            c_vof_large_average,
            c_vof_small_ref_average,
            c_vof_large_ref_average,
            number_of_UDMI_names
        };

#if RP_NODE

    static  int average_counter = 0;
#endif

    #define     diam_large  4.0e-3
    #define     diam_small  2.0e-3

    #define     _batch_run  0

#endif

    /* user set Dicts */
    /*************************************************************************************/

#if 1

    void rCFD_user_set_Solver_Dict(void)
    {
        Solver_Dict.verbose =                              high_verbose;

        Solver_Dict.number_of_phases =                     2;

        Solver_Dict.number_of_layers =                     1;

        Solver_Dict.max_number_of_cells_per_time_step =    5;

        Solver_Dict.global_time_step =                     2.0 * 5.0e-3;

        Solver_Dict.time_steps_per_monitoring_interval =   2;

        Solver_Dict.number_of_frames =                     (1480 / 2);

        Solver_Dict.number_of_runs =                       1;

        Solver_Dict.max_fill_loops =                       5;

        Solver_Dict.data_drifting_on =                     1;

        Solver_Dict.number_of_drift_loops =                1;

        Solver_Dict.face_swap_diffusion_on =               0;

        Solver_Dict.recurrence_process_on =                1;

        Solver_Dict.control_conc_sum_on =                  1;

        Solver_Dict.balance_correction_on =                1;
    }

    void rCFD_user_set_Phase_Dict(void)
    {
#if RP_NODE
        int i_phase;

        loop_phases{

            if(i_phase == ph_gas){

                Phase_Dict[i_phase].number_of_data = 0;

                Phase_Dict[i_phase].time_step = Solver_Dict.global_time_step;

                Phase_Dict[i_phase].density = 1.0;

                Phase_Dict[i_phase].vof_max = 1.0;
            }

            if(i_phase == ph_solid){

                Phase_Dict[i_phase].number_of_data = 2;

                Phase_Dict[i_phase].time_step = Solver_Dict.global_time_step;

                Phase_Dict[i_phase].density = 1450.0;

                Phase_Dict[i_phase].vof_max = 0.65;

                Phase_Dict[i_phase].hindered_drift_on = 1;
            }
        }
#endif
    }

    void rCFD_user_set_Tracer_Dict(void)
    {
#if RP_NODE

        Tracer_Dict.format = guide_by_value_format;
#endif
    }

    void rCFD_user_set_Rec_Dict(void)
    {
        Rec_Dict.format = replay_format;
    }

    void rCFD_user_set_Data_Dict(void)
    {
#if RP_NODE
        int i_phase, i_data;

        loop_phases{

            loop_data{

                if((i_phase == ph_solid) && (i_data == c_drift_y_pos)){

                    Data_Dict[i_phase][i_data].type = concentration_data;

                    Data_Dict[i_phase][i_data].drifting_data = 1;

                    Data_Dict[i_phase][i_data].drift_velocity = (double*)malloc( 3 * sizeof(double));

                    Data_Dict[i_phase][i_data].drift_velocity[0] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[1] = 0.0032;

                    Data_Dict[i_phase][i_data].drift_velocity[2] = 0.0;

                }

                if((i_phase == ph_solid) && (i_data == c_drift_y_neg)){

                    Data_Dict[i_phase][i_data].type = concentration_data;

                    Data_Dict[i_phase][i_data].drifting_data = 1;

                    Data_Dict[i_phase][i_data].drift_velocity = (double*)malloc( 3 * sizeof(double));

                    Data_Dict[i_phase][i_data].drift_velocity[0] = 0.0;

                    Data_Dict[i_phase][i_data].drift_velocity[1] = -0.02;

                    Data_Dict[i_phase][i_data].drift_velocity[2] = 0.0;
                }

                if((i_phase == ph_solid) && (i_data == c_drift_0)){

                    Data_Dict[i_phase][i_data].type = concentration_data;

                    Data_Dict[i_phase][i_data].drifting_data = 0;
                }

            }
        }
#endif
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

                C_UDMI(i_cell, t, i_UDMI) = _C.average_velocity[i_phase][i_cell];

                C_UDMI(i_cell, t, (i_UDMI+1)) = _C.crossing_time[i_phase][i_cell];

                i_UDMI += 2;
            }

        }end_c_loop_int(i_cell, t)}}
#endif
    }

    void rCFD_user_init_Data(const short i_layer)
    {
#if RP_NODE

        int     i_phase, i_cell, i_data;

        loop_phases{

            loop_cells{

                loop_data{

                    if((i_phase == ph_solid) && (i_data == c_drift_y_pos)){

                        _C.data[_i_data] = 3./4.;
                    }

                    if((i_phase == ph_solid) && (i_data == c_drift_y_neg)){

                        _C.data[_i_data] = 1./4.;
                    }
                }
            }
        }

        Message0("\n\n...rCFD_user_init_Data\n");
#endif
    }

#endif

    /* user accesses Data */
    /*************************************************************************************/

#if 1

    void rCFD_user_access_data_before_shift(const short i_phase, const short i_layer)
    {
#if RP_NODE

        if(Solver_Dict.verbose){

            Message0("\n\nUser before Shift: i_phase %d, i_layer %d", i_phase, i_layer);
        }
#endif
    }

    void rCFD_user_post()
    {
#if RP_NODE

#if 1   /* local vars */
        int     i_layer, i_phase, i_cell, i_data, i_face, i_frame, i_diff_loop, i_UDMI;

        int     c0, c1, c, i_frame_c0, i_frame_c1;

        double  mass_of_large_part = 0.0, y_mean_of_large_part = 0.0, mass_of_large_part_ref = 0.0, y_mean_of_large_part_ref = 0.0;

        double  mass_of_small_part = 0.0, y_mean_of_small_part = 0.0, mass_of_small_part_ref = 0.0, y_mean_of_small_part_ref = 0.0;

        double  mass_of_part = 0.0, y_mean_of_part = 0.0, mass_of_part_ref = 0.0, y_mean_of_part_ref = 0.0;

        double  number_of_large_part, number_of_small_part;

        double  volume_large, volume_small, vol_flip;

        double  *number_count_large = NULL, *number_count_small = NULL, *conc_large = NULL, *conc_small = NULL, *swap_large = NULL, *swap_small = NULL;

        Domain  *d=Get_Domain(1);

        Thread  *t = NULL, *t_mix = NULL, *t_gas = NULL;
#endif

        i_layer = 0;

        if(Solver.current_layer > 0){

            rCFD_map_from_to_layer(Solver.current_layer, i_layer);
        }

        /* 1. allocate number_count_large/small, conc_large/small and corresponding swap vars */
        {
            number_count_small = (double*)malloc(_Cell_Dict.number_of_int_cells * sizeof(double));
            number_count_large = (double*)malloc(_Cell_Dict.number_of_int_cells * sizeof(double));
            conc_small = (double*)malloc(_Cell_Dict.number_of_int_cells * sizeof(double));
            conc_large = (double*)malloc(_Cell_Dict.number_of_int_cells * sizeof(double));
            swap_small = (double*)malloc(_Cell_Dict.number_of_int_cells * sizeof(double));
            swap_large = (double*)malloc(_Cell_Dict.number_of_int_cells * sizeof(double));

            thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){

                t_mix = t;

                t_gas = THREAD_SUB_THREAD(t, 0);

                begin_c_loop_int(i_cell, t){

                    number_count_large[i_cell] = C_U(i_cell, t_gas);

                    number_count_small[i_cell] = C_V(i_cell, t_gas);

                    i_frame = Rec.global_frame[_C.island_id[i_cell]];

                    if((_C.volume[i_cell] * _C.vof[i_frame][i_cell][ph_solid]) > 0.0){

                        conc_small[i_cell] = (number_count_small[i_cell] * diam_small * diam_small * diam_small * M_PI/6.0) /

                            (_C.volume[i_cell] * _C.vof[i_frame][i_cell][ph_solid]);

                        conc_large[i_cell] = (number_count_large[i_cell] * diam_large * diam_large * diam_large * M_PI/6.0) /

                            (_C.volume[i_cell] * _C.vof[i_frame][i_cell][ph_solid]);
                    }
                    else{

                        conc_small[i_cell] = 0.0;

                        conc_large[i_cell] = 0.0;
                    }

                    swap_large[i_cell] = 0.0;

                    swap_small[i_cell] = 0.0;

                }end_c_loop_int(i_cell, t)
            }}
        }

        /* 3. eval y_mean_large/small (rCFD data and ref data) */
        {
            /* rCFD data */
            {
                i_phase = ph_solid;

                loop_int_cells{

                    i_frame = Rec.global_frame[_C.island_id[i_cell]];

                    i_data = c_drift_y_neg;

                    mass_of_large_part += _C.data[_i_data] * _C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density;

                    y_mean_of_large_part += _C.x[i_cell][1] * _C.data[_i_data] * _C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density;

                    i_data = c_drift_y_pos;

                    mass_of_small_part += _C.data[_i_data] * _C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density;

                    y_mean_of_small_part += _C.x[i_cell][1] * _C.data[_i_data] * _C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density;

                    mass_of_part += _C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density;

                    y_mean_of_part += _C.x[i_cell][1] * _C.volume[i_cell] * _C.vof[_i_vof] * Phase_Dict[i_phase].density;
                }

                mass_of_large_part =    PRF_GRSUM1(mass_of_large_part);

                y_mean_of_large_part =  PRF_GRSUM1(y_mean_of_large_part);

                mass_of_small_part =    PRF_GRSUM1(mass_of_small_part);

                y_mean_of_small_part =  PRF_GRSUM1(y_mean_of_small_part);

                mass_of_part =          PRF_GRSUM1(mass_of_part);

                y_mean_of_part =        PRF_GRSUM1(y_mean_of_part);

                if(mass_of_large_part > 0.0){

                    y_mean_of_large_part /= mass_of_large_part;
                }

                if(mass_of_small_part > 0.0){

                    y_mean_of_small_part /= mass_of_small_part;
                }

                if(mass_of_part > 0.0){

                    y_mean_of_part /= mass_of_part;
                }
            }

            /* ref data */
            {
                loop_int_cells{

                    mass_of_large_part_ref += number_count_large[i_cell] * diam_large * diam_large * diam_large * M_PI / 6.0 * Phase_Dict[i_phase].density;

                    y_mean_of_large_part_ref += _C.x[i_cell][1] * number_count_large[i_cell] * diam_large * diam_large * diam_large * M_PI / 6.0 * Phase_Dict[i_phase].density;


                    mass_of_small_part_ref += number_count_small[i_cell] * diam_small * diam_small * diam_small * M_PI / 6.0 * Phase_Dict[i_phase].density;

                    y_mean_of_small_part_ref += _C.x[i_cell][1] * number_count_small[i_cell] * diam_small * diam_small * diam_small * M_PI / 6.0 * Phase_Dict[i_phase].density;


                    mass_of_part_ref += M_PI / 6.0 * Phase_Dict[i_phase].density *

                        (number_count_large[i_cell] * diam_large * diam_large * diam_large + number_count_small[i_cell] * diam_small * diam_small * diam_small);

                    y_mean_of_part_ref += _C.x[i_cell][1] *  M_PI / 6.0 * Phase_Dict[i_phase].density *

                        (number_count_large[i_cell] * diam_large * diam_large * diam_large + number_count_small[i_cell] * diam_small * diam_small * diam_small);
                }

                mass_of_large_part_ref =    PRF_GRSUM1(mass_of_large_part_ref);

                y_mean_of_large_part_ref =  PRF_GRSUM1(y_mean_of_large_part_ref);

                mass_of_small_part_ref =    PRF_GRSUM1(mass_of_small_part_ref);

                y_mean_of_small_part_ref =  PRF_GRSUM1(y_mean_of_small_part_ref);

                mass_of_part_ref =          PRF_GRSUM1(mass_of_part_ref);

                y_mean_of_part_ref =        PRF_GRSUM1(y_mean_of_part_ref);

                if(mass_of_large_part_ref > 0.0){

                    y_mean_of_large_part_ref /= mass_of_large_part_ref;
                }

                if(mass_of_small_part_ref > 0.0){

                    y_mean_of_small_part_ref /= mass_of_small_part_ref;
                }

                if(mass_of_part_ref > 0.0){

                    y_mean_of_part_ref /= mass_of_part_ref;
                }
            }

        }

        /* 4. node-0 writes values to monitor file */
        if(myid == 0){

            FILE    *f_out = NULL;

            if(rCFD_monitors_initiated == 0){

                f_out = fopen("./post/monitor_rCFD.out","w");

                rCFD_monitors_initiated = 1;
            }
            else{
                f_out = fopen("./post/monitor_rCFD.out","a");
            }

            if(f_out == NULL){

                Message0("\nERROR: Could not open monitor file\n");

                return;
            }

            fprintf(f_out, "%e %e %e %e %e %e %e %e %e %e %e\n", Solver.global_time,

                y_mean_of_part, y_mean_of_large_part, y_mean_of_small_part,

                y_mean_of_part_ref, y_mean_of_large_part_ref, y_mean_of_small_part_ref,

                mass_of_part_ref, mass_of_large_part_ref, mass_of_small_part_ref, mass_of_part);

            fclose(f_out);
        }

        /* 5. diffuse reference particle information to ref. vof field */
        {
            i_phase = ph_solid;

            for(i_diff_loop = 0; i_diff_loop < 5; i_diff_loop++){

                /* initialize swap's */
                loop_int_cells{

                    swap_large[i_cell] = 0;

                    swap_small[i_cell] = 0;
                }

                /* define swap masses */
                loop_int_faces{

                    c0 = _F.c0[i_face];
                    c1 = _F.c1[i_face];

                    if((c0 < _Cell_Dict.number_of_int_cells) && (c1 < _Cell_Dict.number_of_int_cells)){

                        /* large reference particles */
                        {

                            if(conc_large[c0] > conc_large[c1]){

                                c = c0; c0 = c1; c1 = c;
                            }

                            if(conc_large[c1] > conc_large[c0]){

                                i_frame_c0 = Rec.global_frame[_C.island_id[c0]];
                                i_frame_c1 = Rec.global_frame[_C.island_id[c1]];

                                if((_C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase]) < (_C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase])){

                                    vol_flip = _C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase];
                                }
                                else{
                                    vol_flip = _C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase];
                                }

                                swap_large[c1] -= vol_flip * conc_large[c1] / 2. * (1./8.);

                                swap_large[c0] += vol_flip * conc_large[c1] / 2. * (1./8.);
                            }
                        }

                        /* small reference particles */
                        {

                            if(conc_small[c0] > conc_small[c1]){

                                c = c0; c0 = c1; c1 = c;
                            }

                            if(conc_small[c1] > conc_small[c0]){

                                i_frame_c0 = Rec.global_frame[_C.island_id[c0]];
                                i_frame_c1 = Rec.global_frame[_C.island_id[c1]];

                                if((_C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase]) < (_C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase])){

                                    vol_flip = _C.volume[c0] * _C.vof[i_frame_c0][c0][i_phase];
                                }
                                else{
                                    vol_flip = _C.volume[c1] * _C.vof[i_frame_c1][c1][i_phase];
                                }

                                swap_small[c1] -= vol_flip * conc_small[c1] / 2. * (1./8.);

                                swap_small[c0] += vol_flip * conc_small[c1] / 2. * (1./8.);
                            }
                        }
                    }
                }

                /* update data by swap's */
                loop_int_cells{

                    i_frame = Rec.global_frame[_C.island_id[i_cell]];

                    if((_C.volume[i_cell] * _C.vof[_i_vof]) > 0.0){

                        conc_large[i_cell] += swap_large[i_cell] / (_C.volume[i_cell] * _C.vof[_i_vof]);

                        conc_small[i_cell] += swap_small[i_cell] / (_C.volume[i_cell] * _C.vof[_i_vof]);

                    }
                }
            }
        }

        /* 6. nodes write to UDMI */
        {
            volume_large = M_PI / 6.0 * diam_large * diam_large * diam_large;

            volume_small = M_PI / 6.0 * diam_small * diam_small * diam_small;

            loop_int_cells{

                i_frame = Rec.global_frame[_C.island_id[i_cell]];

                i_phase = ph_solid;

                i_data = c_drift_y_pos;  /* small particles */

                C_UDMI(i_cell, t_mix, c_vof_small) = _C.data[_i_data] * _C.vof[_i_vof];

                i_data = c_drift_y_neg;  /* small particles */

                C_UDMI(i_cell, t_mix, c_vof_large) = _C.data[_i_data] * _C.vof[_i_vof];

                C_UDMI(i_cell, t_mix, c_vof) = _C.vof[_i_vof];

                /* Sauter mean diameter d32 */
                if(_C.vof[_i_vof] > 0.1){

                    i_data = c_drift_y_neg;

                    number_of_large_part = _C.data[_i_data] * _C.volume[i_cell] * _C.vof[_i_vof] / volume_large;

                    i_data = c_drift_y_pos;

                    number_of_small_part = _C.data[_i_data] * _C.volume[i_cell] * _C.vof[_i_vof] / volume_small;

                    if((number_of_large_part + number_of_small_part) >= 1.0e-10){

                        C_UDMI(i_cell, t_mix, c_d23) =

                            (number_of_large_part * diam_large * diam_large * diam_large +

                             number_of_small_part * diam_small * diam_small * diam_small) /

                            (number_of_large_part * diam_large * diam_large + number_of_small_part * diam_small * diam_small);
                    }
                    else{
                        C_UDMI(i_cell, t_mix, c_d23) = 0.0;
                    }
                }
                else{

                    C_UDMI(i_cell, t_mix, c_d23) = 0.0;
                }

                if((conc_large[i_cell] + conc_small[i_cell]) > 0.0){

                    C_UDMI(i_cell, t_mix, c_vof_small_ref) = C_UDMI(i_cell, t_mix, c_vof) *

                        (conc_small[i_cell] / (conc_small[i_cell] + conc_large[i_cell]));

                    C_UDMI(i_cell, t_mix, c_vof_large_ref) = C_UDMI(i_cell, t_mix, c_vof) *

                        (conc_large[i_cell] / (conc_small[i_cell] + conc_large[i_cell]));

                    number_of_large_part = conc_large[i_cell] * _C.volume[i_cell] * _C.vof[_i_vof] / volume_large;

                    number_of_small_part = conc_small[i_cell] * _C.volume[i_cell] * _C.vof[_i_vof] / volume_small;

                    if(_C.vof[_i_vof] > 0.1){

                        C_UDMI(i_cell, t_mix, c_d23_ref) =

                            (number_of_large_part * diam_large * diam_large * diam_large +

                             number_of_small_part * diam_small * diam_small * diam_small) /

                            (number_of_large_part * diam_large * diam_large + number_of_small_part * diam_small * diam_small);
                    }
                    else{

                        C_UDMI(i_cell, t_mix, c_d23_ref) = 0.0;
                    }
                }
                else{

                    C_UDMI(i_cell, t_mix, c_vof_small_ref) = 0.0;

                    C_UDMI(i_cell, t_mix, c_vof_large_ref) = 0.0;

                    C_UDMI(i_cell, t_mix, c_d23_ref) = 0.0;
                }
            }
        }

        /* 2. average vof and conc in between 4s and 7s */
        {
            if((Solver.global_time > 4.0) && (Solver.global_time < 7.0)){

                i_phase = ph_solid;

                loop_int_cells{

                    i_frame = Rec.global_frame[_C.island_id[i_cell]];

                    C_UDMI(i_cell, t_mix, c_vof_average) = (double)average_counter * C_UDMI(i_cell, t_mix, c_vof_average) + _C.vof[_i_vof];

                    C_UDMI(i_cell, t_mix, c_vof_small_average) = (double)average_counter * C_UDMI(i_cell, t_mix, c_vof_small_average) + C_UDMI(i_cell, t_mix, c_vof_small);

                    C_UDMI(i_cell, t_mix, c_vof_large_average) = (double)average_counter * C_UDMI(i_cell, t_mix, c_vof_large_average) + C_UDMI(i_cell, t_mix, c_vof_large);

                    C_UDMI(i_cell, t_mix, c_vof_small_ref_average) = (double)average_counter * C_UDMI(i_cell, t_mix, c_vof_small_ref_average) + C_UDMI(i_cell, t_mix, c_vof_small_ref);

                    C_UDMI(i_cell, t_mix, c_vof_large_ref_average) = (double)average_counter * C_UDMI(i_cell, t_mix, c_vof_large_ref_average) + C_UDMI(i_cell, t_mix, c_vof_large_ref);

                    for(i_UDMI = c_vof_average; i_UDMI <= c_vof_large_ref_average; i_UDMI++){

                        C_UDMI(i_cell, t_mix, i_UDMI) /= (double)(average_counter + 1);
                    }
                }

                average_counter++;
            }
        }

        /* 7. free local vars */
        {
            free(number_count_large);
            free(number_count_small);

            free(swap_large);
            free(swap_small);

            free(conc_large);
            free(conc_small);
        }

        Message0("\n\n...rCFD_monitors\n");
#endif
    }


#endif

    /* register user functions */
    /*************************************************************************************/

    void rCFD_user_register_functions()
    {
        REGISTER_RCFD_UDF( rCFD_user_set_Solver_Dict );
        REGISTER_RCFD_UDF( rCFD_user_set_Phase_Dict );
        REGISTER_RCFD_UDF( rCFD_user_set_Tracer_Dict );
        REGISTER_RCFD_UDF( rCFD_user_set_Rec_Dict );
        REGISTER_RCFD_UDF( rCFD_user_set_Data_Dict );
        REGISTER_RCFD_UDF( rCFD_user_pre_proc );
        REGISTER_RCFD_UDF( rCFD_user_init_Data );
        REGISTER_RCFD_UDF( rCFD_user_access_data_before_shift );
        REGISTER_RCFD_UDF( rCFD_user_post );
    }

#endif
