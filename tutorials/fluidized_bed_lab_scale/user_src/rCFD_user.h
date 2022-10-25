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

        static double   _mean_solid_marker_conc = 0.0;
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

#define _reset_data_fields_interval ((int)(15./Solver_Dict.global_time_step))

#define _turbulence_intensity   0.1

#define _adapt_outgoing_conc_of_gas_A   0.6

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

        Solver_Dict.number_of_frames =                     100;

        Solver_Dict.number_of_runs =                       1;

        Solver_Dict.data_drifting_on =                     1;

        Solver_Dict.face_swap_diffusion_on =               1;

        Solver_Dict.recurrence_process_on =                1;

        Solver_Dict.face_swap_max_per_loop =               (1./8.);

        Solver_Dict.control_conc_sum_on =                  1;
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

                Phase_Dict[i_phase].number_of_data = 1;

                Phase_Dict[i_phase].time_step = Solver_Dict.global_time_step;

                Phase_Dict[i_phase].density = 1000.0;
            }
        }
#endif
    }

    void rCFD_user_set_Tracer_Dict(void)
    {
#if RP_NODE
        Tracer_Dict.format = guide_by_value_format;

        //Tracer_Dict.random_walk[gas] = 1;

        //Tracer_Dict.random_walk_velocity_ratio[gas] = 0.1;

        Tracer_Dict.avoid_information_lock_cells_on = 1;
#endif
    }

    void rCFD_user_set_Rec_Dict(void)
    {
        Rec_Dict.monitor_rec_frames_on = 1;

        Rec_Dict.adapt_vof_stitching_on = 1;

        Rec_Dict.number_of_adapt_vof_loops = 10;

        Rec_Dict.format = off_diagonal_band_format;

        Rec_Dict.off_diagonal_band_width = 30;
    }

    void rCFD_user_set_Data_Dict(void)
    {
#if RP_NODE
        int i_phase, i_data;

        loop_phases{

            loop_data{

                if((i_phase == gas) && (i_data == c_gas_A)){

                    Data_Dict[i_phase][i_data].type = concentration_data;

                    Data_Dict[i_phase][i_data].face_swap_diffusion = -1./8.;
                }

                if((i_phase == gas) && (i_data == c_gas_B)){

                    Data_Dict[i_phase][i_data].type = concentration_data;

                    Data_Dict[i_phase][i_data].face_swap_diffusion = -1./8.;
                }

                if((i_phase == solid) && (i_data == c_drift_0)){

                    Data_Dict[i_phase][i_data].type = concentration_data;

                    //Data_Dict[i_phase][i_data].face_swap_diffusion = 1.0;
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

#if 1
        int     i_phase, i_cell, i_data, i_dim, i_frame;

        double  x[3], conc_volume, volume;
#endif

        /* 1. set concentration values */
        loop_phases{

            loop_cells{

                loop_dim{

                    x[i_dim] = _C.x[i_cell][i_dim];
                }

                loop_data{

                    if(i_phase == solid){

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

        /* 2. calc. initial mean solid concentration (needed for mixing index) */
        {
            volume = 0.0;

            conc_volume = 0.0;

            i_phase = solid;

            i_data = c_drift_0;

            loop_int_cells{

                i_frame = Rec.global_frame[_C.island_id[i_cell]];

                conc_volume += _C.data[_i_data] * _C.vof[_i_vof] * _C.volume[i_cell];

                volume += _C.vof[_i_vof] * _C.volume[i_cell];
            }

            conc_volume = PRF_GRSUM1(conc_volume);

            volume = PRF_GRSUM1(volume);

            if(volume > 0.0){

                _mean_solid_marker_conc = conc_volume / volume;
            }
            else{

                _mean_solid_marker_conc = 0.0;
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

        if(i_phase == gas){

#if 1       // local vars

            int     i_cell, i_frame;

            double  volume_in, volume_in_global, volume_in2, volume_in2_global, volume_out, volume_out_global;

            double  mean_conc_out_A, mean_conc_out_B, sum_of_outgoing_conc;

            double  primary_inflow, secondary_inflow;

#endif

            // initializations
            {
                volume_in = 0.0; volume_in2 = 0.0; volume_out = 0.0;

                mean_conc_out_A = 0.0, mean_conc_out_B = 0.0;

                primary_inflow = 3.28e-3; secondary_inflow = 7.14e-5;      // (kg/s)
            }

            // set conc at inlets, calc volumes, get mean conc at outlet
            {
                loop_cells{

                    i_frame = Rec.global_frame[_C.island_id[i_cell]];

                    if(_C.x[i_cell][2] < 0.02){

                        _C.data[i_phase][i_cell][c_gas_A] = 0.0;

                        _C.data[i_phase][i_cell][c_gas_B] = 1.0;

                        volume_in += _C.volume[i_cell] * _C.vof[_i_vof];
                    }

                    if(_C.x[i_cell][1] < -0.075){

                        _C.data[i_phase][i_cell][c_gas_A] = 1.0;

                        _C.data[i_phase][i_cell][c_gas_B] = 0.0;

                        volume_in2 += _C.volume[i_cell] * _C.vof[_i_vof];
                    }

                    if(_C.x[i_cell][2] > 0.4){

                        mean_conc_out_A += _C.data[gas][i_cell][c_gas_A] * (_C.volume[i_cell] * _C.vof[_i_vof]);

                        mean_conc_out_B += _C.data[gas][i_cell][c_gas_B] * (_C.volume[i_cell] * _C.vof[_i_vof]);

                        volume_out += _C.volume[i_cell] * _C.vof[_i_vof];
                    }
                }

                volume_in_global =  PRF_GRSUM1(volume_in);

                volume_in2_global = PRF_GRSUM1(volume_in2);

                volume_out_global = PRF_GRSUM1(volume_out);
            }

            // inlet mass sources
            {
                Balance[gas][c_gas_A].mass_source += secondary_inflow * Solver.timestep_width[i_layer] * volume_in2 / volume_in2_global; // (kg, per node )

                Balance[gas][c_gas_B].mass_source += primary_inflow * Solver.timestep_width[i_layer] * volume_in / volume_in_global;
            }

            // outlet mass sources
            {
                if((volume_out > 0) && (volume_out_global > 0.0)){

                    mean_conc_out_A /= volume_out;

                    mean_conc_out_B /= volume_out;

                    mean_conc_out_A *= _adapt_outgoing_conc_of_gas_A;

                    sum_of_outgoing_conc = (mean_conc_out_A + mean_conc_out_B);

                    if(sum_of_outgoing_conc > 0.0){

                        mean_conc_out_A /= sum_of_outgoing_conc;

                        mean_conc_out_B /= sum_of_outgoing_conc;
                    }

                    Balance[gas][c_gas_A].mass_source -= mean_conc_out_A * (primary_inflow + secondary_inflow) * Solver.timestep_width[i_layer] *

                         (volume_out / volume_out_global);  // (kg, per node)

                    Balance[gas][c_gas_B].mass_source -= mean_conc_out_B * (primary_inflow + secondary_inflow) * Solver.timestep_width[i_layer] *

                         (volume_out / volume_out_global);
                }
            }
        }

        if(Solver_Dict.verbose){

            Message0("\n\nUser before Shift: i_phase %d, i_layer %d", i_phase, i_layer);
        }
#endif
    }

    void rCFD_user_post()
    {
#if RP_NODE

#if 1   /* local vars */
        int i_layer, i_phase, i_cell, i_dim, i_data, i_UDMI, i_frame;

        double x[3], mass_c_gas_A, mixing_index, mixing_nom, mixing_denom, volume, conc_volume;

        Domain  *d=Get_Domain(1);

        Thread  *t;
#endif

        i_layer = 0;

        if(Solver.current_layer > 0){

            rCFD_map_from_to_layer(Solver.current_layer, i_layer);
        }

        /* 1. event. reset fields, mean conc. & balances */
        {
            if((Solver.global_run_counter % _reset_data_fields_interval) == 0){

                /* data fields */
                loop_cells{

                    loop_dim{

                        x[i_dim] = _C.x[i_cell][i_dim];
                    }

                    i_phase = solid;

                    i_data = c_drift_0;

                    if(x[1] > 0.0){

                        _C.data[_i_data] = 1.0;
                    }
                    else{

                        _C.data[_i_data] = 0.0;
                    }

                    i_phase = gas;

                    if(x[1] < -0.074){

                        _C.data[i_phase][i_cell][c_gas_A] = 1.0;
                        _C.data[i_phase][i_cell][c_gas_B] = 0.0;
                    }
                    else{

                        _C.data[i_phase][i_cell][c_gas_A] = 0.0;
                        _C.data[i_phase][i_cell][c_gas_B] = 1.0;
                    }
                }

                /* mean solid concentration */
                {
                    volume = 0.0;

                    conc_volume = 0.0;

                    i_phase = solid;

                    i_data = c_drift_0;

                    loop_int_cells{

                        i_frame = Rec.global_frame[_C.island_id[i_cell]];

                        conc_volume += _C.data[_i_data] * _C.vof[_i_vof] * _C.volume[i_cell];

                        volume += _C.vof[_i_vof] * _C.volume[i_cell];
                    }

                    conc_volume = PRF_GRSUM1(conc_volume);

                    volume = PRF_GRSUM1(volume);

                    if(volume > 0.0){

                        _mean_solid_marker_conc = conc_volume / volume;
                    }
                    else{

                        _mean_solid_marker_conc = 0.0;
                    }
                }

                /* balances */
                {
                    loop_phases{

                        loop_data{

                            Balance[_i_balance].mass_integral = 0.0;
                        }
                    }

                    loop_int_cells{

                        i_frame = Rec.global_frame[_C.island_id[i_cell]];

                        loop_phases{

                            loop_data{

                                Balance[_i_balance].mass_integral += _C.data[_i_data] *

                                    Phase_Dict[i_phase].density * _C.volume[i_cell] * _C.vof[_i_vof];
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
                }
            }
        }

        /* 2. eval mass_c_gas_A and mixing_index */
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

            loop_int_cells{

                i_frame = Rec.global_frame[_C.island_id[i_cell]];

                mixing_nom += fabs( _C.data[i_phase][i_cell][c_drift_0] - _mean_solid_marker_conc) * _C.volume[i_cell] * _C.vof[_i_vof];

                mixing_denom += _C.volume[i_cell]  * _C.vof[_i_vof];
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

        /* 3. node-0 writes values to ref file */
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

        /* 4. nodes write species conc to UDMI */
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

    /* register user functions */
    /*************************************************************************************/

    void rCFD_user_register_functions()
    {
        REGISTER_RCFD_UDF( rCFD_user_set_Solver_Dict );
        REGISTER_RCFD_UDF( rCFD_user_set_Phase_Dict );
        REGISTER_RCFD_UDF( rCFD_user_set_Tracer_Dict );
        REGISTER_RCFD_UDF( rCFD_user_set_Rec_Dict );
        REGISTER_RCFD_UDF( rCFD_user_set_Data_Dict );
        REGISTER_RCFD_UDF( rCFD_user_set_Balance_Dict );
        REGISTER_RCFD_UDF( rCFD_user_pre_proc );
        REGISTER_RCFD_UDF( rCFD_user_init_Data );
        REGISTER_RCFD_UDF( rCFD_user_access_data_before_shift );
        REGISTER_RCFD_UDF( rCFD_user_post );
    }

#endif
