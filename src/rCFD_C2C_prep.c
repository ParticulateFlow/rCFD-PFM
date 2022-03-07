#include <stdio.h>

#include "rCFD_types.h"
#include "rCFD_globals.h"
#include "rCFD_parallel.h"
#include "rCFD_defaults.h"
#include "rCFD_macros.h"
#include "rCFD_init.h"
#include "rCFD_free.h"

#include "rCFD_user.h"

/* (C)  2021
    Stefan Pirker
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/


/* Contents:

    rCFD_init_all

    rCFD_analyse_CFD

    rCFD_write_Tracer_Positions



    rcfd_init_tracers

    rCFD_update_Tracers

    rcfd_no_standard_drag

    rCFD_guide_Tracers

    rCFD_write_C2Cs


    rCFD_write_Norms

    rCFD_write_Rec

    rCFD_write_Topo (still missing)

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
DEFINE_EXECUTE_AT_END(rCFD_analyse_CFD)
/*************************************************************************************/
{
#if RP_NODE
    Domain      *d=Get_Domain(1);
    Thread      *t, *t_phase;
    cell_t      c;

    int         i_layer = 0, i_phase, i_cell;

    int         time_steps_per_monitoring_interval;

    double      v_mag, w0, L;
    double      *dt_cross_min, *dt_cross_max, dt_cross_min_global;


    dt_cross_min = (double*)malloc(Solver_Dict.number_of_phases * sizeof(double));
    dt_cross_max = (double*)malloc(Solver_Dict.number_of_phases * sizeof(double));

    loop_phases{

        dt_cross_min[i_phase] = 1.e10;
        dt_cross_max[i_phase] = -1.e10;
    }

    thread_loop_c(t,d){begin_c_loop_int(c,t){

        i_cell = (int)c;

        loop_phases{

            if(Solver_Dict.number_of_phases == 1){

                t_phase = t;
            }
            else if(THREAD_SUB_THREAD(t,i_phase) != NULL) {

                t_phase = THREAD_SUB_THREAD(t,i_phase);
            }
            else{

                Message("\nERROR myid %d rCFD_analyse_CFD THREAD_SUB_THREAD(t,i_phase) == NULL", myid);

                return;
            }

            v_mag = sqrt(C_U(c,t_phase)*C_U(c,t_phase) + C_V(c,t_phase)*C_V(c,t_phase) + C_W(c,t_phase)*C_W(c,t_phase));

            w0 = (double)Solver_Dict.analyse_CFD_count;

            _C.average_velocity[i_phase][i_cell] = (w0 * _C.average_velocity[i_phase][i_cell] + v_mag) / (w0 + 1.0);

            L = pow(C_VOLUME(c,t),1./3.);

            if(_C.average_velocity[i_phase][i_cell] > 0.0){

                _C.crossing_time[i_phase][i_cell] = L/_C.average_velocity[i_phase][i_cell];
            }
            else{

                _C.crossing_time[i_phase][i_cell] = 1.e10; /* something large */
            }

            if((dt_cross_max[i_phase] < _C.crossing_time[i_phase][i_cell]) && (_C.crossing_time[i_phase][i_cell] < 1.e10)){

                dt_cross_max[i_phase] = _C.crossing_time[i_phase][i_cell];
            }

            if(dt_cross_min[i_phase] > _C.crossing_time[i_phase][i_cell]){

                dt_cross_min[i_phase] = _C.crossing_time[i_phase][i_cell];
            }
        }

    }end_c_loop_int(c,t)}

    /* propose recurrence time-step */
    {
        dt_cross_min_global = 1.e10;

        loop_phases{

            dt_cross_max[i_phase] = PRF_GRHIGH1(dt_cross_max[i_phase]);
            dt_cross_min[i_phase] = PRF_GRLOW1(dt_cross_min[i_phase]);

            if(dt_cross_min[i_phase] < dt_cross_min_global){

                dt_cross_min_global = dt_cross_min[i_phase];
            }

            Message0("\n...rCFD_analyse_CFD: cell_crossing_time for phase %d (based on %d time-steps): [%e, %e] ... \n",

                i_phase, (Solver_Dict.analyse_CFD_count + 1), dt_cross_min[i_phase], dt_cross_max[i_phase]);
        }

        time_steps_per_monitoring_interval = (int)(dt_cross_min_global * (double)Solver_Dict.max_number_of_cells_per_time_step / CURRENT_TIMESTEP);

        if(time_steps_per_monitoring_interval < 1){

            time_steps_per_monitoring_interval = 1;

            Message0("\n... rCFD_analyse_CFD: WARNING: set monitoring interval = 1 ... \n");
        }

        Message0("\n...rCFD_analyse_CFD: Recommend %d time_steps_per_monitoring_interval (based on %d time-steps) ... \n",

            time_steps_per_monitoring_interval, (Solver_Dict.analyse_CFD_count + 1));
    }

    free(dt_cross_max);
    free(dt_cross_min);

    Solver_Dict.analyse_CFD_count++;

    rCFD_user_pre_proc();
#endif
}

/*************************************************************************************/
DEFINE_ON_DEMAND(rCFD_write_Tracer_Positions)
/*************************************************************************************/
{
#if RP_NODE
    Domain  *d=Get_Domain(1);
    Thread  *t;
    cell_t  c;

    int     number_of_start_positions, number_of_start_position_coords, total_number_of_start_postions;

    double  *start_position_coords = NULL, x[3];
#endif

    /* A: set Tracer_Dict, determine number of start_positions & allocate storage */
    {
#if RP_NODE

        number_of_start_positions = 0;

        thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){begin_c_loop_int(c,t){

            C_CENTROID(x,c,t);

            if  (no_ROI_defined || cell_in_ROI){

                number_of_start_positions++;
            }

        }end_c_loop_int(c,t)}}

        number_of_start_positions *= Solver_Dict.number_of_phases * Tracer_Dict.number_of_Tracers_per_cell;

        total_number_of_start_postions = PRF_GISUM1(number_of_start_positions);

        number_of_start_position_coords = 3. * number_of_start_positions;

        start_position_coords = (double*)malloc( number_of_start_position_coords * sizeof(double));

            Message("\nmyid %d number_of_start_positions %d", myid, number_of_start_positions);

#endif
    }

    /* B: determine position in cell, fill local start_position_coords */
    {
#if RP_NODE
        int i_arr, i_node, i_tracer, i_dim;
        double rand_real;
        int nodes_in_cell, n1_index, n2_index;
        double  x[3],n1[3],n2[3],x_n1[3],x_n2[3],p1[3],p2[3],p1_p2[3],p_start[3];
        Node *n;

        i_arr = 0;

        thread_loop_c(t,d){ if(FLUID_CELL_THREAD_P(t)){ begin_c_loop_int(c,t){

            C_CENTROID(x,c,t);

            if(no_ROI_defined || cell_in_ROI){

                for(i_tracer = 0; i_tracer < Solver_Dict.number_of_phases * Tracer_Dict.number_of_Tracers_per_cell; i_tracer++){

                    /* P1: Chose two nodes */
                    {
                        nodes_in_cell = 0;
                        c_node_loop(c,t,i_node) nodes_in_cell++;

                        rand_real = (double)rand()/(double)RAND_MAX;            /* [0..1] */
                        n1_index = (int)((double)(nodes_in_cell+1)*rand_real);  /* [0..nodes_in_cell] */
                        n2_index = n1_index;

                        while(n2_index == n1_index){

                            rand_real = (double)rand()/(double)RAND_MAX;
                            n2_index = (int)((double)(nodes_in_cell+1)*rand_real);
                        }

                        for(i_dim = 0; i_dim < 3; i_dim++){

                            n1[i_dim]=0.; n2[i_dim]=0.;
                        }

                        c_node_loop(c,t,i_node){
                            if(i_node==n1_index) {n=C_NODE(c,t,i_node);n1[0]=NODE_X(n);n1[1]=NODE_Y(n);n1[2]=NODE_Z(n);}
                            if(i_node==n2_index) {n=C_NODE(c,t,i_node);n2[0]=NODE_X(n);n2[1]=NODE_Y(n);n2[2]=NODE_Z(n);}
                        }
                    }

                    /* P2: Chose two Points between Cell Center and Nodes (closer to nodes) */
                    {
                        NV_VV(x_n1,=,n1,-,x);
                        rand_real=(double)rand()/(double)RAND_MAX;
                        NV_V_VS(p1,=,x,+,x_n1,*,(1.-rand_real*rand_real));
                        NV_VV(x_n2,=,n2,-,x);
                        rand_real=(double)rand()/(double)RAND_MAX;
                        NV_V_VS(p2,=,x,+,x_n2,*,(1.-rand_real*rand_real));
                    }

                    /* P3: Chose final Point between P1 and P2 */
                    {
                        NV_VV(p1_p2,=,p2,-,p1);
                        rand_real=(double)rand()/(double)RAND_MAX;
                        NV_V_VS(p_start,=,p1,+,p1_p2,*,rand_real);

                        start_position_coords[i_arr] = p_start[0]; i_arr++;
                        start_position_coords[i_arr] = p_start[1]; i_arr++;
                        start_position_coords[i_arr] = p_start[2]; i_arr++;
                    }
                }
            }

        }end_c_loop_int(c,t)}}
#endif
    }

    /* C: communicate to Node-0, which writes to file */
    {
#if RP_NODE

        if (myid > 0){

            PRF_CSEND_INT(node_zero, &number_of_start_position_coords, 1, myid);
            PRF_CSEND_REAL(node_zero, start_position_coords, number_of_start_position_coords, myid);

            free(start_position_coords);
        }

        if (myid == 0){

            int i_node, number_of_lines, i_line;
            FILE    *f_out = NULL;

            f_out = fopen(File_Dict.tracer_start_position_filename,"w");

            if(f_out == NULL){

                Message0("\n... ERROR: rCFD_write_tracer_start_pos: fo == NULL ...\n");

                return;
            }

            number_of_lines = number_of_start_position_coords/3;

            for(i_line = 0; i_line < number_of_lines; i_line++){

                fprintf(f_out,"(( ");
                fprintf(f_out,"%f %f %f", start_position_coords[i_line*3], start_position_coords[i_line*3+1], start_position_coords[i_line*3+2]);
                fprintf(f_out," 0. 0. 0. 1.e-6 273. 1.e-10 ))\n");
            }

            free(start_position_coords);

            compute_node_loop_not_zero(i_node){

                PRF_CRECV_INT(i_node, &number_of_start_position_coords, 1, i_node);

                start_position_coords=(double*)malloc(number_of_start_position_coords * sizeof(double));

                PRF_CRECV_REAL(i_node, start_position_coords, number_of_start_position_coords, i_node);

                number_of_lines = number_of_start_position_coords/3;

                for(i_line = 0; i_line < number_of_lines; i_line++){

                    fprintf(f_out,"(( ");
                    fprintf(f_out,"%f %f %f", start_position_coords[i_line*3], start_position_coords[i_line*3+1], start_position_coords[i_line*3+2]);
                    fprintf(f_out," 0. 0. 0. 1.e-6 273. 1.e-10 ))\n");
                }

                free(start_position_coords);
            }

            fclose(f_out);

            Message0("\n...rCFD_write_Tracer_Positions: wrote %d positions to %s ...\n",
                total_number_of_start_postions, File_Dict.tracer_start_position_filename);
        }
#endif
    }
}

/*************************************************************************************/
DEFINE_DPM_INJECTION_INIT(rcfd_init_tracers,I)
/*************************************************************************************/
{
#if RP_NODE
    Particle    *p;
    Domain      *d=Get_Domain(1);
    Thread      *t = NULL, *t_phase = NULL;
    cell_t      c;

    int         i_phase;

    /* A. init particle vars */
    {
        int number_of_initialized_particles = 0;

        int i_tmp_0 = 0, i_tmp_1 = 0;

        loop(p,I->p_init){

            p->user[p_just_started] =   1.;
            p->user[p_start_time] =     CURRENT_TIME - CURRENT_TIMESTEP;
            p->user[p_just_killed] =    0.;

            i_phase = (number_of_initialized_particles % Solver_Dict.number_of_phases);

            if(Solver_Dict.number_of_phases == 1){

                p->user[p_phase_fraction] = 1.0;
            }
            else{

                c = P_CELL(p);

                t = P_CELL_THREAD(p);

                t_phase = THREAD_SUB_THREAD(t, i_phase);

                p->user[p_phase_fraction] = C_VOF(c, t_phase);
            }

            p->user[p_phase_id] = (double)i_phase;

            if(p->user[p_phase_id]<0.5) i_tmp_0++;
            else i_tmp_1++;

            number_of_initialized_particles++;
        }

        /*Message("\nmyid %d i_tmp_0 %d i_tmp_1 %d", myid, i_tmp_0, i_tmp_1);*/
    }


    /* B. upon first call, allocate Tracer.shifts */
    {

        int number_of_cells_in_partition;

        if(Tracer.allocated == 0){

            number_of_cells_in_partition = 0;

            thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){begin_c_loop_int(c,t){

                number_of_cells_in_partition++;

            }end_c_loop_int(c,t)}}

            /* allocate max. number of shifts,
               multiply by 2  because of tracer splitting in slow cells
               by default, each phase has the same number_of_shifts */

            loop_phases{

                Tracer.number_of_shifts[i_phase] = 2. * (node_last + 1) * number_of_cells_in_partition * Tracer_Dict.number_of_Tracers_per_cell;

                Tracer.shifts[i_phase] = (C2C_shift_type*)malloc( Tracer.number_of_shifts[i_phase] * sizeof(C2C_shift_type));
            }

            Tracer.allocated = 1;

            /*Message("\nmyid %d Tracer.number_of_shifts[0] %d Tracer.number_of_shifts[1] %d", myid, Tracer.number_of_shifts[0], Tracer.number_of_shifts[1]);*/

            Message0("\n... rCFD_init_tracers, allocated C2C shifts ...\n");
        }
        else{

            Message0("\n... rCFD_init_tracers, initialized particles  ...\n");
        }
    }
#endif
}

/*************************************************************************************/
DEFINE_DPM_SCALAR_UPDATE(rCFD_update_Tracers,i_cell,t,initialize,p)
/*************************************************************************************/
{
#if RP_NODE
    double  rand_real=0., time_ratio, random_walk_velocity;

    int     i_phase, i_tracer, i_layer = 0;

    Thread  *t_phase = NULL;

    /* U.0 Multiphase settings */
    {
        if(Solver_Dict.number_of_phases == 1){

            i_phase = 0;

            t_phase = t;
        }
        else{

            i_phase = (int)p->user[p_phase_id];

            t_phase = THREAD_SUB_THREAD(t, i_phase);
        }
    }

    /* U.1: Initialize Tracers */

    if(Tracer_just_started){

        rand_real = (double)rand()/(double)RAND_MAX;

        if((not_Tracer_start_interval) || (too_early_to_start_Tracers) || (Tracer_Database_full) || (kill_Tracer_by_coarse_graining)){

            p->user[p_just_killed] = 1.;
            p->state.mass = 0.;
            p->type = 2;        /* kill tracer (evaporates particle) */
        }
        else{

            Tracer.monitoring_started = 1;

            p->user[p_c0]       = (double)i_cell;
            p->user[p_node0]    = (double)myid;
            p->user[p_w0]       = C_VOLUME(i_cell,t) * p->user[p_phase_fraction];
            p->user[p_c_old]    = (double)i_cell;
            p->state.rho        = C_R(i_cell, t_phase);

            /* tracer_splitting in slow cells */
            if(excess_Tracer_cell_crossing_time){

                time_ratio = _C.crossing_time[i_phase][i_cell] / (2. * Phase_Dict[i_phase].time_step);  /* (time_ratio > 1.) */

                p->user[p_time_ratio] = time_ratio;

                /* store c0_c0_tracer */
                if(Tracer_Database_not_full){

                    i_tracer = Tracer.monitor_counter[i_phase];

                    if(i_tracer < (Tracer.number_of_shifts[i_phase] - 1)){

                        Tracer.shifts[i_phase][i_tracer].c0 =       i_cell;
                        Tracer.shifts[i_phase][i_tracer].node0 =    myid;
                        Tracer.shifts[i_phase][i_tracer].w0 =       (1.-1./time_ratio) * p->user[p_w0];
                        Tracer.shifts[i_phase][i_tracer].c1 =       i_cell;
                        Tracer.shifts[i_phase][i_tracer].node1 =    myid;

                        Tracer.monitor_counter[i_phase]++;
                    }
                    else{

                        Message("\n.. WARNING - Lost Tracer in partition %d\n", myid);
                    }
                }

                /* init c0_c1_tracer */
                {
                    p->user[p_w0] = (1./time_ratio) * p->user[p_w0];

                    p->state.V[0] = time_ratio * C_U(i_cell,t_phase);
                    p->state.V[1] = time_ratio * C_V(i_cell,t_phase);
                    p->state.V[2] = time_ratio * C_W(i_cell,t_phase);
                }
            }
            else{

                /* normal tracer initialization */
                {
                    p->state.V[0] = C_U(i_cell,t_phase);
                    p->state.V[1] = C_V(i_cell,t_phase);
                    p->state.V[2] = C_W(i_cell,t_phase);

                    p->user[p_time_ratio] = 1.0;
                }
            }

            if(Tracer_Dict.random_walk[i_phase]){

                random_walk_velocity = rCFD_user_set_random_walk_velocity();

                rand_real =         2.*((double)rand()/(double)RAND_MAX-0.5);
                p->user[p_u_rwm] =  rand_real * random_walk_velocity;

                rand_real =         2.*((double)rand()/(double)RAND_MAX-0.5);
                p->user[p_v_rwm] =  rand_real * random_walk_velocity;

                rand_real =         2.*((double)rand()/(double)RAND_MAX-0.5);
                p->user[p_w_rwm] =  rand_real * random_walk_velocity;

                p->state.V[0] += p->user[p_u_rwm];
                p->state.V[1] += p->user[p_v_rwm];
                p->state.V[2] += p->user[p_w_rwm];

                p->user[p_vel_rwm_old] = random_walk_velocity;
            }
            else{

                p->user[p_u_rwm] = 0.0;
                p->user[p_v_rwm] = 0.0;
                p->user[p_w_rwm] = 0.0;
            }
        }

        p->user[p_just_started] = 0.;
    }


    /* U.2: Consider (i) RWM changes across cells  (ii) fast Tracer from slow cells */

    if(Tracer_has_crossed_cell_border){

        p->user[p_c_old] = (double)i_cell;

        if(Tracer_Dict.random_walk[i_phase]){

            random_walk_velocity = rCFD_user_set_random_walk_velocity();

            p->user[p_u_rwm] *= random_walk_velocity / p->user[p_vel_rwm_old];
            p->user[p_v_rwm] *= random_walk_velocity / p->user[p_vel_rwm_old];
            p->user[p_w_rwm] *= random_walk_velocity / p->user[p_vel_rwm_old];

            p->user[p_vel_rwm_old] = random_walk_velocity;
        }

        if((Tracer_from_slow_cell) && (Tracer_Database_not_full) && (Tracer_not_stored_yet)){

            p->user[p_just_killed] = 1.;
            p->state.mass = 0.;
            p->type = 2;

            i_tracer = Tracer.monitor_counter[i_phase];

            if(i_tracer < (Tracer.number_of_shifts[i_phase] - 1)){

                Tracer.shifts[i_phase][i_tracer].c0 =       (int)p->user[p_c0];
                Tracer.shifts[i_phase][i_tracer].node0 =    (int)p->user[p_node0];
                Tracer.shifts[i_phase][i_tracer].w0 =       p->user[p_w0];
                Tracer.shifts[i_phase][i_tracer].c1 =       i_cell;
                Tracer.shifts[i_phase][i_tracer].node1 =    myid;

                Tracer.monitor_counter[i_phase]++;
            }
            else{

                Message("\n.. WARNING - Lost Tracer in partition %d\n", myid);
            }
        }
    }


    /* U.3: store and kill tracers after lifetime */

    if((Tracer_lifetime_expired) && (Tracer_Database_not_full) && (Tracer_not_stored_yet)){

        p->user[p_just_killed] = 2.;
        p->state.mass = 0.;
        p->type = 2;

        i_tracer = Tracer.monitor_counter[i_phase];

        if(i_tracer < (Tracer.number_of_shifts[i_phase] - 1)){

            Tracer.shifts[i_phase][i_tracer].c0 =       (int)p->user[p_c0];
            Tracer.shifts[i_phase][i_tracer].node0 =    (int)p->user[p_node0];
            Tracer.shifts[i_phase][i_tracer].w0 =       p->user[p_w0];
            Tracer.shifts[i_phase][i_tracer].c1 =       i_cell;
            Tracer.shifts[i_phase][i_tracer].node1 =    myid;

            Tracer.monitor_counter[i_phase]++;
        }
        else{

            Message("\n.. WARNING - Lost Tracer in partition %d\n", myid);
        }

        Tracer.ready2write = 1;
    }

#endif
}

/*************************************************************************************/
DEFINE_DPM_DRAG(rcfd_no_standard_drag,p,i)
/*************************************************************************************/
{
    return 0.;
}

/*************************************************************************************/
DEFINE_DPM_BODY_FORCE(rCFD_guide_Tracers,p,i)
/*************************************************************************************/
{
    Thread* t = P_CELL_THREAD(p);
    cell_t  c = P_CELL(p);

    if(Solver_Dict.number_of_phases == 1){

        if(i==0)    return p->user[p_time_ratio] * (C_U(c,t) + p->user[p_u_rwm]) - p->state.V[0];

        if(i==1)    return p->user[p_time_ratio] * (C_V(c,t) + p->user[p_v_rwm]) - p->state.V[1];

        if(i==2)    return p->user[p_time_ratio] * (C_W(c,t) + p->user[p_w_rwm]) - p->state.V[2];
    }
    else{

        Thread  *t_phase = NULL;

        int i_phase = (int)p->user[p_phase_id];

        t_phase = THREAD_SUB_THREAD(t, i_phase);

        if(i==0)    return p->user[p_time_ratio] * (C_U(c,t_phase) + p->user[p_u_rwm]) - p->state.V[0];

        if(i==1)    return p->user[p_time_ratio] * (C_V(c,t_phase) + p->user[p_v_rwm]) - p->state.V[1];

        if(i==2)    return p->user[p_time_ratio] * (C_W(c,t_phase) + p->user[p_w_rwm]) - p->state.V[2];
    }

    return 0.;
}

/*************************************************************************************/
DEFINE_EXECUTE_AT_END(rCFD_write_C2Cs)
/*************************************************************************************/
{
#if RP_NODE

    /* TODO (9/21) - adapt weights of MPI C2Cs such that they obey to actual mass fluxes */

    if((Tracer.ready2write == 1) && (Tracer_Database_not_full)){

        int     i_phase, i_layer = 0;

        loop_phases{

            if(Tracer.monitor_counter[i_phase] > 0){

                Domain  *d = Get_Domain(1);
                Thread  *t;

                FILE    *f_out = NULL;
                char    file_name[80];

                int     i_state, i_cell, i_tracer;

                int     c0, node0, c1, node1, number_of_lock_cells;

                i_state = Solver.current_state;

                sprintf(file_name,"%s_%d_%d_%d", File_Dict.C2C_filename, i_state, i_phase, myid);

                if(Tracer.frame_counter == 0){

                    f_out = fopen(file_name,"w");
                }
                else{

                    f_out = fopen(file_name,"a");
                }

                if(f_out == NULL){

                    Message("\n... ERROR: rCFD_write_C2Cs: f_out == NULL  ...\n");
                }
                else{

                    /* avoid lock cells (i.e. cells which only point to themselves) */
                    {
                        thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){

                            begin_c_loop(i_cell,t){

                                _C.hit_by_other_cell[i_cell] = 0;

                            }end_c_loop(i_cell,t);
                        }}

                        for(i_tracer = 0; i_tracer < Tracer.monitor_counter[i_phase]; i_tracer++){

                            c0 =        Tracer.shifts[i_phase][i_tracer].c0;
                            node0 =     Tracer.shifts[i_phase][i_tracer].node0;
                            c1 =        Tracer.shifts[i_phase][i_tracer].c1;
                            node1 =     Tracer.shifts[i_phase][i_tracer].node1;

                            if((c0 != c1) || (node0 != node1)){

                                _C.hit_by_other_cell[c1] += 1;
                            }
                        }

                        /* loop Tracers and count lock cells */

                        number_of_lock_cells = 0;

                        for(i_tracer = 0; i_tracer < Tracer.monitor_counter[i_phase]; i_tracer++){

                            c1 = Tracer.shifts[i_phase][i_tracer].c1;

                            if(_C.hit_by_other_cell[c1] == 0){

                                number_of_lock_cells++;
                            }
                        }

                        Message0("\n... rCFD_write_C2Cs: identified %d lock cells\n", PRF_GISUM1(number_of_lock_cells));
                    }

                    fprintf(f_out,"%d \n", (Tracer.monitor_counter[i_phase] - number_of_lock_cells));

                    for(i_tracer = 0; i_tracer < Tracer.monitor_counter[i_phase]; i_tracer++){

                        c1 = Tracer.shifts[i_phase][i_tracer].c1;

                        if(_C.hit_by_other_cell[c1] > 0){

                            fprintf(f_out,"%d ",Tracer.shifts[i_phase][i_tracer].c0);
                            fprintf(f_out,"%d ",Tracer.shifts[i_phase][i_tracer].node0);
                            fprintf(f_out,"%d ",Tracer.shifts[i_phase][i_tracer].c1);
                            fprintf(f_out,"%d ",Tracer.shifts[i_phase][i_tracer].node1);
                            fprintf(f_out,"%e ",Tracer.shifts[i_phase][i_tracer].w0);

                            fprintf(f_out,"\n");
                        }
                    }

                    Tracer.monitor_counter[i_phase] = 0;

                    fclose(f_out);

                    Message0("\n... wrote C2C frame # %d for phase %d ...", Tracer.frame_counter, i_phase);
                }
            }
        }

        Tracer.frame_counter++;

        Tracer.ready2write = 0;
    }
    else{

       Message0("\n... no Tracer frame written; already wrote %d frames ...", Tracer.frame_counter);
    }
#endif
}

/*************************************************************************************/
DEFINE_EXECUTE_AT_END(rCFD_write_Norms)
/*************************************************************************************/
{
#if RP_NODE

    Domain  *d=Get_Domain(1);
    Thread  *t, *t_mix, *t_phase;

    int     i_state, i_norm, i_phase, i_cell, i_layer = 0;

    if ((Tracer.monitoring_started) && (Norms_write_interval) && (Norm_Database_not_full)){

        /* N.1 calc norms */
        {
            i_norm = 0;

            t_mix = NULL;

            if(Norm_Dict.format == standard){

                thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){begin_c_loop_int(i_cell,t){

                    if((i_cell % Norm_Dict.coarse_graining) == 0){

                        if(Solver_Dict.number_of_phases == 1){

                            Norms.norm[i_norm] = sqrt(C_U(i_cell,t)*C_U(i_cell,t) + C_V(i_cell,t)*C_V(i_cell,t) + C_W(i_cell,t)*C_W(i_cell,t));
                        }
                        else{

                            t_mix = t;

                            t_phase = THREAD_SUB_THREAD(t, 0);

                            Norms.norm[i_norm] = C_VOF(i_cell, t_phase);
                        }

                        i_norm++;
                    }

                }end_c_loop_int(i_cell,t)}}

            }
            else{

                rCFD_user_set_Norm();
            }
        }

        /* N.2 write norms */
        {
            FILE    *f_out;
            char    file_name[80];

            i_state = Solver.current_state;

            sprintf(file_name,"%s_%d_%d", File_Dict.Norm_filename, i_state, myid);

            if(Norms.frame_counter == 0){

                f_out = fopen(file_name, "w");
            }
            else{

                f_out = fopen(file_name, "a");
            }

            if(f_out == NULL){

                Message("\n... ERROR: rCFD_write_Norms: f_out == NULL ...\n");

                return;
            }

            fprintf(f_out,"%d \n", Norms.number_of_norms);

            loop_norms{

                fprintf(f_out,"%e \n", Norms.norm[i_norm]);
            }

            fclose(f_out);

            Message0("\n... Norm frame # %d written ...\n", Norms.frame_counter);
        }

        /* N.3 in case of multiphase, write vof Files */
        if(Solver_Dict.number_of_phases > 1){

            FILE    *f_out;
            char    file_name[80];

            i_state = Solver.current_state;

            loop_phases{

                sprintf(file_name,"%s_%d_%d_%d", File_Dict.vof_filename, i_state, i_phase, myid);

                if(Norms.frame_counter == 0){

                    f_out = fopen(file_name, "w");
                }
                else{

                    f_out = fopen(file_name, "a");
                }

                if(f_out == NULL){

                    Message("\n... ERROR: rCFD_write_Norms: vof: f_out == NULL ...\n");

                    return;
                }

                fprintf(f_out,"%d \n", _Cell_Dict.number_of_cells);

                t_phase = THREAD_SUB_THREAD(t_mix, i_phase);

                if(t_phase == NULL){

                    Message("\n... ERROR: rCFD_write_Norms: t_phase == NULL for i_phase = %d\n", i_phase);

                    fclose(f_out);

                    return;
                }

                loop_cells{

                    fprintf(f_out,"%e \n", C_VOF(i_cell, t_phase));
                }

                fclose(f_out);

                Message0("\n... vof frame # %d written for phase %d ...", Norms.frame_counter, i_phase);
            }
        }

        Norms.frame_counter++;
    }
    else{

        Message0("\n... no Norm frame written; already wrote %d frames...", Norms.frame_counter);
    }
#endif
}

/*************************************************************************************/
DEFINE_ON_DEMAND(rCFD_write_Rec)
/*************************************************************************************/
{
#if RP_NODE

    double  **Norm_Database;

    int     i_state, i_norm, i_cell, i_tmp, i_layer = 0;

    FILE    *f_in = NULL;

#endif

    int     i_frame, i_frame2, i_island;

    double  norm, diff;

    char    file_name[80];

#if RP_HOST

    double  Rec_matrix[Solver_Dict.number_of_frames][Solver_Dict.number_of_frames];

    int     Rec_jump[Solver_Dict.number_of_frames];

    FILE    *f_out = NULL;
#endif


    /* read Norms */
    {
#if RP_NODE

        Norm_Database = (double**)malloc(Solver_Dict.number_of_frames * sizeof(double*));

        loop_frames{

            Norm_Database[i_frame] = (double*)malloc(Norms.number_of_norms * sizeof(double));
        }

        i_state = Solver.current_state;

        sprintf(file_name,"%s_%d_%d", File_Dict.Norm_filename, i_state, myid);

        f_in = fopen(file_name,"r");

        if(f_in == NULL){

            Message("\n... ERROR: rCFD_read_Norms: f_in == NULL ...\n");
        }
        else{

            loop_frames{

                fscanf(f_in, "%d \n", &i_tmp);

                loop_norms{

                    fscanf(f_in, "%le \n", &Norm_Database[i_frame][i_norm]);
                }
            }

            fclose(f_in);
        }
#endif
    }

    loop_islands{

        /* calc Rec_Matrix */
        {
#if RP_NODE
            double  mean_norm, sum_of_all_norms = 0.0;

            loop_frames{

                loop_norms{

                    sum_of_all_norms += Norm_Database[i_frame][i_norm];
                }
            }

            mean_norm = PRF_GRSUM1(sum_of_all_norms) /

                (double) PRF_GISUM1(Norms.number_of_norms) /

                (double) Solver_Dict.number_of_frames;
#endif

            loop_frames{

                loop_frames2{
#if RP_NODE
                    norm = 0.0;

                    loop_norms{

                        i_cell = i_norm * Norm_Dict.coarse_graining;

                        if( _C.island_id[i_cell] == i_island){

                            diff = fabs(Norm_Database[i_frame][i_norm] - Norm_Database[i_frame2][i_norm]);

                            if(Norm_Dict.larger_than_mean_norm_only){

                                if(diff > mean_norm){

                                    norm += diff;
                                }
                            }
                            else{

                                norm += diff;
                            }
                        }
                    }

                    norm = PRF_GRSUM1(norm);
#endif
                    node_to_host_real_1(norm);
#if RP_HOST
                    Rec_matrix[i_frame][i_frame2] = norm;
                    Rec_matrix[i_frame2][i_frame] = norm;
#endif
                }
            }
        }

        /* calc jumps */
        {
#if RP_HOST
            if(Rec_Dict.method == quarter_jumps){

                int N1_4=(int)floor((double)Solver_Dict.number_of_frames/4.);
                int N2_4=(int)floor((double)Solver_Dict.number_of_frames/2.);
                int N3_4=(int)floor((double)Solver_Dict.number_of_frames*3./4.);

                int j_frame, j_frame_min = 0;

                loop_frames{

                    diff = 1e10; /* something large */

                    if(i_frame <= N1_4){

                        for(j_frame = N2_4; j_frame < N3_4; j_frame++){

                            if (Rec_matrix[i_frame][j_frame] < diff){

                                diff = Rec_matrix[i_frame][j_frame];

                                j_frame_min = j_frame;
                            }
                        }
                    }

                    else if(i_frame <= N2_4){

                        for(j_frame = N3_4; j_frame < Solver_Dict.number_of_frames; j_frame++){

                            if (Rec_matrix[i_frame][j_frame] < diff){

                                diff = Rec_matrix[i_frame][j_frame];

                                j_frame_min = j_frame;
                            }
                        }
                    }

                    else if(i_frame <= N3_4){

                        for(j_frame = 0; j_frame < N1_4; j_frame++){

                            if (Rec_matrix[i_frame][j_frame] < diff){

                                diff = Rec_matrix[i_frame][j_frame];

                                j_frame_min = j_frame;
                            }
                        }
                    }

                    else{

                        for(j_frame = N1_4; j_frame < N2_4; j_frame++){

                            if (Rec_matrix[i_frame][j_frame] < diff){

                                diff = Rec_matrix[i_frame][j_frame];

                                j_frame_min = j_frame;
                            }
                        }
                    }

                    Rec_jump[i_frame] = j_frame_min;
                }
            }
#endif
        }

        /* write matrix, jumps */
        {
#if RP_HOST
            /* Rec Matrix */
            {
                sprintf(file_name, "%s_%d_%d_%d", File_Dict.Matrix_filename, Solver.current_state, Solver.current_state, i_island);

                f_out = fopen(file_name, "w");

                if(f_out == NULL){

                    Message("\n... ERROR: rCFD_write_Rec: f_out == NULL , Rec_Matrix Filename: %s ...\n", file_name);

                    return;
                }

                loop_frames{

                    loop_frames2{

                        fprintf(f_out, "%e ", Rec_matrix[i_frame][i_frame2]);
                    }

                    fprintf(f_out, "\n");
                }

                Message("\n...rCFD_write_Rec: wrote %s ...\n", file_name);

                fclose(f_out);
            }

            /* Rec Path */
            {
                sprintf(file_name, "%s_%d_%d_%d", File_Dict.Jump_filename, Solver.current_state, Solver.current_state, i_island);

                f_out = fopen(file_name, "w");

                if(f_out == NULL){

                    Message("\n... ERROR: rCFD_write_Rec: f_out == NULL , Rec_Path Filename: %s ...\n", file_name);

                    return;
                }

                loop_frames{

                    fprintf(f_out, "%d \n", Rec_jump[i_frame]);
                }

                Message("\n...rCFD_write_Rec: wrote %s ...\n", file_name);

                fclose(f_out);
            }
#endif
        }

    }

    /* free local vars */
    {
#if RP_NODE
        if(Norm_Database != NULL){

            loop_frames{

                if(Norm_Database[i_frame] != NULL){

                    free(Norm_Database[i_frame]);
                }
            }

            free(Norm_Database);
        }
#endif
    }

    Message0("\n\n...rCFD_write_Rec");
}

/*************************************************************************************/
DEFINE_ON_DEMAND(rCFD_free_all)
/*************************************************************************************/
{
    free_all();

#if RP_NODE
    free_all_parallel();
#endif

    Message0("\n\n...rCFD_free_all");
}
