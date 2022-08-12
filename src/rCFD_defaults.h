#ifndef RCFD_DEFAULTS
#define RCFD_DEFAULTS

#include "udf.h"
#include <sys/stat.h>
#include "rCFD_macros.h"
#include "rCFD_globals.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

/* (C)  2021-22
    Stefan Pirker
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/

    #define     Global_Version_Year     2022
    #define     Global_Version_Month    2

	enum{
		no_verbose,
		low_verbose,
		high_verbose
	};

    #define     Global_Verbose          low_verbose
	

    /* A. default Dicts */

    void rCFD_default_Solver_Dict(void)
    {
        Solver_Dict.version_year =     Global_Version_Year;
        Solver_Dict.version_month =    Global_Version_Month;

        Solver_Dict.verbose =          Global_Verbose;

        Solver_Dict.number_of_frames =     1;
        Solver_Dict.number_of_states =     1;
        Solver_Dict.number_of_phases =     1;
        Solver_Dict.number_of_islands =    1;
        Solver_Dict.number_of_layers =     1;
        Solver_Dict.number_of_runs =       1;

        Solver_Dict.recurrence_process_on =    1;
        Solver_Dict.data_convection_on =       1;
        Solver_Dict.face_diffusion_on =        1;
        Solver_Dict.data_binarization_on =     0;
        Solver_Dict.data_drifting_on =         0;
        Solver_Dict.balance_correction_on =    1;
        Solver_Dict.on_the_fly_post_on =       1;

        Solver_Dict.analyse_CFD_count = 0;

        Solver_Dict.max_number_of_cells_per_time_step = 10;

        Solver_Dict.time_steps_per_monitoring_interval =   10;
        Solver_Dict.start_time_for_monitoring =            0.0;

        Solver_Dict.min_number_of_cells_per_layer =     100;
        Solver_Dict.max_number_of_faces_per_cell =      100;
        Solver_Dict.max_number_of_children =            100;
        Solver_Dict.number_redist_loops_per_layer =     5;

        Solver_Dict.C2C_loading_reduction =    1;

        Solver_Dict.global_time_step = 1.0;

        Solver_Dict.max_fill_loops =   10;

        Solver_Dict.face_swap_max_per_loop =   (1./8.);
        Solver_Dict.face_swap_min =            (1./64.);
        Solver_Dict.face_swap_max_loops =      4;

        Solver_Dict.number_of_drift_loops =    1;

        Solver_Dict.balance_correction_update = 1;
        Solver_Dict.control_conc_sum_on = 1;

        if(Transcript){

            FILE    *f_out = fopen("./Run.trn", "w");

            if(f_out){

                fprintf(f_out,"rCFD_default_Solver_Dict");

                fprintf(f_out,"\n\n   Version: %4d.%02d", Solver_Dict.version_year, Solver_Dict.version_month);

                time_t      current_time = time(NULL);

                char        *c_time_string = ctime(&current_time);

                fprintf(f_out,"\n\n   Run start @ %s", c_time_string);

                fclose(f_out);
            }
        }
    }

    void rCFD_default_File_Dict(void)
    {

        /* Node-0 generates folder structure */
        if(myid == 0){
			
			short directory_error_indicator = 0;

            directory_error_indicator += mkdir("./data",0777);
            directory_error_indicator += mkdir("./data/tmp",0777);
            directory_error_indicator += mkdir("./data/vof",0777);
            directory_error_indicator += mkdir("./data/c2c",0777);
            directory_error_indicator += mkdir("./rec",0777);
            directory_error_indicator += mkdir("./post",0777);

            if(directory_error_indicator != 0){
				 
				Message0("\nWARNING rCFD_default_File_Dict: could not create all sub-directories (maybe because they already exist)\n");
            }
        }

        File_Dict.tracer_start_position_filename   =   "./data/tmp/tracer_start_pos.inj";

        File_Dict.C2C_filename =                       "./data/c2c/c2c";

        File_Dict.Norm_filename =                      "./data/tmp/norm";

        File_Dict.vof_filename =                       "./data/vof/vof";

        File_Dict.Jump_filename =                      "./rec/jump";

        File_Dict.Matrix_filename =                    "./rec/matrix";

        File_Dict.Balance_filename =                   "./post/balance_monitor.out";
    }

    void rCFD_default_Phase_Dict(void)
    {
#if RP_NODE
        int i_phase;

        for(i_phase = 0; i_phase < Solver_Dict.number_of_phases; i_phase++){

            Phase_Dict[i_phase].number_of_data = 1;

            Phase_Dict[i_phase].time_step = Solver_Dict.global_time_step;

            Phase_Dict[i_phase].shift_probability = 1.0;

            Phase_Dict[i_phase].density = 1.0;

            Phase_Dict[i_phase].heat_capacity = 1.0;

            Phase_Dict[i_phase].reference_temperature = 300.0;

            Phase_Dict[i_phase].number_of_user_vars = 0;

            Phase_Dict[i_phase].user = NULL;
        }
#endif
    }

    void rCFD_default_Tracer_Dict(void)
    {
#if RP_NODE

        Tracer_Dict.format = guide_by_force_format;

        Tracer_Dict.number_of_Tracers_per_cell = 1;

        Tracer_Dict.region_of_interest_exists = 0;

        Tracer_Dict.ROI_x_min = 0.0;
        Tracer_Dict.ROI_x_max = 0.0;
        Tracer_Dict.ROI_y_min = 0.0;
        Tracer_Dict.ROI_y_max = 0.0;
        Tracer_Dict.ROI_z_min = 0.0;
        Tracer_Dict.ROI_z_max = 0.0;

        Tracer_Dict.coarse_graining = 1;

        int i_phase;

        loop_phases{

            Tracer_Dict.random_walk[i_phase] = 0;
        }

        Tracer_Dict.C2C_format = c0_n0_c1_n1_w0;
#endif
    }

    void rCFD_default_Norm_Dict(void)
    {
#if RP_NODE
        Norm_Dict.format = standard;

        Norm_Dict.coarse_graining = 1;

        Norm_Dict.larger_than_mean_norm_only = 0;
#endif
    }

    void rCFD_default_Rec_Dict(void)
    {

        Rec_Dict.format = quarter_jumps_format;

        Rec_Dict.min_seq_length = (int)((double)Solver_Dict.number_of_frames/25.);

        if(Rec_Dict.min_seq_length < 1){

            Rec_Dict.min_seq_length = 1;
        }

        Rec_Dict.max_seq_length = (int)((double)Solver_Dict.number_of_frames/4.);

        if(Rec_Dict.max_seq_length < 1){

            Rec_Dict.max_seq_length = 1;
        }

        Rec_Dict.off_diagonal_band_width = (int)((double)Solver_Dict.number_of_frames/10.);

        if(Rec_Dict.off_diagonal_band_width < 1){

            Rec_Dict.off_diagonal_band_width = 1;
        }
    }

    void rCFD_default_Data_Dict(void)
    {
#if RP_NODE
        int i_phase, i_data;

        loop_phases{

            loop_data{

                Data_Dict[i_phase][i_data].type = generic_data;

                Data_Dict[i_phase][i_data].physical_diff = 1./8.;

                Data_Dict[i_phase][i_data].binarization_art_diff = 1./8.;

                Data_Dict[i_phase][i_data].drifting_data = 0;

                Data_Dict[i_phase][i_data].drift_velocity = NULL;

                Data_Dict[i_phase][i_data].number_of_user_vars = 0;

                Data_Dict[i_phase][i_data].user = NULL;
            }
        }
#endif
    }

    void rCFD_default_Balance_Dict(void)
    {
#if RP_NODE
        int i_phase, i_data;

        loop_phases{

            loop_data{

                Balance_Dict[i_phase][i_data].type = global_balancing;

                Balance_Dict[i_phase][i_data].max_correction_loops = 10;

                Balance_Dict[i_phase][i_data].accuracy_level = 0.01;

                Balance_Dict[i_phase][i_data].write_balance_to_file = 0;

                Balance_Dict[i_phase][i_data].write_balance_to_file_interval = 1;
            }
        }
#endif
    }

    void rCFD_default_Topo_Dict(void)
    {
        Topo_Dict.Cell_Dict = NULL;

        Topo_Dict.Face_Dict = NULL;
    }

    void rCFD_default_Cell_Dict_L0(void)
    {
        int i_layer = 0;

        Domain  *d=Get_Domain(1);
        Thread  *t;
        cell_t  c;

        int     number_of_local_cells = 0, number_of_local_internal_cells = 0;

        int     only_one_fluid_cell_thread = 0;

        thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){

            only_one_fluid_cell_thread++;

            begin_c_loop_all(c,t){

                number_of_local_cells++;

            }end_c_loop_all(c,t)

            begin_c_loop_int(c,t){

                number_of_local_internal_cells++;

            }end_c_loop_int(c,t)

        }}

        only_one_fluid_cell_thread = PRF_GIHIGH1(only_one_fluid_cell_thread);

        if(only_one_fluid_cell_thread > 1){

            Message0("\n... WARNING: found more than one fluid cell thread ...\n");
        }

        _Cell_Dict.number_of_cells = number_of_local_cells;

        _Cell_Dict.number_of_int_cells = number_of_local_internal_cells;

        _Cell_Dict.number_of_ext_cells = number_of_local_cells - number_of_local_internal_cells;

        _Cell_Dict.number_of_user_vars = 0;
    }

    void rCFD_default_Face_Dict_L0(void)
    {
        int i_layer = 0;

        Domain  *d=Get_Domain(1);
        Thread  *t;
        face_t  f;
        cell_t  c0, c1;

        int     number_of_local_int_faces = 0, number_of_local_ext_faces = 0;

        thread_loop_f(t,d){if(THREAD_TYPE(t)==THREAD_F_INTERIOR){

            begin_f_loop(f,t){

                c0 = (int)F_C0(f, t);

                c1 = (int)F_C1(f, t);

                if((c0 >= 0)&&(c1 >= 0)){

                    number_of_local_int_faces++;
                }

            }end_f_loop(f,t)

            begin_f_loop_ext(f,t){

                c0 = (int)F_C0(f, t);

                c1 = (int)F_C1(f, t);

                if((c0 >= 0)&&(c1 >= 0)){

                    number_of_local_ext_faces++;
                }

            }end_f_loop_ext(f,t)
        }}

        _Face_Dict.number_of_int_faces = number_of_local_int_faces;

        _Face_Dict.number_of_ext_faces = number_of_local_ext_faces;

        _Face_Dict.number_of_faces = (number_of_local_int_faces + number_of_local_ext_faces);
    }

    /* B. default vars settings */

    void rCFD_default_Solver(void)
    {
        Solver.current_state = 0;

        Solver.current_layer = 0;

        Solver.global_run_counter = 0;

        Solver.timestep_width = (double*)malloc(Solver_Dict.number_of_layers * sizeof(double));

        int i_layer;

        loop_layers{

            Solver.timestep_width[i_layer] = Solver_Dict.global_time_step * pow(2.0, (double) i_layer);
        }

        Solver.global_time = 0.0;
    }

    void rCFD_default_Topo(void)
    {
#if RP_NODE
        Topo.Cell = NULL;

        Topo.Face = NULL;
#endif
    }

    void rCFD_default_Cell_L0(void)
    {
#if RP_NODE
        int i_layer = 0;

        int i_phase, i_cell, i_frame, i_user, i_data, i_dim;

        Domain  *d=Get_Domain(1);
        Thread  *t;
        double  x[3];

        thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){begin_c_loop_all(i_cell,t){

                C_CENTROID(x, i_cell, t);

                loop_dim{

                    _C.x[i_cell][i_dim] = x[i_dim];
                }

                _C.volume[i_cell] = C_VOLUME(i_cell, t);

        }end_c_loop_all(i_cell,t)}}

        loop_phases{

            loop_cells{

                _C.average_velocity[i_phase][i_cell] =  0.0;
                _C.crossing_time[i_phase][i_cell] =     0.0;
            }
        }

        loop_cells{

            _C.hit_by_other_cell[i_cell] = 0;
            _C.island_id[i_cell] = 0;

            _C.weight_after_shift[i_cell] = 0.0;
            _C.weight_after_swap[i_cell] =  0.0;
        }

        if(_C.vof != NULL){

            FILE    *f_in = NULL;
            char    file_name[80];

            int     i_state, i_tmp;

            if(Solver_Dict.number_of_phases == 1){

                i_phase = 0;

                loop_frames{

                    loop_cells{

                        _C.vof[_i_vof] = 1.0;
                    }
                }
            }
            else{

				short	vof_file_exisits = 1;
				
                i_state = 0;

                loop_phases{

                    sprintf(file_name,"%s_%d_%d_%d", File_Dict.vof_filename, i_state, i_phase, myid);

                    f_in = fopen(file_name, "r");

                    if(f_in == NULL){

                        loop_frames{

                            loop_cells{

                                if(i_phase == 0){

                                    _C.vof[_i_vof] = 1.0;
                                }
                                else{

                                    _C.vof[_i_vof] = 0.0;
                                }
                            }
                        }

						vof_file_exisits = 0;						
                    }
                    else{
						
                        loop_frames{

                            fscanf(f_in,"%d\n", &i_tmp);

                            if(i_tmp != _Cell_Dict.number_of_cells){

                                Message("\nERROR: rCFD_default_Cell: _C.vof: i_tmp != _Cell_Dict.number_of_cells ...\n");

                                return;
                            }

                            loop_cells{

                                fscanf(f_in,"%le\n", &_C.vof[_i_vof]);

                                if((_C.vof[_i_vof] < 0.0) || (_C.vof[_i_vof] > 1.0)){

                                    Message("\nERROR myid %d i_frame %d i_cell %d i_phase %d vof %e", myid, i_frame, i_cell, i_phase, _C.vof[_i_vof]);

                                    return;
                                }
                            }
                        }

                        fclose(f_in);
                    }
                }
				
				vof_file_exisits = PRF_GILOW1(vof_file_exisits);
				
				if( ! vof_file_exisits){
					
					Message0("\nWARNING rCFD_default_Cell_L0: at least one vof file doesn't exist (for the preparation phase, this is normal)\n");
				}
            }
        }

        loop_phases{

            loop_cells{

                loop_data{

                    _C.data[i_phase][i_cell][i_data] =      0.0;
                    _C.data_shift[i_phase][i_cell][i_data] = 0.0;
                    _C.data_swap[i_phase][i_cell][i_data] =     0.0;
                }
            }
        }

        if(_C.drift_exchange != NULL){

            loop_cells{

                _C.drift_exchange[i_cell] = 0.0;
            }
        }

        if(_C.user != NULL){

            loop_cells{

                for(i_user = 0; i_user < Topo_Dict.Cell_Dict[i_layer].number_of_user_vars; i_user++){

                    _C.user[i_cell][i_user] = 0.0;
                }
            }
        }

#endif
    }

#if 0
    void rCFD_default_Cell(const short i_layer)
    {
#if RP_NODE
        if(i_layer == 0){

            int i_phase, i_cell, i_frame, i_user, i_data, i_dim;

            Domain  *d=Get_Domain(1);
            Thread  *t;
            double  x[3];

            thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){begin_c_loop_all(i_cell,t){

                    C_CENTROID(x, i_cell, t);

                    loop_dim{

                        C.x[i_cell][i_dim] = x[i_dim];
                    }

                    C.volume[i_cell] = C_VOLUME(i_cell, t);

            }end_c_loop_all(i_cell,t)}}

            loop_phases{

                loop_cells{

                    C.average_velocity[i_phase][i_cell] =  0.0;
                    C.crossing_time[i_phase][i_cell] =     0.0;
                }
            }

            loop_cells{

                C.hit_by_other_cell[i_cell] = 0;
                C.island_id[i_cell] = 0;

                C.weight_after_shift[i_cell] = 0.0;
                C.weight_after_swap[i_cell] =  0.0;
            }

            if(C.vof != NULL){

                FILE    *f_in = NULL;
                char    file_name[80];

                int     i_state = 0, i_tmp;

                loop_phases{

                    if(Solver_Dict.number_of_phases == 1){

                        loop_frames{

                            loop_int_cells{

                                C.vof[i_frame][i_cell][i_phase] = 1.0;
                            }
                        }
                    }
                    else{

                        sprintf(file_name,"%s_%d_%d_%d", File_Dict.vof_filename, i_state, i_phase, myid);

                        f_in = fopen(file_name, "r");

                        if(f_in == NULL){

                            loop_frames{

                                loop_int_cells{

                                    if(i_phase == 0){

                                        C.vof[i_frame][i_cell][i_phase] = 1.0;
                                    }
                                    else{

                                        C.vof[i_frame][i_cell][i_phase] = 0.0;
                                    }
                                }
                            }

                            Message0("\nWARNING: rCFD_default_Cell: C.vof: f_in == NULL for i_phase %d ...\n", i_phase);
                        }
                        else{

                            loop_frames{

                                fscanf(f_in,"%d\n", &i_tmp);

                                if(i_tmp != _Cell_Dict.number_of_int_cells){

                                    Message("\nERROR: rCFD_default_Cell: C.vof: i_tmp != _Cell_Dict.number_of_int_cells ...\n");

                                    return;
                                }

                                loop_int_cells{

                                    fscanf(f_in,"%le\n", &C.vof[_i_vof]);

                                    if((C.vof[_i_vof] < 0.0) || (C.vof[_i_vof] > 1.0)){

                                        Message("\nERROR myid %d i_frame %d i_cell %d i_phase %d vof %e", myid, i_frame, i_cell, i_phase, C.vof[_i_vof]);

                                        return;
                                    }
                                }
                            }

                            fclose(f_in);
                        }
                    }
                }
            }

            loop_phases{

                loop_cells{

                    loop_data{

                        C.data[i_phase][i_cell][i_data] =      0.0;
                        C.data_shift[i_phase][i_cell][i_data] = 0.0;
                        C.data_swap[i_phase][i_cell][i_data] =     0.0;
                    }
                }
            }

            if(C.drift_exchange != NULL){

                loop_cells{

                    C.drift_exchange[i_cell] = 0.0;
                }
            }

            if(C.user != NULL){

                loop_cells{

                    for(i_user = 0; i_user < Topo_Dict.Cell_Dict[i_layer].number_of_user_vars; i_user++){

                        C.user[i_cell][i_user] = 0.0;
                    }
                }
            }
        }
#endif
    }
#endif
    void rCFD_default_Face_L0(void)
    {
#if RP_NODE
        int i_layer = 0;

        Domain  *d=Get_Domain(1);
        Thread  *t;
        int     f, i_face, c0, c1, i_dim;

        double  A[3];

        i_face = 0;

        thread_loop_f(t,d){if(THREAD_TYPE(t)==THREAD_F_INTERIOR){

            begin_f_loop(f,t){

                c0 = (int)F_C0(f, t);

                c1 = (int)F_C1(f, t);

                if((c0 >= 0)&&(c1 >= 0)){

                    _F.c0[i_face] = c0;

                    _F.c1[i_face] = c1;

                    F_AREA(A, f, t);

                    loop_dim{

                        _F.area[i_face][i_dim] = A[i_dim];
                    }

                    i_face++;
                }

            }end_f_loop(f,t)

            begin_f_loop_ext(f,t){

                c0 = (int)F_C0(f, t);

                c1 = (int)F_C1(f, t);

                if((c0 >= 0)&&(c1 >= 0)){

                    _F.c0[i_face] = c0;

                    _F.c1[i_face] = c1;

                    F_AREA(A, f, t);

                    loop_dim{

                        _F.area[i_face][i_dim] = A[i_dim];
                    }

                    i_face++;
                }

            }end_f_loop_ext(f,t)

        }}
#endif
    }

#if 0
    void rCFD_default_Face(const short i_layer)
    {
#if RP_NODE
        if(i_layer == 0){

            Domain  *d=Get_Domain(1);
            Thread  *t;
            int     f, i_face, c0, c1, i_dim;

            double  A[3];

            i_face = 0;

            thread_loop_f(t,d){if(THREAD_TYPE(t)==THREAD_F_INTERIOR){

                begin_f_loop(f,t){

                    c0 = (int)F_C0(f, t);

                    c1 = (int)F_C1(f, t);

                    if((c0 >= 0)&&(c1 >= 0)){

                        F.c0[i_face] = c0;

                        F.c1[i_face] = c1;

                        F_AREA(A, f, t);

                        loop_dim{

                            F.area[i_face][i_dim] = A[i_dim];
                        }

                        i_face++;
                    }

                }end_f_loop(f,t)

                begin_f_loop_ext(f,t){

                    c0 = (int)F_C0(f, t);

                    c1 = (int)F_C1(f, t);

                    if((c0 >= 0)&&(c1 >= 0)){

                        F.c0[i_face] = c0;

                        F.c1[i_face] = c1;

                        F_AREA(A, f, t);

                        loop_dim{

                            F.area[i_face][i_dim] = A[i_dim];
                        }

                        i_face++;
                    }

                }end_f_loop_ext(f,t)

            }}
        }
#endif
    }
#endif
    void rCFD_default_Tracer(void)
    {
#if RP_NODE
        Tracer.allocated =     0;
        Tracer.initialized =   0;
        Tracer.ready2write =   0;
        Tracer.monitoring_started  = 0;

        Tracer.frame_counter =     0;

        int     i_phase;

        loop_phases{

            Tracer.monitor_counter[i_phase] =  0;
            Tracer.number_of_shifts[i_phase] = 0;

            Tracer.shifts[i_phase] = NULL;
        }
#endif
    }

    void rCFD_default_Norms(void)
    {
#if RP_NODE
        Norms.frame_counter = 0;
#endif
    }

    void rCFD_default_Rec(void)
    {
        int i_state, i_state2, i_island, i_frame;

        int number_of_jump_files = 0;

        loop_islands{

            Rec.global_frame[i_island] = 0;
        }

        Rec.frame_in_sequence = 0;

        Rec.sequence_length = Rec_Dict.min_seq_length;

        loop_states{

            loop_states2{

                loop_islands{

                    loop_frames{

                        Rec.jumps[i_state][i_state2][i_island][i_frame] = -1;
                    }

#if RP_HOST
                    char    filename[40];
                    FILE    *fi;

                    sprintf(filename,"%s_%d_%d_%d", File_Dict.Jump_filename, i_state, i_state2, i_island);

                    fi=fopen(filename,"r");

                    if(fi!=NULL){

                        fopen(filename,"r");

                        loop_frames{

                            fscanf(fi,"%d \n", &Rec.jumps[i_state][i_state2][i_island][i_frame]);
                        }

                        fclose(fi);

                        number_of_jump_files++;
                    }
#endif
                    loop_frames{

                        host_to_node_int_1(Rec.jumps[i_state][i_state2][i_island][i_frame]);
                    }
                }
            }
        }

        host_to_node_int_1(number_of_jump_files);

        if(Transcript){

            FILE    *f_out = fopen("./Run.trn", "a");

            if(f_out){

                fprintf(f_out,"\n\nrCFD_default_Rec");

                fprintf(f_out,"\n\n   Read %d jump files from %s\n", number_of_jump_files, File_Dict.Jump_filename);

                fclose(f_out);
            }

        }

#if RP_HOST
        if(number_of_jump_files > 0){

            Message("\n\n...rCFD_default_Rec . Read %d jump files from %s", number_of_jump_files, File_Dict.Jump_filename);
        }
#endif

    }
#endif
