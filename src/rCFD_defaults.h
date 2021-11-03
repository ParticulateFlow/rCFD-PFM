#ifndef RCFD_DEFAULTS
#define RCFD_DEFAULTS

#include "udf.h"
#include <sys/stat.h>
#include "rCFD_macros.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

/* (C) 	2021 
	Stefan Pirker
	Particulate Flow Modelling
	Johannes Kepler University, Linz, Austria
	www.particulate-flow.at
*/	

	#define		Global_Version_Year		2021
	#define		Global_Version_Month	11
	
	#define		Global_Verbal			1
	
	/* A. default Dicts */
	
	void rCFD_default_Solver_Dict(Solver_Dict_type *Solver_Dict)
	{	
		Solver_Dict->version_year = 	Global_Version_Year;
		Solver_Dict->version_month = 	Global_Version_Month;
		
		Solver_Dict->verbal = 			Global_Verbal;

		Solver_Dict->number_of_frames = 	1;			
		Solver_Dict->number_of_states = 	1;			
		Solver_Dict->number_of_phases = 	1;			
		Solver_Dict->number_of_islands = 	1;			
		Solver_Dict->number_of_runs = 		1;	
		
		Solver_Dict->recurrence_process_on = 	1;
		Solver_Dict->data_convection_on = 		1;
		Solver_Dict->face_diffusion_on = 		1;
		Solver_Dict->data_binarization_on = 	0;
		Solver_Dict->data_drifting_on = 		0;
		Solver_Dict->balance_correction_on = 	1;
		Solver_Dict->on_the_fly_post_on = 		1;
	
		Solver_Dict->analyse_CFD_count = 0;
		
		Solver_Dict->max_number_of_cells_per_time_step = 10;
		
		Solver_Dict->time_steps_per_monitoring_interval	=	10;
		Solver_Dict->start_time_for_monitoring = 			0.0;
		
		Solver_Dict->C2C_loading_reduction = 	1;
		
		Solver_Dict->global_time_step = 1.0;
			
		Solver_Dict->max_fill_loops = 	10;

		Solver_Dict->face_swap_max_per_loop = 	(1./8.);		
		Solver_Dict->face_swap_min = 			(1./64.);
		Solver_Dict->face_swap_max_loops = 		4;	
		
		Solver_Dict->number_of_drift_loops = 	1;
		
		Solver_Dict->balance_correction_update = 1;
		Solver_Dict->control_conc_sum_on = 1;
		
		if(Transcript){
							
			FILE	*f_out = fopen("./Run.trn", "w");
			
			if(f_out){
				
				fprintf(f_out,"rCFD_default_Solver_Dict");
			
				fprintf(f_out,"\n\n   Version: %4d.%02d", Solver_Dict->version_year, Solver_Dict->version_month);
				
				time_t		current_time = time(NULL);
				
				char 		*c_time_string = ctime(&current_time); 

				fprintf(f_out,"\n\n   Run start @ %s", c_time_string);

				fclose(f_out);
			}	
		}
	}	

	void rCFD_default_File_Dict(Solver_Dict_type *Solver_Dict, File_Dict_type *File_Dict)
	{	

		/* Node-0 generates folder structure */
		if(myid == 0){
			
			if(	(mkdir("./DATA",0777) != 0) ||
				(mkdir("./DATA/TMP",0777) != 0) ||
				(mkdir("./REC",0777) != 0) ||
				(mkdir("./DICT",0777) != 0) ||
				(mkdir("./POST",0777) != 0) ){
				
				/*Message0("\nERROR rCFD_default_File_Dict");*/
			}
		}
		
		File_Dict->tracer_start_position_filename	= 	"./DATA/TMP/Tracer_start_pos.inj";
		
		File_Dict->C2C_filename = 						"./DATA/C2C/C2C";
		
		File_Dict->Norm_filename = 						"./DATA/TMP/Norm";
		
		File_Dict->Jump_filename = 						"./REC/Jump";

		File_Dict->Matrix_filename = 					"./REC/Matrix";

		File_Dict->Balance_filename = 					"./POST/Balance_monitor.out";		
		
		File_Dict->Dict_filename = 						"./DICT/Dict_protocol";
	}
	
	void rCFD_default_Phase_Dict(Solver_Dict_type *Solver_Dict, Phase_Dict_type *Phase_Dict)
	{
		int i_phase;
		
		for(i_phase = 0; i_phase < Solver_Dict->number_of_phases; i_phase++){
											
			Phase_Dict[i_phase].number_of_data = 1;

			Phase_Dict[i_phase].time_step = 1.0;

			Phase_Dict[i_phase].shift_probability = 1.0;
			
			Phase_Dict[i_phase].density = 1.0;

			Phase_Dict[i_phase].heat_capacity = 1.0;

			Phase_Dict[i_phase].reference_temperature = 300.0;
			
			Phase_Dict[i_phase].number_of_user_vars = 0;
			
			Phase_Dict[i_phase].user = NULL;
		}			
	}
	
	void rCFD_default_Tracer_Dict(Solver_Dict_type *Solver_Dict, Tracer_Dict_type *Tracer_Dict)
	{
		
		Tracer_Dict->number_of_Tracers_per_cell	= 1;
		
		Tracer_Dict->region_of_interest_exists = 0;
		
		Tracer_Dict->ROI_x_min = 0.0;
		Tracer_Dict->ROI_x_max = 0.0;
		Tracer_Dict->ROI_y_min = 0.0;
		Tracer_Dict->ROI_y_max = 0.0;
		Tracer_Dict->ROI_z_min = 0.0;
		Tracer_Dict->ROI_z_max = 0.0;
		
		Tracer_Dict->coarse_graining = 1;
		
		int	i_phase;
		
		loop_phases_ptr{
			
			Tracer_Dict->random_walk[i_phase] = 0;
		}
		
		Tracer_Dict->C2C_format = c0_n0_c1_n1_w0;
	}
	
	void rCFD_default_Norm_Dict(Solver_Dict_type *Solver_Dict, Norm_Dict_type *Norm_Dict)
	{
		
		Norm_Dict->format = standard;
		
		Norm_Dict->coarse_graining = 1;
	}
	
	void rCFD_default_Rec_Dict(Solver_Dict_type *Solver_Dict, Rec_Dict_type *Rec_Dict)
	{
		
		Rec_Dict->method = quarter_jumps;
		
		Rec_Dict->min_seq_length = (int)Solver_Dict->number_of_frames/25;
		
		if(Rec_Dict->min_seq_length < 1){
			
			Rec_Dict->min_seq_length = 1;
		}
		
		Rec_Dict->max_seq_length = (int)Solver_Dict->number_of_frames/4;
		
		if(Rec_Dict->max_seq_length < 1){
			
			Rec_Dict->max_seq_length = 1;
		}		
	}
	
	void rCFD_default_Data_Dict(Solver_Dict_type *Solver_Dict, Phase_Dict_type *Phase_Dict, Data_Dict_type **Data_Dict)
	{
		int	i_phase, i_data;

		loop_phases_ptr{

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
	}
	
	void rCFD_default_Balance_Dict(Solver_Dict_type *Solver_Dict, Phase_Dict_type *Phase_Dict, Balance_Dict_type **Balance_Dict)
	{
		int	i_phase, i_data;

		loop_phases_ptr{

			loop_data{
			
				Balance_Dict[i_phase][i_data].type = global_balancing;
				
				Balance_Dict[i_phase][i_data].write_balance_to_file = 0;
				
				Balance_Dict[i_phase][i_data].write_balance_to_file_interval = 1;
			}
		}
	}	
	
	void rCFD_default_Cell_Dict(Solver_Dict_type *Solver_Dict, Cell_Dict_type *Cell_Dict)
	{
		Domain 	*d=Get_Domain(1);
		Thread	*t;
		cell_t 	c;

		int		number_of_local_cells = 0, number_of_local_internal_cells = 0;
		
		int 	only_one_fluid_cell_thread = 0;
		
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
			
		Cell_Dict->number_of_cells = number_of_local_cells;	

		Cell_Dict->number_of_int_cells = number_of_local_internal_cells;

		Cell_Dict->number_of_ext_cells = number_of_local_cells - number_of_local_internal_cells;
		
		Cell_Dict->number_of_user_vars = 0;
	}
		
	void rCFD_default_Face_Dict(Solver_Dict_type *Solver_Dict, Face_Dict_type *Face_Dict)
	{
		Domain 	*d=Get_Domain(1);
		Thread	*t;
		face_t 	f;

		int		number_of_local_faces = 0;
		
		thread_loop_f(t,d){if(THREAD_TYPE(t)==THREAD_F_INTERIOR){begin_f_loop(f,t){
				
				number_of_local_faces++;

		}end_f_loop(f,t)}}
		
		Face_Dict->number_of_faces = number_of_local_faces;
	}

	/* B. default vars settings */
	
	void rCFD_default_Solver(Solver_type *Solver)
	{	
		Solver->current_state = 0;
		
		Solver->global_run_counter = 	0;			
	}	

	void rCFD_default_Cell(Solver_Dict_type *Solver_Dict, Phase_Dict_type *Phase_Dict, Cell_Dict_type *Cell_Dict, Cell_type *C)
	{
		int	i_phase, i_cell, i_frame, i_user, i_data, i_dim;
		
		Domain 	*d=Get_Domain(1);
		Thread	*t;
		double	x[3];
		
		thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){begin_c_loop_all(i_cell,t){
				
				C_CENTROID(x, i_cell, t);
				
				loop_dim{
					
					C->x[i_cell][i_dim] = x[i_dim];
				}
				
				C->volume[i_cell] = C_VOLUME(i_cell, t);
				
		}end_c_loop_all(i_cell,t)}}
	
		loop_phases_ptr{
			
			loop_cells_ptr{
				
				C->average_velocity[i_phase][i_cell] = 	0.0;				
				C->crossing_time[i_phase][i_cell] = 	0.0;
			}
		}
		
		loop_cells_ptr{
			
			C->hit_by_other_cell[i_cell] = 0;
			C->island_id[i_cell] = 0;
			
			C->weight_after_shift[i_cell] = 0.0;
			C->weight_after_swap[i_cell] = 	0.0;
		}
		
		if(C->vof != NULL){
				
			loop_frames_ptr{
				
				loop_cells_ptr{
					
					loop_phases_ptr{
					
						if(i_phase == 0){
							
							C->vof[i_frame][i_cell][i_phase] = 1.0;
						}
						else{
						
							C->vof[i_frame][i_cell][i_phase] = 0.0;
						}
					}
				}
			}		
		}
		
		loop_phases_ptr{
			
			loop_cells_ptr{
				
				loop_data{
				
					C->data[i_phase][i_cell][i_data] = 		0.0;
					C->data_shift[i_phase][i_cell][i_data] = 0.0;
					C->data_swap[i_phase][i_cell][i_data] = 	0.0;
				}
			}
		}		

		if(C->drift_exchange != NULL){
			
			loop_cells_ptr{
				
				C->drift_exchange[i_cell] = 0.0;
			}
		}
		
		if(C->user != NULL){
			
			loop_cells_ptr{
				
				for(i_user = 0; i_user < Cell_Dict->number_of_user_vars; i_user++){
					
					C->user[i_cell][i_user] = 0.0;
				}
			}
		}
	}
	
	void rCFD_default_Face(Solver_Dict_type *Solver_Dict, Cell_Dict_type *Cell_Dict, Face_Dict_type *Face_Dict, Face_type *F)
	{
		Domain 	*d=Get_Domain(1);
		Thread	*t;
		int		f, i_face, c0, c1, i_dim;
		
		double	A[3];
		
		i_face = 0;
		
		thread_loop_f(t,d){if(THREAD_TYPE(t)==THREAD_F_INTERIOR){begin_f_loop(f,t){
				
				c0 = F_C0(f, t);

				c1 = F_C1(f, t);
					
				F->c0[i_face] = c0;
				
				F->c1[i_face] = c1;
				
				F_AREA(A, f, t);
				
				loop_dim{
					
					F->area[i_face][i_dim] = A[i_dim];
				}
				
				i_face++;

		}end_f_loop(f,t)}}
	}
	
	void rCFD_default_Tracer(Solver_Dict_type *Solver_Dict, Tracer_type *Tracer)
	{
		Tracer->allocated = 	0;
		Tracer->initialized = 	0;
		Tracer->ready2write = 	0;
		Tracer->monitoring_started  = 0;		
		
		Tracer->frame_counter = 	0;
		
		int 	i_phase;
		
		loop_phases_ptr{
			
			Tracer->monitor_counter[i_phase] = 	0;
			Tracer->number_of_shifts[i_phase] = 0;
			
			Tracer->shifts[i_phase] = NULL;
		}		
	}
	
	void rCFD_default_Norms(Solver_Dict_type *Solver_Dict, Norm_type *Norms)
	{
		Norms->frame_counter = 0;
	}	
	
	void rCFD_default_Rec(Solver_Dict_type * Solver_Dict, File_Dict_type *File_Dict, Rec_Dict_type *Rec_Dict, Rec_type *Rec)
	{
		int	i_state, i_state2, i_island, i_frame;
		
		int number_of_jump_files = 0;
		
		loop_islands_ptr{
			
			Rec->global_frame[i_island] = 0;
		}
		
		Rec->frame_in_sequence = 0;
		
		Rec->sequence_length = Rec_Dict->min_seq_length;
		
		loop_states_ptr{
			
			loop_states2_ptr{
				
				loop_islands_ptr{
					
					loop_frames_ptr{
						
						Rec->jumps[i_state][i_state2][i_island][i_frame] = -1;
					}

#if RP_HOST 
					char	filename[40];
					FILE 	*fi;

					sprintf(filename,"%s_%d_%d_%d", File_Dict->Jump_filename, i_state, i_state2, i_island);
					
					fi=fopen(filename,"r");
					
					if(fi!=NULL){
						
						fopen(filename,"r");
						
						loop_frames_ptr{ 
						
							fscanf(fi,"%d \n", &Rec->jumps[i_state][i_state2][i_island][i_frame]);
						}
						
						fclose(fi);
						
						number_of_jump_files++;
					}
#endif				
					loop_frames_ptr{ 
					
						host_to_node_int_1(Rec->jumps[i_state][i_state2][i_island][i_frame]);
					}
				}
			}
		}
		
		if(Transcript){

			FILE	*f_out = fopen("./Run.trn", "a");
			
			if(f_out){
				
				fprintf(f_out,"\n\nrCFD_default_Rec");
			
				fprintf(f_out,"\n\n   Read %d jump files from %s", number_of_jump_files, File_Dict->Jump_filename);

				fclose(f_out);
			}
		
		}
		
		Message0("\n\n...rCFD_default_Rec -> Read %d jump files from %s", 
			
			number_of_jump_files, File_Dict->Jump_filename);
	}
#endif	