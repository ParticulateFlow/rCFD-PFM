#ifndef RCFD_INIT
#define RCFD_INIT

#include "rCFD_types.h"
#include "rCFD_globals.h"
#include "rCFD_macros.h"
#include "rCFD_parallel.h"
#include "rCFD_user.h"

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
        
        Topo_Dict.Face_Dict = (Face_Dict_type*)malloc(Solver_Dict.number_of_layers * sizeof(Face_Dict_type));

        rCFD_default_Cell_Dict_L0();

        rCFD_default_Face_Dict_L0();
#endif   
    }    
    
    /* G1. Topo */
    {
        rCFD_default_Topo();
        
#if RP_NODE     
        int i_layer, i_phase, i_frame, i_cell, i_face;
        
        Topo.Cell = (Cell_type*)malloc(Solver_Dict.number_of_layers * sizeof(Cell_type));
        
        Topo.Face = (Face_type*)malloc(Solver_Dict.number_of_layers * sizeof(Face_type));

        i_layer = 0;
                
        /* T.1. allocate cells for L0 */
        {
            _C.x = (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
            
            loop_cells{
                
                _C.x[i_cell] = (double*)malloc( 3 * sizeof(double));
            }
            
            _C.volume = (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));
            
            _C.average_velocity =    (double**)malloc(Solver_Dict.number_of_phases * sizeof(double*));
            _C.crossing_time =       (double**)malloc(Solver_Dict.number_of_phases * sizeof(double*));
            
            loop_phases{
                
                _C.average_velocity[i_phase] =   (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));            
                _C.crossing_time[i_phase] =      (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));
            }

            _C.hit_by_other_cell =   (short*)malloc(_Cell_Dict.number_of_cells * sizeof(short));

            _C.island_id =           (short*)malloc(_Cell_Dict.number_of_cells * sizeof(short));
            
            _C.weight_after_shift =  (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));
            _C.weight_after_swap =   (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));

            _C.vof = (double***)malloc(Solver_Dict.number_of_frames * sizeof(double**));
            
            loop_frames{
                
                _C.vof[i_frame] = (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
                
                loop_cells{
                    
                    _C.vof[i_frame][i_cell] = (double*)malloc(Solver_Dict.number_of_phases * sizeof(double));
                }
            }

            _C.data =        (double***)malloc(Solver_Dict.number_of_phases * sizeof(double**));
            _C.data_shift =  (double***)malloc(Solver_Dict.number_of_phases * sizeof(double**));
            _C.data_swap =   (double***)malloc(Solver_Dict.number_of_phases * sizeof(double**));
            
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
                
                _C.user = (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
                
                loop_cells{
                    
                    _C.user[i_cell] = (double*)malloc(_Cell_Dict.number_of_user_vars * sizeof(double));
                }
            }
            else{
                
                _C.user = NULL;
            }                        
        
			/* allocation is done in init_layer */
			_C.marked				= NULL;
			_C.parent_cell			= NULL;
			_C.number_of_children	= NULL;
			_C.children				= NULL;						
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
				
				_C.x					= NULL;
				_C.volume				= NULL;

				_C.average_velocity		= NULL;
				_C.crossing_time		= NULL;

				_C.hit_by_other_cell	= NULL;
				_C.island_id			= NULL;

				_C.weight_after_shift	= NULL;
				_C.weight_after_swap	= NULL;

				_C.vof					= NULL;

				_C.data					= NULL;
				_C.data_shift			= NULL;
				_C.data_swap			= NULL;

				_C.drift_exchange		= NULL;

				_C.user					= NULL;
				_C.rec_user				= NULL;
						
				_C.marked				= NULL;
				_C.parent_cell			= NULL;
				_C.number_of_children	= NULL;
				_C.children				= NULL;				
			
				_F.c0 					= NULL;
				_F.c1 					= NULL;
				
				_F.area					= NULL;			
			}
		}

#endif      
    }
    
    /* G1. Cells (first initialization) */
    {
#if RP_NODE     
        int     i_phase, i_frame, i_cell, i_layer;
        
        loop_layers{
            
            if(i_layer == 0){
            
                C.x = (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
                
                loop_cells{
                    
                    C.x[i_cell] = (double*)malloc( 3 * sizeof(double));
                }
                
                C.volume = (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));
                
                C.average_velocity =    (double**)malloc(Solver_Dict.number_of_phases * sizeof(double*));
                C.crossing_time =       (double**)malloc(Solver_Dict.number_of_phases * sizeof(double*));
                
                loop_phases{
                    
                    C.average_velocity[i_phase] =   (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));            
                    C.crossing_time[i_phase] =      (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));
                }

                C.hit_by_other_cell =   (short*)malloc(_Cell_Dict.number_of_cells * sizeof(short));

                C.island_id =           (short*)malloc(_Cell_Dict.number_of_cells * sizeof(short));
                
                C.weight_after_shift =  (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));
                C.weight_after_swap =   (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));

                C.vof = (double***)malloc(Solver_Dict.number_of_frames * sizeof(double**));
                
                loop_frames{
                    
                    C.vof[i_frame] = (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
                    
                    loop_cells{
                        
                        C.vof[i_frame][i_cell] = (double*)malloc(Solver_Dict.number_of_phases * sizeof(double));
                    }
                }

                C.data =        (double***)malloc(Solver_Dict.number_of_phases * sizeof(double**));
                C.data_shift =  (double***)malloc(Solver_Dict.number_of_phases * sizeof(double**));
                C.data_swap =   (double***)malloc(Solver_Dict.number_of_phases * sizeof(double**));
                
                loop_phases{
                    
                    C.data[i_phase] =       (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
                    C.data_shift[i_phase] = (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
                    C.data_swap[i_phase] =  (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
                    
                    loop_cells{
                        
                        C.data[i_phase][i_cell] =       (double*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(double));
                        C.data_shift[i_phase][i_cell] = (double*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(double));
                        C.data_swap[i_phase][i_cell] =  (double*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(double));
                    }
                }       
                
                if(Solver_Dict.data_drifting_on){
                    
                    C.drift_exchange = (double*)malloc(_Cell_Dict.number_of_cells * sizeof(double));
                }
                else{
                    
                    C.drift_exchange = NULL;
                }       
                
                if(_Cell_Dict.number_of_user_vars > 0){
                    
                    C.user = (double**)malloc(_Cell_Dict.number_of_cells * sizeof(double*));
                    
                    loop_cells{
                        
                        C.user[i_cell] = (double*)malloc(_Cell_Dict.number_of_user_vars * sizeof(double));
                    }
                }
                else{
                    
                    C.user = NULL;
                }
                        
                rCFD_default_Cell(i_layer);
            }
        }
#endif  
    }   

    /* G2. Faces */
    {
#if RP_NODE
        int i_face, i_layer;
        
        loop_layers{
            
            if(i_layer == 0){
                
                F.c0 = (int*)malloc(_Face_Dict.number_of_faces * sizeof(int));
        
                F.c1 = (int*)malloc(_Face_Dict.number_of_faces * sizeof(int));
                
                F.area = (double**)malloc(_Face_Dict.number_of_faces * sizeof(double*));
                
                loop_faces{
                    
                    F.area[i_face] = (double*)malloc( 3 * sizeof(double));
                }

                rCFD_default_Face(i_layer);
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
        
        rCFD_user_post(i_layer);   /* post-process initialization */
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
		int		i_layer, i_cell, i_face, i_cell_face, i_face2, i_cell_face2, i_child, i_dim, i_tmp;
		
		int		cell_index_in_upper_layer, number_of_cells_in_upper_layer, upper_cell, c0, c1, c1_c0, c1_c1, min_dist_parent_cell;
		
		double	dist, min_dist;
		
		short	debug_this_code = 1;
		
		int		**tmp_faces_of_cell = NULL;	/* [i_cell][i_cell_face] */
		
		i_layer = 0;
		
		while(	((i_layer + 1) < Solver_Dict.number_of_layers) && 
		
				(_Cell_Dict.number_of_cells > Solver_Dict.min_number_of_cells_per_layer))
		{

			/* L.1 create tmp_faces_of_cell */
			{

				tmp_faces_of_cell = (int**)malloc(_Cell_Dict.number_of_cells * sizeof(int*));
				
				loop_cells{
					
					tmp_faces_of_cell[i_cell] = (int*)malloc(Solver_Dict.max_number_of_faces_per_cell * sizeof(int));
					
					for(i_cell_face = 0; i_cell_face < Solver_Dict.max_number_of_faces_per_cell; i_cell_face++){
						
						tmp_faces_of_cell[i_cell][i_cell_face] = -1;
					}
				}
				
				loop_faces{
					
					c0 = _F.c0[i_face];
					c1 = _F.c1[i_face];
					
					i_cell_face = 0;
					
					while((tmp_faces_of_cell[c0][i_cell_face] > 0) && (i_cell_face < Solver_Dict.max_number_of_faces_per_cell)){
						
						i_cell_face++;
					}
					
					if(i_cell_face < Solver_Dict.max_number_of_faces_per_cell){
						
						tmp_faces_of_cell[c0][i_cell_face] = i_face;
					}
					else{
						
						Message("\nWARNING: init layers, L1: myid %d i_cell_face %d\n", myid, i_cell_face);
					}

					i_cell_face = 0;
					
					while((tmp_faces_of_cell[c1][i_cell_face] > 0) && (i_cell_face < Solver_Dict.max_number_of_faces_per_cell)){
						
						i_cell_face++;
					}
					
					if(i_cell_face < Solver_Dict.max_number_of_faces_per_cell){
						
						tmp_faces_of_cell[c1][i_cell_face] = i_face;
					}
					else{
						
						Message("\nWARNING: init layers, L1: myid %d i_cell_face %d\n", myid, i_cell_face);
					}
				}
			}				
					
			/* L.2 create cell clusters, set number of upper layer cells */
			{				
				_C.marked = (short*)malloc(_Cell_Dict.number_of_cells * sizeof(short));
				_C.parent_cell = (int*)malloc(_Cell_Dict.number_of_cells * sizeof(int));
				
				enum{
					not_yet_visited,
					center_of_cluster,
					member_of_cluster,
					neighbor_of_cluster
				};
					
				loop_cells{
					
					_C.marked[i_cell] = not_yet_visited;
					_C.parent_cell[i_cell] = -1;
				}
				
				cell_index_in_upper_layer = 0;
				
				loop_cells{
					
					if(_C.marked[i_cell] == not_yet_visited){
						
						_C.marked[i_cell] = center_of_cluster;
						
						_C.parent_cell[i_cell] = cell_index_in_upper_layer;
						
						/* also mark neighbors */
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

								/* also mark neighbors of neighbor */
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

				if(debug_this_code){
										
					i_tmp = 0;
				
					loop_cells{
						
						if(_C.marked[i_cell] == neighbor_of_cluster){
							
							i_tmp++;
						}
					}
								
					Message("\nDEBUG myid %d Topo_Dict.Cell_Dict[%d].number_of_cells %d; unused neighbors %d\n", 
					
						myid, upper_layer, Topo_Dict.Cell_Dict[upper_layer].number_of_cells, i_tmp);
				}

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
							
							if(_C.marked[c1] != neighbor_of_cluster){
								
								dist = (_C.x[c0][0] - _C.x[c1][0])*(_C.x[c0][0] - _C.x[c1][0]) +
									   (_C.x[c0][1] - _C.x[c1][1])*(_C.x[c0][1] - _C.x[c1][1]) +
									   (_C.x[c0][2] - _C.x[c1][2])*(_C.x[c0][2] - _C.x[c1][2]);
								
								if(dist < min_dist){
									
									min_dist = dist;
									
									min_dist_parent_cell = _C.parent_cell[c1];
								}
							}
							
							i_cell_face++;
						}
						
						_C.marked[i_cell] = member_of_cluster;
						
						_C.parent_cell[i_cell] = min_dist_parent_cell;						
					}
				}
				
				if(debug_this_code){
										
					i_tmp = 0;
				
					loop_cells{
						
						if(_C.marked[i_cell] == neighbor_of_cluster){
							
							i_tmp++;
						}
					}
								
					Message("\nDEBUG myid %d unused neighbors %d\n", 
					
						myid, i_tmp);
				}
				
			}
			
			/* L.2 allocate and initialize upper layer grid vars */
			{
				Topo.Cell[upper_layer].x = 					(double**)malloc(number_of_cells_in_upper_layer * sizeof(double*));
				
				Topo.Cell[upper_layer].volume = 			(double*)malloc(number_of_cells_in_upper_layer * sizeof(double));

				Topo.Cell[upper_layer].number_of_children = (short*)malloc(number_of_cells_in_upper_layer * sizeof(short));

				Topo.Cell[upper_layer].children = 			(int**)malloc(number_of_cells_in_upper_layer * sizeof(int*));
				
				loop_cells_of_upper_layer{
				
					Topo.Cell[upper_layer].x[i_cell] = (double*)malloc(3 * sizeof(double));
					
					loop_dim{
						
						Topo.Cell[upper_layer].x[i_cell][i_dim] = 0.0;
					}

					Topo.Cell[upper_layer].volume[i_cell] = 0.0;
					
					Topo.Cell[upper_layer].number_of_children[i_cell] = 0;
				}
			}
			
			/* L.3 fill upper layer grid vars */
			{
				loop_cells{
					
					if(_C.marked[i_cell] == 1){
						
						upper_cell = _C.parent_cell[i_cell];

						loop_dim{
							
							Topo.Cell[upper_layer].x[upper_cell][i_dim] += _C.x[i_cell][i_dim];
						}
					
						Topo.Cell[upper_layer].volume[upper_cell] += _C.volume[i_cell];
					
						Topo.Cell[upper_layer].number_of_children[upper_cell] += 1;
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
			}
			
			/* L.X fill upper layer list of children */
			{
				/* allocate upper layer list of children */
				loop_cells_of_upper_layer{
					
					Topo.Cell[upper_layer].children[i_cell] = 
					
						(int*)malloc(Topo.Cell[upper_layer].number_of_children[i_cell] * sizeof(int));
					
					for(i_child = 0; i_child < Topo.Cell[upper_layer].number_of_children[i_cell]; i_child++){
						
						Topo.Cell[upper_layer].children[i_cell][i_child] = -1;
					}					
				}

				/* fill upper layer list of children */
				loop_cells{
					
					if(_C.parent_cell[i_cell] > -1){
						
						upper_cell = _C.parent_cell[i_cell];
						
						{
							i_child = 0;
							
							while(	(Topo.Cell[upper_layer].children[upper_cell][i_child] < 0) && 
							
									(i_child < Solver_Dict.max_number_of_children)){
								
								i_child++;
							}
							
							if(i_child < Solver_Dict.max_number_of_children){
								
								Topo.Cell[upper_layer].children[upper_cell][i_child] = i_cell;
							}
							else{
								
								Message("ERROR: init_layers, LX: myid %d, i_layer %d i_child %d\n", myid, i_layer, i_child);
								
								return;
							}
						}
					}
				}
			}
			
			/* L.Y free local vars */
			{
				free(tmp_faces_of_cell);
			}
			
			i_layer++;
		}
#endif		
	}
}

#endif