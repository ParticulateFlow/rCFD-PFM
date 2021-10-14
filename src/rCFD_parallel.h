#ifndef RCFD_PARALLEL
#define RCFD_PARALLEL

#include "rCFD_types.h"
#include "rCFD_macros.h"

/* (C) 	2021 
	Stefan Pirker
	Particulate Flow Modelling
	Johannes Kepler University, Linz, Austria
	www.particulate-flow.at
*/	

#if 1 /* macros */

#define valid_parallel_face 	MPI_faces.principal_face[i_face]

#endif

#if 1 /* typedef */
#if RP_NODE

	typedef struct MPI_cells_struct
	{
		int		number_of_ext_cells;

		int		*number_of_ext_cells_per_node;

		int		*number_of_host_cells_per_node;
							
		int		*host_of_cell;				/* all nodes, list of host of all cells */

		int		*host_of_ext_cells;			/* node-0, list of all ext. cells */
		
		int		*hosting_cell_index;		/* all nodes, list of cells, which host ext. cells of other nodes */

		int		*host2ext_index;	
		
		double	*data;
		
	} MPI_cells_type;
	
	typedef struct MPI_faces_struct
	{
		int		*principal_face;
		
	} MPI_faces_type;
	
#endif
#endif
	
#if 1 /* global vars */
#if RP_NODE

	static	MPI_cells_type		MPI_cells;
	
	static	MPI_faces_type		MPI_faces;
	
#endif
#endif	

#if 1 /* defaults */
#if RP_NODE

	void rCFD_default_MPI_Faces(Solver_Dict_type *Solver_Dict, MPI_faces_type *MPI_faces)
	{
		Domain 	*d=Get_Domain(1);
		Thread	*t;
		int		i_face;
		
		thread_loop_f(t,d){if(THREAD_TYPE(t)==THREAD_F_INTERIOR){begin_f_loop(i_face,t){
				
			if(PRINCIPAL_FACE_P(i_face,t)){
				
				MPI_faces->principal_face[i_face] = 1;
			}
			else{
				
				MPI_faces->principal_face[i_face] = 0;
			}

		}end_f_loop(i_face,t)}}
	}
		
	void rCFD_default_MPI_Cells(Solver_Dict_type *Solver_Dict, MPI_cells_type *MPI_cells)
	{
		int	i_cell;
		
		Domain 	*d=Get_Domain(1);
		Thread	*t;
		
		thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){begin_c_loop_all(i_cell,t){
				
				MPI_cells->host_of_cell[i_cell] = C_PART(i_cell,t);

		}end_c_loop_all(i_cell,t)}}		
	}
		
#endif
#endif

#if 1 /* init */
#if RP_NODE

	void init_parallel_grid(Solver_Dict_type *Solver_Dict, Cell_Dict_type *Cell_Dict, Face_Dict_type *Face_Dict)
	{
		int 		i_node, i_cell, i_list, i_host;
		int 		size;
		int 		*list_of_int = NULL;
		double 		*list_of_real = NULL;
		
		/* P1: allocate MPI_cells */
		{
			MPI_cells.host_of_cell = (int*)malloc(Cell_Dict->number_of_cells * sizeof(int));

			MPI_cells.hosting_cell_index = NULL;

			if(myid > 0){
				
				PRF_CSEND_INT(node_zero, &Cell_Dict->number_of_ext_cells, 1, myid);
				
				MPI_cells.number_of_ext_cells = -1;
				
				MPI_cells.number_of_ext_cells_per_node = NULL;
				
				MPI_cells.number_of_host_cells_per_node = NULL;
				
				MPI_cells.host_of_ext_cells = NULL;
						
				MPI_cells.host2ext_index = NULL;	
				
				MPI_cells.data = NULL;			
			}
			
			if(myid == 0){
				
				MPI_cells.number_of_ext_cells_per_node = (int*)malloc((node_last+1) * sizeof(int));
				
				MPI_cells.number_of_host_cells_per_node = (int*)malloc((node_last+1) * sizeof(int));	
				
				MPI_cells.number_of_ext_cells_per_node[0] = Cell_Dict->number_of_ext_cells;
				
				/* get number_of_ext_cells (sum up all nodes) */
				MPI_cells.number_of_ext_cells = Cell_Dict->number_of_ext_cells;
				
				for(i_node = 1; i_node < (node_last+1); i_node++){
					
					PRF_CRECV_INT(i_node, &MPI_cells.number_of_ext_cells_per_node[i_node], 1, i_node);
					
					MPI_cells.number_of_ext_cells += MPI_cells.number_of_ext_cells_per_node[i_node];
				}
				
				MPI_cells.host_of_ext_cells = (int*)malloc(MPI_cells.number_of_ext_cells * sizeof(int));

				MPI_cells.host2ext_index = (int*)malloc(MPI_cells.number_of_ext_cells * sizeof(int));	
				
				MPI_cells.data = (double*)malloc(MPI_cells.number_of_ext_cells * sizeof(double));
			}
			
			rCFD_default_MPI_Cells(Solver_Dict, &MPI_cells);	
		}
			
		/* P2. MPI_cells communication */
		{			
			double *tmp_x;
			
			/* P2.1 fill MPI_cells.host_of_ext_cells, tmp_x */
			{
				if(myid > 0){
					
					list_of_int = (int*)malloc( Cell_Dict->number_of_ext_cells * sizeof(int));
					
					list_of_real = (double*)malloc( 3 * Cell_Dict->number_of_ext_cells * sizeof(double));
					
					i_list = 0;
					
					loop_ext_cells_ptr{
						
						list_of_int[i_list] = MPI_cells.host_of_cell[i_cell];
						
						list_of_real[(3 * i_list)] = 		C.x[i_cell][0];
						list_of_real[(3 * i_list + 1)] = 	C.x[i_cell][1];
						list_of_real[(3 * i_list + 2)] = 	C.x[i_cell][2];
						
						i_list++;
					}
					
					PRF_CSEND_INT(node_zero, list_of_int, Cell_Dict->number_of_ext_cells, myid);

					PRF_CSEND_REAL(node_zero, list_of_real, (3 * Cell_Dict->number_of_ext_cells), myid);
					
					free(list_of_int);
					
					free(list_of_real);
				}
				
				if(myid == 0){
								
					tmp_x = (double*)malloc(3 * MPI_cells.number_of_ext_cells * sizeof(double));

					i_list = 0;
					
					loop_ext_cells_ptr{
						
						MPI_cells.host_of_ext_cells[i_list] = MPI_cells.host_of_cell[i_cell];
						
						tmp_x[(3 * i_list)] = 		C.x[i_cell][0];
						tmp_x[(3 * i_list + 1)] = 	C.x[i_cell][1];
						tmp_x[(3 * i_list + 2)] = 	C.x[i_cell][2];
						
						i_list++;
					}

					/* get ext cell info from Nodes */
					for(i_node = 1; i_node < (node_last+1); i_node++){
						
						list_of_int = (int*)malloc(MPI_cells.number_of_ext_cells_per_node[i_node] * sizeof(int));
					
						list_of_real = (double*)malloc( 3 * MPI_cells.number_of_ext_cells_per_node[i_node] * sizeof(double));
						
						PRF_CRECV_INT(i_node, list_of_int, MPI_cells.number_of_ext_cells_per_node[i_node], i_node);

						PRF_CRECV_REAL(i_node, list_of_real, (3 * MPI_cells.number_of_ext_cells_per_node[i_node]), i_node);
						
						for(i_cell = 0; i_cell < MPI_cells.number_of_ext_cells_per_node[i_node]; i_cell++){
							
							MPI_cells.host_of_ext_cells[i_list] = list_of_int[i_cell];
							
							tmp_x[(3 * i_list)] = 		list_of_real[(3 * i_cell)];
							tmp_x[(3 * i_list + 1)] = 	list_of_real[(3 * i_cell + 1)];
							tmp_x[(3 * i_list + 2)] = 	list_of_real[(3 * i_cell + 2)];
										
							i_list++;
						}
						
						free(list_of_int);
						
						free(list_of_real);
					}
				}
			}
			
			/* P2.2 analyze MPI_cells.host_of_ext_cells, set MPI_cells.host2ext_index */
			{
				if(myid == 0){
				
					for(i_host = 1; i_host < (node_last+1); i_host++){
						
						MPI_cells.number_of_host_cells_per_node[i_host] = 0;
					}
					
					for(i_cell = 0; i_cell < MPI_cells.number_of_ext_cells; i_cell++){
						
						i_host = MPI_cells.host_of_ext_cells[i_cell];
						
						MPI_cells.number_of_host_cells_per_node[i_host]++;
					}
					
					int	host_start_index[(node_last+1)], host_index[(node_last+1)];
						
					host_start_index[0] = 0; host_index[0] = 0;
					
					for(i_host = 1; i_host < (node_last+1); i_host++){
						
						host_start_index[i_host] = host_start_index[(i_host-1)] + MPI_cells.number_of_host_cells_per_node[(i_host-1)];
						
						host_index[i_host] = 0;
					}	
					
					for(i_cell = 0; i_cell < MPI_cells.number_of_ext_cells; i_cell++){
						
						i_host = MPI_cells.host_of_ext_cells[i_cell];
						
						MPI_cells.host2ext_index[(host_start_index[i_host] + host_index[i_host])] = i_cell;
						
						host_index[i_host]++;
					}
				}
			}
			
			/* P2.3 set MPI_cells.hosting_cell_index */
			{
				if(myid == 0){
					
					list_of_real = (double*)malloc( 3 * MPI_cells.number_of_host_cells_per_node[0] * sizeof(double));
					
					i_list = 0;
					
					for(i_cell = 0; i_cell < MPI_cells.number_of_ext_cells; i_cell++){
						
						i_host = MPI_cells.host_of_ext_cells[i_cell];
						
						if(i_host == 0){
							
							list_of_real[(3 * i_list)] = 		tmp_x[(3 * i_cell)]; 
							list_of_real[(3 * i_list + 1)] = 	tmp_x[(3 * i_cell + 1)]; 
							list_of_real[(3 * i_list + 2)] = 	tmp_x[(3 * i_cell + 2)]; 
							
							i_list++;
						}
					}
					
					MPI_cells.hosting_cell_index = (int*)malloc(MPI_cells.number_of_host_cells_per_node[0] * sizeof(int));
					
					for(i_list = 0; i_list < MPI_cells.number_of_host_cells_per_node[0]; i_list++){
						
						loop_cells_ptr{
							
							if( (C.x[i_cell][0] == list_of_real[(3 * i_list + 0)]) && 
								(C.x[i_cell][1] == list_of_real[(3 * i_list + 1)]) &&
								(C.x[i_cell][2] == list_of_real[(3 * i_list + 2)])){
								   
								MPI_cells.hosting_cell_index[i_list] = i_cell;
							}
						}
					}				
					
					free(list_of_real);
					
					/* send tmp_x to Nodes */
					for(i_node = 1; i_node < (node_last+1); i_node++){
						
						PRF_CSEND_INT(i_node, &MPI_cells.number_of_host_cells_per_node[i_node], 1, node_zero);
						
						list_of_real = (double*)malloc( 3 * MPI_cells.number_of_host_cells_per_node[i_node] * sizeof(double));
						
						i_list = 0;
						
						for(i_cell = 0; i_cell < MPI_cells.number_of_ext_cells; i_cell++){
							
							i_host = MPI_cells.host_of_ext_cells[i_cell];
							
							if(i_host == i_node){
								
								list_of_real[(3 * i_list)] = 		tmp_x[(3 * i_cell)]; 
								list_of_real[(3 * i_list + 1)] = 	tmp_x[(3 * i_cell + 1)]; 
								list_of_real[(3 * i_list + 2)] = 	tmp_x[(3 * i_cell + 2)]; 
								
								i_list++;
							}
						}
						
						PRF_CSEND_REAL(i_node, list_of_real, (3 * MPI_cells.number_of_host_cells_per_node[i_node]), node_zero);
						
						free(list_of_real);
					}
				
					free(tmp_x);
				}		

				if(myid > 0){
					
					PRF_CRECV_INT(node_zero, &size, 1, node_zero);
					
					MPI_cells.hosting_cell_index = (int*)malloc(size * sizeof(int));
					
					list_of_real = (double*)malloc(3 * size * sizeof(double));
					
					PRF_CRECV_REAL(node_zero, list_of_real, 3 * size, node_zero);
					
					for(i_list = 0; i_list < size; i_list++){
						
						loop_cells_ptr{
							
							if( (C.x[i_cell][0] == list_of_real[(3 * i_list + 0)]) && 
								(C.x[i_cell][1] == list_of_real[(3 * i_list + 1)]) &&
								(C.x[i_cell][2] == list_of_real[(3 * i_list + 2)])){
								   
								MPI_cells.hosting_cell_index[i_list] = i_cell;
							}
						}
					}					
					
					free(list_of_real);
				}
			}
		}	

		/* P3. allocate and set MPI_faces */
		{
			MPI_faces.principal_face = (int*)malloc(Face_Dict->number_of_faces * sizeof(int));
			
			rCFD_default_MPI_Faces(Solver_Dict, &MPI_faces);
		}		
	}
#endif
#endif

#if 1 /* corona cell communication */
#if RP_NODE

	void sum_up_parallel_corona_cells(Solver_Dict_type *Solver_Dict, Cell_Dict_type *Cell_Dict, double *cell_data)
	{
		int		i_cell, i_list, i_node;
		
		int		size;
		
		double	*list_of_real;
				
		/* MPI.1 Node-0 fills MPI_cells.data */
		{
			if(myid > 0){
			
				list_of_real = (double*)malloc(Cell_Dict->number_of_ext_cells * sizeof(double));
				
				i_list = 0;
				
				loop_ext_cells_ptr{
					
					list_of_real[i_list] = cell_data[i_cell];
					
					i_list++;
				}
				
				PRF_CSEND_INT(node_zero, &Cell_Dict->number_of_ext_cells, 1, myid);
				
				PRF_CSEND_REAL(node_zero, list_of_real, Cell_Dict->number_of_ext_cells, myid);
				
				free(list_of_real);
			}
			
			if(myid == 0){
				
				i_list = 0;
				
				loop_ext_cells_ptr{
					
					MPI_cells.data[i_list] = cell_data[i_cell];
					
					i_list++;
				}
				
				for(i_node = 1; i_node < (node_last+1); i_node++){
					
					PRF_CRECV_INT(i_node, &size, 1, i_node);
					
					list_of_real = (double*)malloc(size * sizeof(double));
					
					PRF_CRECV_REAL(i_node, list_of_real, size, i_node);
					
					for(i_cell = 0; i_cell < size; i_cell++){
						
						MPI_cells.data[i_list] = list_of_real[i_cell];
						
						i_list++;
					}
					
					free(list_of_real);
				}
			}
		}
		
		/* MPI.2 Node-0 sends MPI_cells.data to host-Nodes */
		{
			if(myid == 0){
				
				i_list = MPI_cells.number_of_host_cells_per_node[0];
				
				for(i_node = 1; i_node < (node_last+1); i_node++){
					
					PRF_CSEND_INT(i_node, &MPI_cells.number_of_host_cells_per_node[i_node], 1, node_zero);
					
					list_of_real = (double*)malloc(MPI_cells.number_of_host_cells_per_node[i_node] * sizeof(double));
					
					for(i_cell = 0; i_cell < MPI_cells.number_of_host_cells_per_node[i_node]; i_cell++){
						
						list_of_real[i_cell] = MPI_cells.data[MPI_cells.host2ext_index[i_list]];
						
						i_list++;
					}
					
					PRF_CSEND_REAL(i_node, list_of_real, MPI_cells.number_of_host_cells_per_node[i_node], node_zero);
					
					free(list_of_real);
				}			
			}
		}
		
		
		/* MPI.3 Node-0 collects updated data into MPI_cells.data */
		{
			if(myid == 0){
				
				/* Node-0 fills his own MPI_cells.data */
				
				size = MPI_cells.number_of_host_cells_per_node[0];
				
				list_of_real = (double*)malloc(size * sizeof(double));
				
				for(i_cell = 0; i_cell < size; i_cell++){
					
					list_of_real[i_cell] = MPI_cells.data[MPI_cells.host2ext_index[i_cell]];
				}
				
				for(i_list = 0; i_list < size; i_list++){
					
					cell_data[MPI_cells.hosting_cell_index[i_list]] += list_of_real[i_list];
				}
				
				for(i_list = 0; i_list < size; i_list++){
					
					list_of_real[i_list] = cell_data[MPI_cells.hosting_cell_index[i_list]];
				}
				
				for(i_cell = 0; i_cell < size; i_cell++){
						
					MPI_cells.data[MPI_cells.host2ext_index[i_cell]] = list_of_real[i_cell];
				}
				
				free(list_of_real);
				
				/* Node-0 collects MPI_cells.data from other nodes */
				
				i_list = MPI_cells.number_of_host_cells_per_node[0];
				
				for(i_node = 1; i_node < (node_last+1); i_node++){
				
					PRF_CRECV_INT(i_node, &size, 1, i_node);
					
					list_of_real = (double*)malloc(size * sizeof(double));
					
					PRF_CRECV_REAL(i_node, list_of_real, size, i_node);
					
					for(i_cell = 0; i_cell < size; i_cell++){
						
						MPI_cells.data[MPI_cells.host2ext_index[i_list]] = list_of_real[i_cell];
						
						i_list++;
					}
					
					free(list_of_real);
				}
			}

			/* other nodes sum up data on their host cells */
			if(myid > 0){

				PRF_CRECV_INT(node_zero, &size, 1, node_zero);
				
				list_of_real = (double*)malloc(size * sizeof(double));
				
				PRF_CRECV_REAL(node_zero, list_of_real, size, node_zero);
				
				for(i_list = 0; i_list < size; i_list++){
				
					cell_data[MPI_cells.hosting_cell_index[i_list]] += list_of_real[i_list];
				}
			
				for(i_list = 0; i_list < size; i_list++){
				
					list_of_real[i_list] = cell_data[MPI_cells.hosting_cell_index[i_list]];
				}

				PRF_CSEND_INT(node_zero, &size, 1, myid);
				
				PRF_CSEND_REAL(node_zero, list_of_real, size, myid);

				free(list_of_real);
			}
		}
		
		/* MPI.4 Node-0 sends updated data to nodes */
		{
			if(myid == 0){
				
				i_list = MPI_cells.number_of_ext_cells_per_node[0];
				
				for(i_node = 1; i_node < (node_last+1); i_node++){
					
					PRF_CSEND_INT(i_node, &MPI_cells.number_of_ext_cells_per_node[i_node], 1, node_zero);
					
					list_of_real = (double*)malloc(MPI_cells.number_of_ext_cells_per_node[i_node] * sizeof(double));
					
					for(i_cell = 0; i_cell < MPI_cells.number_of_ext_cells_per_node[i_node]; i_cell++){
						
						list_of_real[i_cell] = MPI_cells.data[i_list];
						
						i_list++;
					}
					
					PRF_CSEND_REAL(i_node, list_of_real, MPI_cells.number_of_ext_cells_per_node[i_node], node_zero);
					
					free(list_of_real);
				}
				
				i_list = 0;
				
				loop_ext_cells_ptr{
					
					cell_data[i_cell] = MPI_cells.data[i_list];
					
					i_list++;
				}
				
			}
			
			if(myid > 0){
				
				PRF_CRECV_INT(node_zero, &size, 1, node_zero);
				
				list_of_real = (double*)malloc(size * sizeof(double));
				
				PRF_CRECV_REAL(node_zero, list_of_real, size, node_zero);
				
				i_list = 0;
				
				loop_ext_cells_ptr{
					
					cell_data[i_cell] = list_of_real[i_list];
					
					i_list++;
				}
				
				free(list_of_real);
			}
		}
		
	}
#endif
#endif
		
#if 1 /* free */
#if RP_NODE

	void free_all_parallel(void)
	{
		if(MPI_cells.number_of_ext_cells_per_node != NULL){
			
			free(MPI_cells.number_of_ext_cells_per_node);
		}
		
		if(MPI_cells.number_of_host_cells_per_node != NULL){
			
			free(MPI_cells.number_of_host_cells_per_node);
		}
		
		if(MPI_cells.host_of_cell != NULL){
			
			free(MPI_cells.host_of_cell);
		}
		
		if(MPI_cells.host_of_ext_cells != NULL){
			
			free(MPI_cells.host_of_ext_cells);
		}
		
		if(MPI_cells.hosting_cell_index != NULL){
			
			free(MPI_cells.hosting_cell_index);
		}

		if(MPI_cells.host2ext_index != NULL){

			free(MPI_cells.host2ext_index);
		}			
		
		if(MPI_cells.data != NULL){
			
			free(MPI_cells.data);
		}

		if(MPI_faces.principal_face != NULL){
			
			free(MPI_faces.principal_face);
		}
	}

#endif
#endif

#endif