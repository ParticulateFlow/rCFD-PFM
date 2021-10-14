#include <stdio.h>

#include "rCFD_types.h"
#include "rCFD_globals.h"
#include "rCFD_parallel.h"
#include "rCFD_defaults.h"
#include "rCFD_macros.h"
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
    
    rCFD_read_C2Cs
    
    rCFD_prepare_C2Cs_for_MP
    
    rCFD_run
    
    rCFD_free_all
*/

/*************************************************************************************/
DEFINE_ON_DEMAND(rCFD_init_all)
/*************************************************************************************/
{
    /* D1. Solver_Dict & Solver */
    {
        rCFD_default_Solver_Dict(&Solver_Dict);
        
        rCFD_user_set_Solver_Dict(&Solver_Dict);
        
#if RP_NODE     
        rCFD_default_Solver(&Solver);
#endif      
    }

    /* D2. File_Dict */
    {
        rCFD_default_File_Dict(&Solver_Dict, &File_Dict);
        
        rCFD_user_set_File_Dict(&Solver_Dict, &File_Dict);
    }

    /* D3. Phase_Dict */
    {
#if RP_NODE     
        Phase_Dict = (Phase_Dict_type*)malloc(Solver_Dict.number_of_phases * sizeof(Phase_Dict_type));
        
        rCFD_default_Phase_Dict(&Solver_Dict, Phase_Dict);
        
        rCFD_user_set_Phase_Dict(&Solver_Dict, Phase_Dict); 
#endif      
    }       
    
    /* D4. Tracer_Dict */
    {
#if RP_NODE
        Tracer_Dict.random_walk = (short*)malloc(Solver_Dict.number_of_phases * sizeof(short));
        
        rCFD_default_Tracer_Dict(&Solver_Dict, &Tracer_Dict);

        rCFD_user_set_Tracer_Dict(&Solver_Dict, &Tracer_Dict);
#endif      
    }
    
    /* D5. Norm_Dict */
    {
#if RP_NODE
        rCFD_default_Norm_Dict(&Solver_Dict, &Norm_Dict);

        rCFD_user_set_Norm_Dict(&Solver_Dict, &Norm_Dict);
#endif      
    }   

    /* D6. Rec_Dict */
    {
        rCFD_default_Rec_Dict(&Solver_Dict, &Rec_Dict);

        rCFD_user_set_Rec_Dict(&Solver_Dict, &Rec_Dict);    
    }   

    /* D7. Data_Dict */	
	{
#if RP_NODE     
        int i_phase;
        
        Data_Dict = (Data_Dict_type**)malloc(Solver_Dict.number_of_phases * sizeof(Data_Dict_type*));
        
        loop_phases{
            
            Data_Dict[i_phase] = (Data_Dict_type*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(Data_Dict_type));
        }   
        
        rCFD_default_Data_Dict(&Solver_Dict, Phase_Dict, Data_Dict);
        
        rCFD_user_set_Data_Dict(&Solver_Dict, Phase_Dict, Data_Dict);
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
        
        rCFD_default_Balance_Dict(&Solver_Dict, Phase_Dict, Balance_Dict);
        
        rCFD_user_set_Balance_Dict(&Solver_Dict, Phase_Dict, Balance_Dict);
#endif      
    }

    /* D9. Cell_Dict */
    {
#if RP_NODE
    
        rCFD_default_Cell_Dict(&Solver_Dict, &Cell_Dict);

        rCFD_user_set_Cell_Dict(&Solver_Dict, &Cell_Dict);
#endif      
    }   

    /* D10. Face_Dict */
    {
#if RP_NODE
        rCFD_default_Face_Dict(&Solver_Dict, &Face_Dict);

        rCFD_user_set_Face_Dict(&Solver_Dict, &Face_Dict);
#endif      
    }   
        
    /* G1. Cells (first initialization) */
    {
#if RP_NODE     
        int     i_phase, i_frame, i_cell;
        
        C.x = (double**)malloc(Cell_Dict.number_of_cells * sizeof(double*));
        
        loop_cells{
            
            C.x[i_cell] = (double*)malloc( 3 * sizeof(double));
        }
        
        C.volume = (double*)malloc(Cell_Dict.number_of_cells * sizeof(double));
        
        C.average_velocity =    (double**)malloc(Solver_Dict.number_of_phases * sizeof(double*));
        C.crossing_time =       (double**)malloc(Solver_Dict.number_of_phases * sizeof(double*));
        
        loop_phases{
            
            C.average_velocity[i_phase] =   (double*)malloc(Cell_Dict.number_of_cells * sizeof(double));            
            C.crossing_time[i_phase] =      (double*)malloc(Cell_Dict.number_of_cells * sizeof(double));
        }

        C.hit_by_other_cell =   (short*)malloc(Cell_Dict.number_of_cells * sizeof(short));

        C.island_id =           (short*)malloc(Cell_Dict.number_of_cells * sizeof(short));
        
        C.weight_after_shift =  (double*)malloc(Cell_Dict.number_of_cells * sizeof(double));
        C.weight_after_swap =   (double*)malloc(Cell_Dict.number_of_cells * sizeof(double));

        C.vof = (double***)malloc(Solver_Dict.number_of_frames * sizeof(double**));
        
        loop_frames{
            
            C.vof[i_frame] = (double**)malloc(Cell_Dict.number_of_cells * sizeof(double*));
            
            loop_cells{
                
                C.vof[i_frame][i_cell] = (double*)malloc(Solver_Dict.number_of_phases * sizeof(double));
            }
        }

        C.data =        (double***)malloc(Solver_Dict.number_of_phases * sizeof(double**));
        C.data_shift =  (double***)malloc(Solver_Dict.number_of_phases * sizeof(double**));
        C.data_swap =   (double***)malloc(Solver_Dict.number_of_phases * sizeof(double**));
        
        loop_phases{
            
            C.data[i_phase] =       (double**)malloc(Cell_Dict.number_of_cells * sizeof(double*));
            C.data_shift[i_phase] = (double**)malloc(Cell_Dict.number_of_cells * sizeof(double*));
            C.data_swap[i_phase] =  (double**)malloc(Cell_Dict.number_of_cells * sizeof(double*));
            
            loop_cells{
                
                C.data[i_phase][i_cell] =       (double*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(double));
                C.data_shift[i_phase][i_cell] = (double*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(double));
                C.data_swap[i_phase][i_cell] =  (double*)malloc(Phase_Dict[i_phase].number_of_data * sizeof(double));
            }
        }       
        
        if(Solver_Dict.data_drifting_on){
            
            C.drift_exchange = (double*)malloc(Cell_Dict.number_of_cells * sizeof(double));
        }
        else{
            
            C.drift_exchange = NULL;
        }       
        
        if(Cell_Dict.number_of_user_vars > 0){
            
            C.user = (double**)malloc(Cell_Dict.number_of_cells * sizeof(double*));
            
            loop_cells{
                
                C.user[i_cell] = (double*)malloc(Cell_Dict.number_of_user_vars * sizeof(double));
            }
        }
        else{
            
            C.user = NULL;
        }
                
        rCFD_default_Cell(&Solver_Dict, Phase_Dict, &Cell_Dict, &C);
#endif  
    }   

    /* G2. Faces */
    {
#if RP_NODE
        int i_face;
        
        F.c0 = (int*)malloc(Face_Dict.number_of_faces * sizeof(int));
        
        F.c1 = (int*)malloc(Face_Dict.number_of_faces * sizeof(int));
        
        F.area = (double**)malloc(Face_Dict.number_of_faces * sizeof(double*));
        
        loop_faces{
            
            F.area[i_face] = (double*)malloc( 3 * sizeof(double));
        }

        rCFD_default_Face(&Solver_Dict, &Cell_Dict, &Face_Dict, &F);
#endif      
    }
    
    /* G3. Tracer */
    {
#if RP_NODE

        Tracer.monitor_counter = (int*)malloc(Solver_Dict.number_of_phases * sizeof(int));

        Tracer.number_of_shifts = (int*)malloc(Solver_Dict.number_of_phases * sizeof(int));
        
        Tracer.shifts = (C2C_shift_type**)malloc(Solver_Dict.number_of_phases * sizeof(C2C_shift_type*));
        
        rCFD_default_Tracer(&Solver_Dict, &Tracer);

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
        
        rCFD_default_Norms(&Solver_Dict, &Norms);
        
#endif      
    }

    /* G5. Rec */
    {
        int     i_state, i_state2, i_island;
        
        Rec.global_frame = (int*)malloc(Solver_Dict.number_of_islands * sizeof(int));
        
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
        
        rCFD_default_Rec(&Solver_Dict, &File_Dict, &Rec_Dict, &Rec);
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
                
                Balance[i_phase][i_data].mass_error = 0.0;

                Balance[i_phase][i_data].mass_error_prev = 0.0;

                Balance[i_phase][i_data].face_swap = Solver_Dict.face_swap_min;

                Balance[i_phase][i_data].face_swap_loops = 1;               
            }           
        }
#endif      
    }

    /* G1 + G6: Init data fields and update balances */
    {
#if RP_NODE     
        
        rCFD_user_init_Data(&Solver_Dict, Balance, Phase_Dict, Data_Dict, &Cell_Dict, &C);
        
        int i_frame = 0, i_phase, i_data, i_cell;
        
        loop_int_cells{
            
            loop_phases{
            
                loop_data{
                    
                    if(Data_Dict[i_phase][i_data].type == temperature_data){
                        
                        Balance[i_phase][i_data].mass_integral += (C.data[i_phase][i_cell][i_data] - Phase_Dict[i_phase].reference_temperature) * 
                    
                            Phase_Dict[i_phase].heat_capacity * Phase_Dict[i_phase].density * C.volume[i_cell] * C.vof[i_frame][i_cell][i_phase];
                    }
                    else{
                        
                        Balance[i_phase][i_data].mass_integral += C.data[i_phase][i_cell][i_data] * 
                    
                            Phase_Dict[i_phase].density * C.volume[i_cell] * C.vof[i_frame][i_cell][i_phase];
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
        
        rCFD_user_post(&Solver_Dict, Phase_Dict, &Cell_Dict, &C);
#endif
    }       

    /* P: Parallel grid communication */
    {
#if RP_NODE
        init_parallel_grid(&Solver_Dict, &Cell_Dict, &Face_Dict);
#endif
    }       
    
}

/*************************************************************************************/
DEFINE_ON_DEMAND(rCFD_read_C2Cs)
/*************************************************************************************/
{
#if RP_NODE 
    int     i_state, i_phase, i_frame, i_C2C, i_read, i_read_int, i_read_incoming, i_island;
    
    int     C2C_index;
    
    int     *i_island_int, *i_island_incoming;

    int     total_number_of_C2Cs_read = 0;

    C2C_type    tmp_C2Cs;
    
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
                        C2Cs[current_pattern].format = c0_n0_c1_n1_w0;
                        
                        C2Cs[current_pattern].number_of_shifts =        0;
                        C2Cs[current_pattern].number_of_shifts_in =     0;
                        C2Cs[current_pattern].number_of_shifts_out =    0;
                        
                        C2Cs[current_pattern].shifts = NULL;
                        C2Cs[current_pattern].shifts_in = NULL;
                        C2Cs[current_pattern].shifts_out = NULL;
                        
                        C2Cs[current_pattern].island_offsets =      (int*)malloc((Solver_Dict.number_of_islands+1) * sizeof(int));
                        C2Cs[current_pattern].island_offsets_in =   (int*)malloc((Solver_Dict.number_of_islands+1) * sizeof(int));

                        loop_islands{
                            
                            C2Cs[current_pattern].island_offsets[i_island] =        0;
                            C2Cs[current_pattern].island_offsets_in[i_island] =     0;
                        }
                        
                        C2Cs[current_pattern].island_offsets[(i_island + 1)] =      0;
                        C2Cs[current_pattern].island_offsets_in[(i_island + 1)] =   0;
                        
                        C2Cs[current_pattern].number_of_MPI_shifts = 0;
                        C2Cs[current_pattern].number_of_shifts_to_node_zero = NULL;
                        C2Cs[current_pattern].number_of_shifts_from_node_zero = NULL;
                        C2Cs[current_pattern].in2out = NULL;                        
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
        int     i_tmp, c0, node0, c1, node1, my_island_id;
        double  tmp, w0;

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
                                
                                for(i_C2C = 0; i_C2C < tmp_C2Cs.number_of_shifts; i_C2C++){
                                    
                                    tmp_C2Cs.shifts[i_C2C].c0 =     -1;
                                    tmp_C2Cs.shifts[i_C2C].node0 =  -1;
                                    tmp_C2Cs.shifts[i_C2C].c1 =     -1;
                                    tmp_C2Cs.shifts[i_C2C].node1 =  -1;
                                    tmp_C2Cs.shifts[i_C2C].w0 =     0.0;
                                }
                            
                                tmp_C2Cs.island_offsets =       (int*)malloc((Solver_Dict.number_of_islands+1) * sizeof(int));
                                tmp_C2Cs.island_offsets_in =    (int*)malloc((Solver_Dict.number_of_islands+1) * sizeof(int));

                                loop_islands{
                                    
                                    tmp_C2Cs.island_offsets[i_island] =     0;
                                    tmp_C2Cs.island_offsets_in[i_island] =  0;
                                }
                            
                            }
                            
                            /* fill tmp_C2Cs and monitor indices */
                            for(i_C2C = 0; i_C2C < tmp_C2Cs.number_of_shifts; i_C2C++){
              
                                if((i_C2C % Solver_Dict.C2C_loading_reduction) == 0){
                      
                                    if(tmp_C2Cs.format == c0_n0_c1_n1_w0){

                                        fscanf(fi,"%d %d %d %d %le\n", &c0, &node0, &c1, &node1, &w0);
                                    }
                                
                                    tmp_C2Cs.shifts[i_read].c0 = c0;
                                    tmp_C2Cs.shifts[i_read].node0 = node0;
                                    tmp_C2Cs.shifts[i_read].c1 = c1;
                                    tmp_C2Cs.shifts[i_read].node1 = node1;
                                    tmp_C2Cs.shifts[i_read].w0 = w0;
                                                                        
                                    i_read++;
                                
                                    my_island_id = C.island_id[c1];
                                    
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
                                C2Cs[current_pattern].number_of_shifts = i_read_int;
                                
                                C2Cs[current_pattern].shifts = (C2C_shift_type*)malloc(i_read_int * sizeof(C2C_shift_type));
                                
                                for(i_C2C=0; i_C2C < C2Cs[current_pattern].number_of_shifts; i_C2C++){
                                    
                                    C2Cs[current_pattern].shifts[i_C2C].c0 =        -1;
                                    C2Cs[current_pattern].shifts[i_C2C].node0 =     -1;
                                    C2Cs[current_pattern].shifts[i_C2C].c1 =        -1;
                                    C2Cs[current_pattern].shifts[i_C2C].node1 =     -1;
                                    C2Cs[current_pattern].shifts[i_C2C].w0 =        0.0;
                                }

                                C2Cs[current_pattern].number_of_shifts_in = i_read_incoming;

                                C2Cs[current_pattern].shifts_in = (C2C_shift_type*)malloc(i_read_incoming * sizeof(C2C_shift_type));

                                for(i_C2C=0; i_C2C < C2Cs[current_pattern].number_of_shifts_in; i_C2C++){
                                    
                                    C2Cs[current_pattern].shifts_in[i_C2C].c0 =     -1;
                                    C2Cs[current_pattern].shifts_in[i_C2C].node0 =  -1;
                                    C2Cs[current_pattern].shifts_in[i_C2C].c1 =     -1;
                                    C2Cs[current_pattern].shifts_in[i_C2C].node1 =  -1;
                                    C2Cs[current_pattern].shifts_in[i_C2C].w0 =     0.0;
                                }
                                
                                total_number_of_C2Cs_read += 
                                
                                    C2Cs[current_pattern].number_of_shifts + C2Cs[current_pattern].number_of_shifts_in;
                            }

                            /* B.4.3. fill C2Cs.shifts/shifts_in by tmp_C2Cs*/
                            {
                                for(i_C2C = 0; i_C2C < tmp_C2Cs.number_of_shifts; i_C2C++){ 
 
                                    my_island_id = C.island_id[tmp_C2Cs.shifts[i_C2C].c1];
                                    
                                    node0 = tmp_C2Cs.shifts[i_C2C].node0;
                                    node1 = tmp_C2Cs.shifts[i_C2C].node1;
                                    
                            
                                    if(node0 == node1){ 
                                    
                                        C2C_index = tmp_C2Cs.island_offsets[my_island_id] + i_island_int[my_island_id];
                                        
                                        C2Cs[current_pattern].shifts[C2C_index].c0 =    tmp_C2Cs.shifts[i_C2C].c0;
                                        C2Cs[current_pattern].shifts[C2C_index].node0 = tmp_C2Cs.shifts[i_C2C].node0;
                                        C2Cs[current_pattern].shifts[C2C_index].c1 =    tmp_C2Cs.shifts[i_C2C].c1;
                                        C2Cs[current_pattern].shifts[C2C_index].node1 = tmp_C2Cs.shifts[i_C2C].node1;
                                        C2Cs[current_pattern].shifts[C2C_index].w0 =    tmp_C2Cs.shifts[i_C2C].w0;
                                    
                                        i_island_int[my_island_id]++;                       
                                    }
                                    else{  
                                
                                        C2C_index = tmp_C2Cs.island_offsets_in[my_island_id] + i_island_incoming[my_island_id];
                                    
                                        C2Cs[current_pattern].shifts_in[C2C_index].c0 =     tmp_C2Cs.shifts[i_C2C].c0;
                                        C2Cs[current_pattern].shifts_in[C2C_index].node0 =  tmp_C2Cs.shifts[i_C2C].node0;
                                        C2Cs[current_pattern].shifts_in[C2C_index].c1 =     tmp_C2Cs.shifts[i_C2C].c1;
                                        C2Cs[current_pattern].shifts_in[C2C_index].node1 =  tmp_C2Cs.shifts[i_C2C].node1;
                                        C2Cs[current_pattern].shifts_in[C2C_index].w0 =     tmp_C2Cs.shifts[i_C2C].w0;
                                        
                                        i_island_incoming[my_island_id]++;
                                    }                               
                                }
                                
                                C2Cs[current_pattern].island_offsets[0] = 0;
                                C2Cs[current_pattern].island_offsets_in[0] = 0;
                                
                                loop_islands{
                                
                                    C2Cs[current_pattern].island_offsets[(i_island + 1)] = tmp_C2Cs.island_offsets[(i_island + 1)];
                                    C2Cs[current_pattern].island_offsets_in[(i_island + 1)] = tmp_C2Cs.island_offsets_in[(i_island + 1)];
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

    /* D. Message & Transcript */
    {
#if RP_NODE
    
        total_number_of_C2Cs_read = PRF_GISUM1(total_number_of_C2Cs_read);
            
        if((myid == 0) && (Solver_Dict.verbal == 1)){
                
            FILE    *f_out = fopen("./Run.trn", "a");
            
            if(f_out){
                
                fprintf(f_out,"\n\nrCFD_read_C2Cs");
            
                fprintf(f_out,"\n\n   Read %d C2Cs from %s in %f seconds", 
                
                    total_number_of_C2Cs_read, File_Dict.C2C_filename, 
                    
                    (double)(clock() - Solver.clock)/(double)CLOCKS_PER_SEC);                                           

                fclose(f_out);
            }       
        }
        
        Message0("\n...rCFD_read_C2Cs ->  Read %d C2Cs from %s in %f seconds", 
        
            total_number_of_C2Cs_read, File_Dict.C2C_filename, 
            
            (double)(clock() - Solver.clock)/(double)CLOCKS_PER_SEC);
#endif
    }
}
    
/*************************************************************************************/
DEFINE_ON_DEMAND(rCFD_prepare_C2Cs_for_MPI)
/*************************************************************************************/
{
#if RP_NODE

    int                 i_state, i_phase, i_frame, i_out, i_node, i_C2C;
    int                 size, max_size, index_in_total_list;
    int                 *list_of_int = NULL;
    double              *list_of_real = NULL;
    
    C2C_type            tmp_C2Cs;
    
    Solver.clock = clock();
    
    /* A. allocate C2C_MPI.shifts (shifts.data is allocated later) */
    {
        
        max_size = 0;
        
        loop_states{
        
            loop_phases{
                
                loop_frames{
                    
                    size = PRF_GISUM1(C2Cs[current_pattern].number_of_shifts_in);
                    
                    if(size > max_size){
                        
                        max_size = size;
                    }
                }
            }
        }
        
        if(myid == 0){
            
            C2Cs_MPI.shifts_in = (C2C_MPI_shift_type*)malloc(max_size * sizeof(C2C_MPI_shift_type));
        }
        else{
            C2Cs_MPI.shifts_in = NULL;
        }           
        
        C2Cs_MPI.max_number_of_MPI_shifts = max_size;

        C2Cs_MPI.number_of_shifts_in = 0;
    }

    /* B. determine C2Cs.shift_out, C2Cs.in2out, C2Cs.number_of_shifts_to_node_zero */
    {
        
        loop_states{
            
            loop_phases{
        
                loop_frames{
                
                    /* B.1: send local rCFD_C2C_incoming entries to Node-0,      */
                    /*      which collects them into rCFD_C2C_MPI_buffer         */
                    {
                    
                        if (myid > 0){
                
                            size = C2Cs[current_pattern].number_of_shifts_in;
                            
                            PRF_CSEND_INT(node_zero, &size, 1, myid);

                            if(size > 0){
                                
                                list_of_int=(int*)malloc(size * sizeof(int));
                                list_of_real=(real*)malloc(size * sizeof(double));
                                
                                loop_C2Cs_size{ 
                                
                                    list_of_int[i_C2C] = C2Cs[current_pattern].shifts_in[i_C2C].c0;
                                }
                                
                                PRF_CSEND_INT(node_zero, list_of_int, size, myid);
                                
                                loop_C2Cs_size{ 
                                
                                    list_of_int[i_C2C] = C2Cs[current_pattern].shifts_in[i_C2C].node0;
                                }
                                
                                PRF_CSEND_INT(node_zero, list_of_int, size, myid);
                                
                                loop_C2Cs_size{ 
                                
                                    list_of_int[i_C2C] = C2Cs[current_pattern].shifts_in[i_C2C].c1;
                                }
                                
                                PRF_CSEND_INT(node_zero, list_of_int, size, myid);
                                
                                loop_C2Cs_size{ 
                                
                                    list_of_int[i_C2C] = C2Cs[current_pattern].shifts_in[i_C2C].node1;
                                }
                                
                                PRF_CSEND_INT(node_zero, list_of_int, size, myid);
                                
                                loop_C2Cs_size{ 
                                
                                    list_of_real[i_C2C] = C2Cs[current_pattern].shifts_in[i_C2C].w0;
                                }
                                
                                PRF_CSEND_REAL(node_zero, list_of_real, size, myid);
                                
                                free(list_of_int);
                                free(list_of_real);
                            }
                        }
              
                        if (myid == 0){
                
                            /* initiate for each frame */
                            for(i_C2C = 0; i_C2C < C2Cs[current_pattern].number_of_shifts_in; i_C2C++){
                        
                                C2Cs_MPI.shifts_in[i_C2C].c0=    C2Cs[current_pattern].shifts_in[i_C2C].c0;
                                C2Cs_MPI.shifts_in[i_C2C].node0= C2Cs[current_pattern].shifts_in[i_C2C].node0;
                                C2Cs_MPI.shifts_in[i_C2C].c1=    C2Cs[current_pattern].shifts_in[i_C2C].c1;
                                C2Cs_MPI.shifts_in[i_C2C].node1= C2Cs[current_pattern].shifts_in[i_C2C].node1;
                                C2Cs_MPI.shifts_in[i_C2C].w0=    C2Cs[current_pattern].shifts_in[i_C2C].w0;
                            }
                  
                            index_in_total_list = C2Cs[current_pattern].number_of_shifts_in;    /* first index from Node-0 */
                            
                            C2Cs[current_pattern].number_of_shifts_to_node_zero = (int*)malloc((node_last + 1) * sizeof(int));
                            
                            C2Cs[current_pattern].number_of_shifts_to_node_zero[0] = C2Cs[current_pattern].number_of_shifts_in;
                            
                            for(i_node = 1; i_node < (node_last+1); i_node++){
                                
                                PRF_CRECV_INT(i_node, &size, 1, i_node);
                        
                                if(size > 0){
                      
                                    list_of_int=(int*)malloc(size * sizeof(int));
                                    list_of_real=(real*)malloc(size * sizeof(real));
                            
                                    PRF_CRECV_INT(i_node, list_of_int, size, i_node); 
                                    
                                    loop_C2Cs_size{ 
                                    
                                        C2Cs_MPI.shifts_in[index_in_total_list + i_C2C].c0 = list_of_int[i_C2C]; 
                                    }
                            
                                    PRF_CRECV_INT(i_node, list_of_int, size, i_node); 
                                    
                                    loop_C2Cs_size{ 
                                    
                                        C2Cs_MPI.shifts_in[index_in_total_list + i_C2C].node0 = list_of_int[i_C2C]; 
                                    }
                            
                                    PRF_CRECV_INT(i_node, list_of_int, size, i_node); 
                                    
                                    loop_C2Cs_size{ 
                                    
                                        C2Cs_MPI.shifts_in[index_in_total_list + i_C2C].c1 = list_of_int[i_C2C]; 
                                    }
                            
                                    PRF_CRECV_INT(i_node, list_of_int, size, i_node); 
                                    
                                    loop_C2Cs_size{ 
                                    
                                        C2Cs_MPI.shifts_in[index_in_total_list + i_C2C].node1 = list_of_int[i_C2C]; 
                                    }
                            
                                    PRF_CRECV_REAL(i_node, list_of_real, size, i_node); 
                                    
                                    loop_C2Cs_size{ 
                                    
                                        C2Cs_MPI.shifts_in[index_in_total_list + i_C2C].w0 = list_of_real[i_C2C]; 
                                    }
                                    
                                    free(list_of_int);
                                    free(list_of_real);
                                }
                                
                                C2Cs[current_pattern].number_of_shifts_to_node_zero[i_node] = size;

                                index_in_total_list += size;
                            }
                            
                            C2Cs_MPI.number_of_shifts_in = index_in_total_list;
                            
                            C2Cs[current_pattern].number_of_MPI_shifts = index_in_total_list;
                        }
                    }       
                                
                    /* B.2: Node-0 defines and sends C2Cs.shifts_out for each Node */
                    {

                        if (myid == 0){
                        
                            i_out = 0;          
                        
                            for(i_C2C = 0; i_C2C < C2Cs_MPI.number_of_shifts_in; i_C2C++){
                                    
                                if(C2Cs_MPI.shifts_in[i_C2C].node0 == node_zero){
                            
                                    i_out++;
                                }
                            }
                            
                            C2Cs[current_pattern].number_of_shifts_out = i_out;
                            
                            C2Cs[current_pattern].shifts_out = (C2C_shift_type*)malloc(i_out * sizeof(C2C_shift_type));
                            
                            C2Cs[current_pattern].number_of_shifts_from_node_zero = (int*)malloc((node_last + 1) * sizeof(int));
                            
                            C2Cs[current_pattern].number_of_shifts_from_node_zero[0] = i_out;

                            
                            i_out = 0;
                            
                            /* fill C2C.shifts_out */
                            for(i_C2C = 0; i_C2C < C2Cs_MPI.number_of_shifts_in; i_C2C++){
                                    
                                if(C2Cs_MPI.shifts_in[i_C2C].node0 == node_zero){
                            
                                    C2Cs[current_pattern].shifts_out[i_out].c0 =    C2Cs_MPI.shifts_in[i_C2C].c0;
                                    C2Cs[current_pattern].shifts_out[i_out].node0 = C2Cs_MPI.shifts_in[i_C2C].node0;
                                    C2Cs[current_pattern].shifts_out[i_out].c1 =    C2Cs_MPI.shifts_in[i_C2C].c1;
                                    C2Cs[current_pattern].shifts_out[i_out].node1 = C2Cs_MPI.shifts_in[i_C2C].node1;
                                    C2Cs[current_pattern].shifts_out[i_out].w0 =    C2Cs_MPI.shifts_in[i_C2C].w0;
                                    
                                    i_out++;
                                }
                            }
        
                            /* loop other nodes */
                            for(i_node = 1; i_node < (node_last+1); i_node++){
                            
                                i_out = 0;
                                
                                /* get number of outgoing tracers */                            
                                for(i_C2C = 0; i_C2C < C2Cs_MPI.number_of_shifts_in; i_C2C++){
                        
                                    if(C2Cs_MPI.shifts_in[i_C2C].node0 == i_node){
                            
                                        i_out++;
                                    }
                                }
                                
                                C2Cs[current_pattern].number_of_shifts_from_node_zero[i_node] = i_out;
                        
                                /* wrong !! C2Cs[current_pattern].number_of_shifts_out = i_out;*/
                                
                                size = i_out;
                                
                                /* send outgoing tracers to individual nodes */ 
                                PRF_CSEND_INT(i_node, &size, 1, node_zero);
                                
                                if (size > 0){
                            
                                    /* allocate tmp storage for outgoing tracers */
                                    tmp_C2Cs.shifts_out = (C2C_shift_type*)malloc(size * sizeof(C2C_shift_type));
                                
                                    /* fill tmp storage for outgoing tracers */
                                    i_out = 0;
                                    
                                    for(i_C2C = 0; i_C2C < C2Cs_MPI.number_of_shifts_in; i_C2C++){
                            
                                        if(C2Cs_MPI.shifts_in[i_C2C].node0 == i_node){
                            
                                            tmp_C2Cs.shifts_out[i_out].c0 =    C2Cs_MPI.shifts_in[i_C2C].c0;
                                            tmp_C2Cs.shifts_out[i_out].node0 = C2Cs_MPI.shifts_in[i_C2C].node0;
                                            tmp_C2Cs.shifts_out[i_out].c1 =    C2Cs_MPI.shifts_in[i_C2C].c1;
                                            tmp_C2Cs.shifts_out[i_out].node1 = C2Cs_MPI.shifts_in[i_C2C].node1;
                                            tmp_C2Cs.shifts_out[i_out].w0 =    C2Cs_MPI.shifts_in[i_C2C].w0;

                                            i_out++;
                                        }
                                    }
                            
                                    list_of_int=(int*)malloc(size * sizeof(int));
                                    list_of_real=(real*)malloc(size * sizeof(real));
                            
                                    loop_C2Cs_size{ 
                                    
                                        list_of_int[i_C2C] = tmp_C2Cs.shifts_out[i_C2C].c0;
                                    }
                                    
                                    PRF_CSEND_INT(i_node, list_of_int, size, node_zero);
                            
                                    loop_C2Cs_size{ 
                                    
                                        list_of_int[i_C2C] = tmp_C2Cs.shifts_out[i_C2C].node0;
                                    }
                                    
                                    PRF_CSEND_INT(i_node, list_of_int, size, node_zero);
                            
                                    loop_C2Cs_size{ 
                                    
                                        list_of_int[i_C2C] = tmp_C2Cs.shifts_out[i_C2C].c1;
                                    }
                                    
                                    PRF_CSEND_INT(i_node, list_of_int, size, node_zero);
                            
                                    loop_C2Cs_size{ 
                                    
                                        list_of_int[i_C2C] = tmp_C2Cs.shifts_out[i_C2C].node1;
                                    }
                            
                                    PRF_CSEND_INT(i_node, list_of_int, size, node_zero);
                            
                                    loop_C2Cs_size{
                                        
                                        list_of_real[i_C2C] = tmp_C2Cs.shifts_out[i_C2C].w0;
                                    }
                                    
                                    PRF_CSEND_REAL(i_node, list_of_real, size, node_zero);

                                    free(list_of_int); 
                                    free(list_of_real);
                                
                                    free(tmp_C2Cs.shifts_out);  
                                }
                                        
                            } /* loop other nodes */
                    
                        } /* node == 0 */
                    
                        if (myid > 0){

                            /* receive rCFD_Tracer_outgoing from Node-0 */
             
                            PRF_CRECV_INT(node_zero, &size, 1, node_zero);
                            
                            C2Cs[current_pattern].number_of_shifts_out = size;
                            
                            if(size > 0){
                        
                                C2Cs[current_pattern].shifts_out = (C2C_shift_type*)malloc(size * sizeof(C2C_shift_type));
                        
                                list_of_int=(int*)malloc(size * sizeof(int));
                                list_of_real=(real*)malloc(size * sizeof(real));
                                
                                PRF_CRECV_INT(node_zero, list_of_int, size, node_zero);

                                loop_C2Cs_size{
                                    
                                    C2Cs[current_pattern].shifts_out[i_C2C].c0=list_of_int[i_C2C];
                                }
                        
                                PRF_CRECV_INT(node_zero, list_of_int, size, node_zero);

                                loop_C2Cs_size{
                                    
                                    C2Cs[current_pattern].shifts_out[i_C2C].node0=list_of_int[i_C2C];
                                }
                        
                                PRF_CRECV_INT(node_zero, list_of_int, size, node_zero);

                                loop_C2Cs_size{
                                    
                                    C2Cs[current_pattern].shifts_out[i_C2C].c1=list_of_int[i_C2C];
                                }
                        
                                PRF_CRECV_INT(node_zero, list_of_int, size, node_zero);

                                loop_C2Cs_size{
                                    
                                    C2Cs[current_pattern].shifts_out[i_C2C].node1=list_of_int[i_C2C];
                                }
                        
                                PRF_CRECV_REAL(node_zero, list_of_real, size, node_zero);

                                loop_C2Cs_size{
                                    
                                    C2Cs[current_pattern].shifts_out[i_C2C].w0=list_of_real[i_C2C];
                                }
                                
                                free(list_of_int);
                                free(list_of_real); 
                                
                            } 
                            
                        } /* node > 0 */

                    }

                    /* B.3: Node-0 link outgoing data to list of incoming data */
                    /***********************************************************/
                    {
                        
                        if(myid == 0){
                            
                            C2Cs[current_pattern].in2out = (int**)malloc((node_last+1) * sizeof(int*));
                            
                            for(i_node = 0; i_node < (node_last+1); i_node++){
                                
                                C2Cs[current_pattern].in2out[i_node] = (int*)malloc(C2Cs[current_pattern].number_of_shifts_from_node_zero[i_node] * sizeof(int));
                            }
                            
                        
                            for(i_node = 0; i_node < (node_last+1); i_node++){
                
                                i_out = 0;
                            
                                for(i_C2C = 0; i_C2C < C2Cs_MPI.number_of_shifts_in; i_C2C++){
                                
                                    if(C2Cs_MPI.shifts_in[i_C2C].node0 == i_node){
                                                    
                                        C2Cs[current_pattern].in2out[i_node][i_out] = i_C2C;
                            
                                        i_out++;
                                    }
                                }
                            }
                        }
                    }
                    
                } /* loop frames */

            } /* loop phases */
        
        } /* loop states */
        
    }

    /* C. Allocate C2C_MPI_shifts.data */
    {
#if RP_NODE     
        int     max_number_of_data = 0;
        
        loop_phases{
            
            if(max_number_of_data < Phase_Dict[i_phase].number_of_data){

                max_number_of_data = Phase_Dict[i_phase].number_of_data;
            }           
        }
        
        if(myid == 0){
            
            size = C2Cs_MPI.max_number_of_MPI_shifts;
            
            loop_C2Cs_size{
                
                C2Cs_MPI.shifts_in[i_C2C].data = (double*)malloc(max_number_of_data * sizeof(double));
            }
        }
#endif      
    }   

#if 0
    /* D. allocate MPI_cells */
    {
#if RP_NODE
    
    MPI_cells.host_of_cell = (int*)malloc(Cell_Dict.number_of_cells * sizeof(int));

    MPI_cells.hosting_cell_index = NULL;

    if(myid > 0){
        
        PRF_CSEND_INT(node_zero, &Cell_Dict.number_of_ext_cells, 1, myid);
        
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
        
        MPI_cells.number_of_ext_cells_per_node[0] = Cell_Dict.number_of_ext_cells;
        
        /* get number_of_ext_cells (sum up all nodes) */
        MPI_cells.number_of_ext_cells = Cell_Dict.number_of_ext_cells;
        
        for(i_node = 1; i_node < (node_last+1); i_node++){
            
            PRF_CRECV_INT(i_node, &MPI_cells.number_of_ext_cells_per_node[i_node], 1, i_node);
            
            MPI_cells.number_of_ext_cells += MPI_cells.number_of_ext_cells_per_node[i_node];
        }
        
        MPI_cells.host_of_ext_cells = (int*)malloc(MPI_cells.number_of_ext_cells * sizeof(int));

        MPI_cells.host2ext_index = (int*)malloc(MPI_cells.number_of_ext_cells * sizeof(int));   
        
        MPI_cells.data = (double*)malloc(MPI_cells.number_of_ext_cells * sizeof(double));
    }
    
    rCFD_default_MPI_Cells(&Solver_Dict, &MPI_cells);   
    
#endif      
    }
    
    /* E. MPI_cells communication */
    {
#if RP_NODE     
        double *tmp_x;
        
        /* E.1 fill MPI_cells.host_of_ext_cells, tmp_x */
        {
            if(myid > 0){
                
                list_of_int = (int*)malloc( Cell_Dict.number_of_ext_cells * sizeof(int));
                
                list_of_real = (double*)malloc( 3 * Cell_Dict.number_of_ext_cells * sizeof(double));
                
                i_list = 0;
                
                loop_ext_cells{
                    
                    list_of_int[i_list] = MPI_cells.host_of_cell[i_cell];
                    
                    list_of_real[(3 * i_list)] =        C.x[i_cell][0];
                    list_of_real[(3 * i_list + 1)] =    C.x[i_cell][1];
                    list_of_real[(3 * i_list + 2)] =    C.x[i_cell][2];
                    
                    i_list++;
                }
                
                PRF_CSEND_INT(node_zero, list_of_int, Cell_Dict.number_of_ext_cells, myid);

                PRF_CSEND_REAL(node_zero, list_of_real, (3 * Cell_Dict.number_of_ext_cells), myid);
                
                free(list_of_int);
                
                free(list_of_real);
            }
            
            if(myid == 0){
                            
                tmp_x = (double*)malloc(3 * MPI_cells.number_of_ext_cells * sizeof(double));

                i_list = 0;
                
                loop_ext_cells{
                    
                    MPI_cells.host_of_ext_cells[i_list] = MPI_cells.host_of_cell[i_cell];
                    
                    tmp_x[(3 * i_list)] =       C.x[i_cell][0];
                    tmp_x[(3 * i_list + 1)] =   C.x[i_cell][1];
                    tmp_x[(3 * i_list + 2)] =   C.x[i_cell][2];
                    
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
                        
                        tmp_x[(3 * i_list)] =       list_of_real[(3 * i_cell)];
                        tmp_x[(3 * i_list + 1)] =   list_of_real[(3 * i_cell + 1)];
                        tmp_x[(3 * i_list + 2)] =   list_of_real[(3 * i_cell + 2)];
                                    
                        i_list++;
                    }
                    
                    free(list_of_int);
                    
                    free(list_of_real);
                }
            }
        }
        
        /* E.2 analyze MPI_cells.host_of_ext_cells, set MPI_cells.host2ext_index */
        {
            if(myid == 0){
            
                for(i_host = 1; i_host < (node_last+1); i_host++){
                    
                    MPI_cells.number_of_host_cells_per_node[i_host] = 0;
                }
                
                for(i_cell = 0; i_cell < MPI_cells.number_of_ext_cells; i_cell++){
                    
                    i_host = MPI_cells.host_of_ext_cells[i_cell];
                    
                    MPI_cells.number_of_host_cells_per_node[i_host]++;
                }
                
                int host_start_index[(node_last+1)], host_index[(node_last+1)];
                    
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
        
        /* E.3 set MPI_cells.hosting_cell_index */
        {
            if(myid == 0){
                
                list_of_real = (double*)malloc( 3 * MPI_cells.number_of_host_cells_per_node[0] * sizeof(double));
                
                i_list = 0;
                
                for(i_cell = 0; i_cell < MPI_cells.number_of_ext_cells; i_cell++){
                    
                    i_host = MPI_cells.host_of_ext_cells[i_cell];
                    
                    if(i_host == 0){
                        
                        list_of_real[(3 * i_list)] =        tmp_x[(3 * i_cell)]; 
                        list_of_real[(3 * i_list + 1)] =    tmp_x[(3 * i_cell + 1)]; 
                        list_of_real[(3 * i_list + 2)] =    tmp_x[(3 * i_cell + 2)]; 
                        
                        i_list++;
                    }
                }
                
                MPI_cells.hosting_cell_index = (int*)malloc(MPI_cells.number_of_host_cells_per_node[0] * sizeof(int));
                
                for(i_list = 0; i_list < MPI_cells.number_of_host_cells_per_node[0]; i_list++){
                    
                    loop_cells{
                        
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
                            
                            list_of_real[(3 * i_list)] =        tmp_x[(3 * i_cell)]; 
                            list_of_real[(3 * i_list + 1)] =    tmp_x[(3 * i_cell + 1)]; 
                            list_of_real[(3 * i_list + 2)] =    tmp_x[(3 * i_cell + 2)]; 
                            
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
                    
                    loop_cells{
                        
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
#endif      
    }   

    /* F. allocate and set MPI_faces */
    {
#if RP_NODE

        MPI_faces.principal_face = (int*)malloc(Face_Dict.number_of_faces * sizeof(int));
        
        rCFD_default_MPI_Faces(&Solver_Dict, &MPI_faces);
#endif      
    }

#endif  
    /* G. Message & Transcript */
    {
        if((myid == 0) && (Solver_Dict.verbal == 1)){
                
            FILE    *f_out = fopen("./Run.trn", "a");
            
            if(f_out){
                
                fprintf(f_out,"\n\nrCFD_prepare_C2Cs_for_MPI");
            
                fprintf(f_out,"\n\n   Prepared C2Cs in %f seconds", 
                    
                    (double)(clock() - Solver.clock)/(double)CLOCKS_PER_SEC);                                           

                fclose(f_out);
            }       
        }
        
        Message0("\n...rCFD_prepare_C2Cs_for_MPI -> Prepared C2Cs in %f seconds", 
        
            (double)(clock() - Solver.clock)/(double)CLOCKS_PER_SEC);
    }
#endif
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
    
    int         i_run, i_island, i_phase, i_state, i_state2;
    
    double      rand_per_phase_loop;
    
#if RP_NODE 
    int         i_frame, i_data, i_C2C, i_node, i_node2, i_swap, i_list, i_cell, i_face, i_drift, i_dim;

    int         loop_offset0, loop_offset1, loop_offset_in0, loop_offset_in1;
    int         size, MPI_size, MPI_index, index_in_total_list;
    int         number_of_fill_loops, number_of_unhit_cells;
    int         c, c0, c1, n0, n1;
    int         max_swap_loops;
    int         i_frame_c0, i_frame_c1;
    int         i_tmp;
    
    double      w0, data0, vol_flip;
    double      *list_of_real = NULL;
    double      drift_volume, drift_mass, local_mass_c0, local_mass_c1, local_mass;
    double      flux_in, flux_out, data_in_mean, data_out_mean, flux_mean;

    FILE        *f_out = NULL;
#else

    double      rand_real;
#endif

#endif

    Solver.clock = clock();
    
    loop_runs{
        
        /* ST: select state in RAM for all phases */
        {   
            i_state = 0, i_state2 = 0;      
        }   
        
        /* N+1: Get next frames[islands] for all phases */
        {
            
            if(Solver_Dict.recurrence_process_on){
            
                if(Rec.frame_in_sequence < Rec.sequence_length){ 
                
                    Rec.frame_in_sequence++;
                }
                else{
                    
                    /* new rCFD_seq */
                    
                    Rec.frame_in_sequence = 0; 
        
#if RP_HOST
                    rand_real = (double)rand()/(double)RAND_MAX; 
                    
                    Rec.sequence_length = Rec_Dict.min_seq_length + (int)(rand_real*(double)(Rec_Dict.max_seq_length - Rec_Dict.min_seq_length));
#endif    
                    host_to_node_int_1(Rec.sequence_length);

                    loop_islands{
                        
                        Rec.global_frame[i_island] = Rec.jumps[i_state][i_state2][i_island][Rec.global_frame[i_island]];
                    }
                }
          
                /* set next frame */
                loop_islands{
            
                    if(Rec.global_frame[i_island] < (Solver_Dict.number_of_frames - 1)){
                        
                        Rec.global_frame[i_island]++;
                    }
                    else{ 
                        
                        Rec.global_frame[i_island] = Rec.jumps[i_state][i_state2][i_island][(Solver_Dict.number_of_frames-1)];
                    }
                }
            }
            else{
                
                loop_islands{
                    
                    Rec.global_frame[i_island] = 0;
                }
            }       
        }
                    
        loop_phases{
            
            /* I: Initialization (weights, data_shift, data_swap, mass_drift) */
            {
#if RP_HOST
                rand_per_phase_loop = (double)rand()/(double)RAND_MAX;  /* [0..1] */
#endif
                host_to_node_real_1(rand_per_phase_loop);
#if RP_NODE

                loop_cells{
                    
                    C.weight_after_shift[i_cell] =  0.0;
                    C.weight_after_swap[i_cell] =   0.0;
                    
                    loop_data{
                        
                        C.data_shift[i_phase][i_cell][i_data] = 0.0;
                        C.data_swap[i_phase][i_cell][i_data] =  0.0;
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
#endif          
            }
    
            /* AD1: Access data before shift */     
            {
#if RP_NODE
                rCFD_user_access_data_before_shift(Balance, Phase_Dict, &Cell_Dict, &C, &Rec, i_phase);
#endif                  
            }
            
            /* C: Convection, cell-to-cell shifts */
            {
#if RP_NODE         
                if(Solver_Dict.data_convection_on){
                    
                    if(rand_per_phase_loop <= Phase_Dict[i_phase].shift_probability){
                    
                        /* C1,2 C2C shifts (local, MPI) */
                        loop_islands{
                            
                            i_frame = Rec.global_frame[i_island];
                            
                            /* C1: local cell shifts */
                            {
                                loop_offset0    = C2Cs[current_pattern].island_offsets[i_island];
                                loop_offset1    = C2Cs[current_pattern].island_offsets[(i_island + 1)];
                                
                                /*Message0("\n\n loop_offset0 %d, loop_offset1 %d\\", loop_offset0, loop_offset1);*/

                                i_tmp = 0;
                                
                                for(i_C2C = loop_offset0; i_C2C < loop_offset1; i_C2C++){
                                    
                                    c0 = C2Cs[current_pattern].shifts[i_C2C].c0;
                                    c1 = C2Cs[current_pattern].shifts[i_C2C].c1;
                                    w0 = C2Cs[current_pattern].shifts[i_C2C].w0;
                                
                                    if(w0 > 0.0){
                             
                                        i_tmp++;
                                        
                                        loop_data{
                                 
                                            data0 = C.data[i_phase][c0][i_data];

                                            C.data_shift[i_phase][c1][i_data] =
                                            
                                                (w0 * data0 + C.weight_after_shift[c1] * C.data_shift[i_phase][c1][i_data])/
                                                
                                                (w0 +  C.weight_after_shift[c1]);                                       
                                        }
                             
                                        C.weight_after_shift[c1] += w0;
                                    }
                                }
                                
                                /*Message0("\n\n i_tmp = %d\n\n", i_tmp);*/
                            }
                        
                            /* C2: MPI shifts */
                       
                            /* C2.a: fill local C2Cs.shifts_out and send it to Node-0 */
                            {
                                if (myid > 0){

                                    size = C2Cs[current_pattern].number_of_shifts_out;
                                    
                                    PRF_CSEND_INT(node_zero, &size, 1, myid);
                              
                                    if(size > 0){
                              
                                        MPI_size = size * Phase_Dict[i_phase].number_of_data;
                                        
                                        list_of_real = (double*)malloc(MPI_size * sizeof(double));
                                        
                                        loop_C2Cs_size{

                                            c0 = C2Cs[current_pattern].shifts_out[i_C2C].c0;
                                            
                                            loop_data{
                                                
                                                i_list =  Phase_Dict[i_phase].number_of_data * i_C2C + i_data;
                                   
                                                list_of_real[i_list] = C.data[i_phase][c0][i_data];
                                                
                                                if(Balance_Dict[i_phase][i_data].type == per_node_balancing){
                                                
                                                    w0 = C2Cs[current_pattern].shifts_out[i_C2C].w0;
                                                    n1 = C2Cs[current_pattern].shifts_out[i_C2C].node1;
                                                    
                                                    if((n1<0) ||(n1>node_last)){
                                                        
                                                        Message("\n ERROR");
                                                    }
                                                    if((w0 < 0) || (w0 > 1e5)){
                                                        
                                                        Message("\n ERROR w0");
                                                    }
                                                        
                                                    Balance[i_phase][i_data].node2node_flux[myid][n1] += w0;
                                                    Balance[i_phase][i_data].node2node_data_flux[myid][n1] += w0 * C.data[i_phase][c0][i_data];
                                                }
                                            }
                                        }
                                        
                                        PRF_CSEND_REAL(node_zero, list_of_real, MPI_size, myid);
                                  
                                        free(list_of_real);
                                    }
                                }
                            
                                if (myid == 0){
                            
                                    size = C2Cs[current_pattern].number_of_shifts_out;

                                    loop_C2Cs_size{
                                        
                                        c0 = C2Cs[current_pattern].shifts_out[i_C2C].c0;
                                        
                                        MPI_index = C2Cs[current_pattern].in2out[node_zero][i_C2C];
                                        
                                        loop_data{
                                   
                                            C2Cs_MPI.shifts_in[MPI_index].data[i_data] = C.data[i_phase][c0][i_data];
                                            
                                            if(Balance_Dict[i_phase][i_data].type == per_node_balancing){
                                            
                                                w0 = C2Cs[current_pattern].shifts_out[i_C2C].w0;
                                                n1 = C2Cs[current_pattern].shifts_out[i_C2C].node1;
                                                
                                                Balance[i_phase][i_data].node2node_flux[myid][n1] += w0;
                                                Balance[i_phase][i_data].node2node_data_flux[myid][n1] += w0 * C.data[i_phase][c0][i_data];

                                                if((n1<0) ||(n1>node_last)){
                                                    
                                                    Message("\n ERROR");
                                                }
                                            }
                                            
                                        }
                                    } 
                                    
                                    /* fill rCFD_TRACER_COMM by Node>0 rCFD_TRACER_OUTGOING data */
                                    for(i_node = 1; i_node < (node_last+1); i_node++){
                                        
                                        PRF_CRECV_INT(i_node, &size, 1, i_node);
                                        
                                        if(size > 0){

                                            MPI_size = size * Phase_Dict[i_phase].number_of_data;   
                                            
                                            list_of_real = (real*)malloc(MPI_size * sizeof(double));
                                            
                                            PRF_CRECV_REAL(i_node, list_of_real, MPI_size, i_node);
                                        
                                            loop_C2Cs_size{
                                           
                                                MPI_index = C2Cs[current_pattern].in2out[i_node][i_C2C];
                                                
                                                loop_data{
                                                    
                                                    i_list = Phase_Dict[i_phase].number_of_data * i_C2C + i_data;

                                                    C2Cs_MPI.shifts_in[MPI_index].data[i_data] = list_of_real[i_list]; 
                                                }
                                            }

                                            free(list_of_real);                 
                                        }
                                    }
                                }
                            }
                        
                            /* C2.b: Node-0 sends C2Cs.shifts_in to other Nodes */
                            {
                                if (myid == 0){
                                    
                                    size = C2Cs[current_pattern].number_of_shifts_in;
                                    
                                    index_in_total_list = size;             

                                    for(i_node = 1; i_node < (node_last+1); i_node++){
                                    
                                        size = C2Cs[current_pattern].number_of_shifts_to_node_zero[i_node]; 
                                        
                                        PRF_CSEND_INT(i_node, &size, 1, node_zero);
                                        
                                        if(size > 0){
                                            
                                            MPI_size = size * Phase_Dict[i_phase].number_of_data;
                                            
                                            list_of_real=(real*)malloc(MPI_size * sizeof(double));
                                            
                                            loop_C2Cs_size{
                                            
                                                loop_data{
                                                    
                                                    i_list = Phase_Dict[i_phase].number_of_data * i_C2C + i_data;
                                            
                                                    list_of_real[i_list] = C2Cs_MPI.shifts_in[(index_in_total_list + i_C2C)].data[i_data];
                                                }
                                            }

                                            PRF_CSEND_REAL(i_node, list_of_real, MPI_size, node_zero);
                                            
                                            free(list_of_real);
                                            
                                            index_in_total_list += size;
                                        }
                                    }               
                                }

                                if (myid > 0){
                            
                                    PRF_CRECV_INT(node_zero, &size, 1, node_zero);
                                    
                                    if(size != C2Cs[current_pattern].number_of_shifts_in){
                                        
                                        Message("\n myid %d, size %d, number_of_shifts_in %d", myid, size, C2Cs[current_pattern].number_of_shifts_in);
                                    }
                              
                                    if(size > 0){
                                  
                                        MPI_size = size * Phase_Dict[i_phase].number_of_data;
                                        
                                        list_of_real=(real*)malloc(MPI_size * sizeof(double));
                                        
                                        PRF_CRECV_REAL(node_zero, list_of_real, MPI_size, node_zero); 
                                    }
                                }
                            }

                            /* C2.c: redistribute on Node (data_shift, w_after_shift) */
                            {

                                loop_offset_in0 = C2Cs[current_pattern].island_offsets_in[i_island];
                                loop_offset_in1 = C2Cs[current_pattern].island_offsets_in[(i_island + 1)];

                                for(i_C2C = loop_offset_in0; i_C2C < loop_offset_in1; i_C2C++){
                                                                
                                    c0 = C2Cs[current_pattern].shifts_in[i_C2C].c0;
                                    c1 = C2Cs[current_pattern].shifts_in[i_C2C].c1;
                                    w0 = C2Cs[current_pattern].shifts_in[i_C2C].w0;
                                    
                                    if(w0 > 0.0){
                                    
                                        loop_data{
                             
                                            if (myid == 0){
                                                
                                                data0 = C2Cs_MPI.shifts_in[i_C2C].data[i_data];
                                            }
                                            else{
                                                
                                                i_list = i_C2C * Phase_Dict[i_phase].number_of_data + i_data;
                                                
                                                data0 = list_of_real[i_list];
                                            }
                                                            
                                            C.data_shift[i_phase][c1][i_data] =  
                                            
                                                (w0 * data0 + C.weight_after_shift[c1] * C.data_shift[i_phase][c1][i_data])/
                                                
                                                (w0 + C.weight_after_shift[c1]);
                                                
                                                
                                            if(Balance_Dict[i_phase][i_data].type == per_node_balancing){
                                            
                                                n0 = C2Cs[current_pattern].shifts_in[i_C2C].node0;
                                                
                                                Balance[i_phase][i_data].node2node_flux[n0][myid] += w0;
                                                Balance[i_phase][i_data].node2node_data_flux[n0][myid] += w0 * data0;
                                            }
                                        }
                                
                                        C.weight_after_shift[c1] += w0; 
                                    }
                                }
                                
                                if(myid > 0){
                                    
                                    free(list_of_real);
                                }
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

                                    c0 = F.c0[i_face];
                                    c1 = F.c1[i_face];
                                
                                    if((C.weight_after_shift[c0] > 0.)&&(C.weight_after_shift[c1] == 0.0)){
                                        
                                        c = c0; c0 = c1; c1 = c;
                                    }
                                  
                                    if((C.weight_after_shift[c1] > 0.)&&(C.weight_after_shift[c0] == 0.0)){
                                  
                                        loop_data{
                                            
                                            C.data_swap[i_phase][c0][i_data] = 
                                            
                                                (C.weight_after_shift[c1] * C.data_shift[i_phase][c1][i_data] + 
                                                
                                                 C.weight_after_swap[c0] * C.data_swap[i_phase][c0][i_data]) / 
                                                
                                                (C.weight_after_shift[c1] + C.weight_after_swap[c0]);
                                        }
                                                                  
                                        C.weight_after_swap[c0] += C.weight_after_shift[c1];
                                    }                       
                                }
                        
                                number_of_unhit_cells = 0;
                                
                                loop_cells{
                               
                                    if(C.weight_after_swap[i_cell] > 0.0){
                                        
                                        loop_data{
                                            
                                            C.data_shift[i_phase][i_cell][i_data] = C.data_swap[i_phase][i_cell][i_data];
                                        }

                                        C.weight_after_shift[i_cell] = C.weight_after_swap[i_cell];
                                        
                                        C.weight_after_swap[i_cell] = 0.0;
                                    }
                                    
                                    if(C.weight_after_swap[i_cell] == 0.0){
                                        
                                        number_of_unhit_cells++;
                                    }
                                     
                                }

                                number_of_unhit_cells = PRF_GISUM1(number_of_unhit_cells);
                            }

                            /* data = data_update */
                            loop_cells{
                            
                                loop_data{
                                    
                                    C.data[i_phase][i_cell][i_data] = C.data_shift[i_phase][i_cell][i_data];
                                }   
                            }   
                        }
                    }
                }
#endif      
            }

            /* D: Diffusion, face swaps */
            {               
                /* D.1 Drifting */
                if(Solver_Dict.data_drifting_on){
                    
                    /* Existing Limitations (10/21)
                    
                        L1: drifting of concentration might lead to conc > 1. 
                    */
#if RP_NODE
                    loop_data{
                            
                        if(Data_Dict[i_phase][i_data].drifting_data){
                            
                            for(i_drift = 0; i_drift < Solver_Dict.number_of_drift_loops; i_drift++){
                                                
                                /* init C.drift_exchange */
                                loop_cells{

                                    C.drift_exchange[i_cell] = 0.0;
                                }
                                
                                /* set C.drift_exchange */
                                loop_faces{
                                    
                                    if(valid_parallel_face){
                                        
                                        c0 = F.c0[i_face];
                                        c1 = F.c1[i_face];
                                        
                                        i_frame_c0 = Rec.global_frame[C.island_id[c0]]; 
                                        i_frame_c1 = Rec.global_frame[C.island_id[c1]]; 
                                        
                                        local_mass_c0 = C.volume[c0] * C.vof[i_frame_c0][c0][i_phase] * Phase_Dict[i_phase].density;
                                        local_mass_c1 = C.volume[c1] * C.vof[i_frame_c1][c1][i_phase] * Phase_Dict[i_phase].density;
                                            
                                        /* don't drift across interfaces */
                                        if(local_mass_c0 * local_mass_c1 > 0.0){
                                            
                                            drift_volume = 0.0;
                                            
                                            loop_dim{
                                                
                                                drift_volume += Data_Dict[i_phase][i_data].drift_velocity[i_dim] * F.area[i_face][i_dim] * Phase_Dict[i_phase].time_step / 
                                            
                                                    (double)Solver_Dict.number_of_drift_loops;
                                            }
                                            
                                            /* flip cells, such that flux is from c0 to c1 */
                                            if(drift_volume < 0.0){
                                                
                                                c = c0; c0 = c1; c1 = c;
                                                
                                                i_frame = i_frame_c0; i_frame_c0 = i_frame_c1; i_frame_c1 = i_frame;
                                                
                                                local_mass = local_mass_c0; local_mass_c0 = local_mass_c1; local_mass_c1 = local_mass;
                                                
                                                drift_volume *= -1.0;
                                            }
                                                                                                                                        
                                            drift_mass = C.data[i_phase][c0][i_data] * drift_volume * C.vof[i_frame_c0][c0][i_phase] * Phase_Dict[i_phase].density;
                                                
                                            C.drift_exchange[c0] -= drift_mass; 
                                            C.drift_exchange[c1] += drift_mass;     
                                        }                                           
                                    }                               
                                }
                                
                                /* parallel exchange of drift_mass */
                                {
                                    sum_up_parallel_corona_cells(&Solver_Dict, &Cell_Dict, C.drift_exchange); 
                                }
                                
                                /* adjust data by C.drift_exchange */
                                loop_cells{

                                    i_frame = Rec.global_frame[C.island_id[i_cell]];        
                                    
                                    local_mass = C.volume[i_cell] * C.vof[i_frame][i_cell][i_phase] * Phase_Dict[i_phase].density;

                                    if(local_mass > 0.0){
                                        
                                        C.data[i_phase][i_cell][i_data] = ((C.data[i_phase][i_cell][i_data] * local_mass) + C.drift_exchange[i_cell]) / 
                                        
                                            local_mass;
                                    }
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
                    
                    /* Init data_swap */
                    loop_cells{
                        
                        loop_data{
                            
                            C.data_swap[i_phase][i_cell][i_data] = 0.0;
                        }           
                    }   

                    /* define swap masses */
                    loop_faces{

                        c0 = F.c0[i_face];
                        c1 = F.c1[i_face];
                        
                        loop_data{
                            
                            if(C.data[i_phase][c0][i_data] > C.data[i_phase][c1][i_data]){
                                
                                c = c0; c0 = c1; c1 = c;
                            }
                          
                            if(C.data[i_phase][c1][i_data] > C.data[i_phase][c0][i_data]){
                                
                                if(C.volume[c0] < C.volume[c1]){
                                    
                                    vol_flip = C.volume[c0] * Data_Dict[i_phase][i_data].physical_diff;
                                }
                                else{
                                    vol_flip = C.volume[c1] * Data_Dict[i_phase][i_data].physical_diff;
                                }
                                
                                C.data_swap[i_phase][c1][i_data] -= vol_flip * 
                                
                                    (C.data[i_phase][c1][i_data] - C.data[i_phase][c0][i_data]) / 2.;

                                C.data_swap[i_phase][c0][i_data] += vol_flip * 
                                
                                    (C.data[i_phase][c1][i_data] - C.data[i_phase][c0][i_data]) / 2.;
                            }
                        }   
                    }
        
                    /* update data by data_swap */
                    loop_cells{
                        
                        loop_data{
                  
                            C.data[i_phase][i_cell][i_data] += C.data_swap[i_phase][i_cell][i_data] / C.volume[i_cell];

                            C.data_swap[i_phase][i_cell][i_data] = 0.0;
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

                                C.data_swap[i_phase][i_cell][i_data] = 0.0;
                            }   

                            /* Artificial Diffusion (mass neutral) */
                            loop_faces{

                                c0 = F.c0[i_face];
                                c1 = F.c1[i_face];
                                    
                                if(C.data[i_phase][c0][i_data] > C.data[i_phase][c1][i_data]){
                                    
                                    c = c0; c0 = c1; c1 = c;
                                }
                              
                                if(C.data[i_phase][c1][i_data] > C.data[i_phase][c0][i_data]){
                                    
                                    if(C.volume[c0] < C.volume[c1]){
                                        
                                        vol_flip = C.volume[c0] * Data_Dict[i_phase][i_data].binarization_art_diff;
                                    }
                                    else{
                                        vol_flip = C.volume[c1] * Data_Dict[i_phase][i_data].binarization_art_diff;
                                    }
                                    
                                    C.data_swap[i_phase][c1][i_data] -= vol_flip * 
                                    
                                        (C.data[i_phase][c1][i_data] - C.data[i_phase][c0][i_data]) / 2.;

                                    C.data_swap[i_phase][c0][i_data] += vol_flip * 
                                    
                                        (C.data[i_phase][c1][i_data] - C.data[i_phase][c0][i_data]) / 2.;
                                }
                            }   
                            
                            /* Binarization (mass acting) */
                            loop_cells{

                                if(C.data[i_phase][i_cell][i_data] > 0.5){
                                        
                                    C.data[i_phase][i_cell][i_data] = 1.0;
                                }
                                else{ 
                                    
                                    C.data[i_phase][i_cell][i_data] = 0.0;
                                }   
                            }
                        }           
                    }       
#endif          
                }       
            }
            
            /* AD2: Access data after swap */
            {
#if RP_NODE                 
                rCFD_user_access_data_after_swap(Balance, Phase_Dict, &Cell_Dict, &C, i_phase); 
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

                            Balance[i_phase][i_data].mass_integral = 0.0;
                        }               
                            
                        loop_int_cells{
                        
                            i_frame = Rec.global_frame[C.island_id[i_cell]];
                            
                            loop_data{
                                
                                if(Data_Dict[i_phase][i_data].type == temperature_data){
                                    
                                    Balance[i_phase][i_data].mass_integral += (C.data[i_phase][i_cell][i_data] - Phase_Dict[i_phase].reference_temperature) * 
                                    
                                        C.volume[i_cell] * C.vof[i_frame][i_cell][i_phase] * Phase_Dict[i_phase].density * Phase_Dict[i_phase].heat_capacity;
                                }
                                else{
                                    
                                    Balance[i_phase][i_data].mass_integral += C.data[i_phase][i_cell][i_data] * 
                                    
                                        C.volume[i_cell] * C.vof[i_frame][i_cell][i_phase] * Phase_Dict[i_phase].density;
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
                            
                            Balance[i_phase][i_data].mass_integral_target += Balance[i_phase][i_data].mass_source;
                            
                            Balance[i_phase][i_data].mass_integral_target_global += PRF_GRSUM1(Balance[i_phase][i_data].mass_source);
                            
                            Balance[i_phase][i_data].mass_source = 0.0;
                        }
#endif                      
                    }
                    
                    /* B3: mass_integral_target += node2node fluxes */
                    {
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

                    /* B4: mass_integral, mass_integral_target */
                    {
#if RP_NODE
                        loop_data{ 
                        
                            if(Balance_Dict[i_phase][i_data].type == global_balancing){
                                
                                Balance[i_phase][i_data].mass_integral = Balance[i_phase][i_data].mass_integral_global;
                                
                                Balance[i_phase][i_data].mass_integral_target = Balance[i_phase][i_data].mass_integral_target_global;
                            }
                        }
#endif                      
                    }

                    /* B5:  data_int_error  */
                    {
#if RP_NODE                     
                        if(Solver.global_run_counter == 0){
                            
                            loop_data{       
                                
                                Balance[i_phase][i_data].mass_error =       
                                
                                    Balance[i_phase][i_data].mass_integral_target - Balance[i_phase][i_data].mass_integral;
                                
                                Balance[i_phase][i_data].mass_error_prev = Balance[i_phase][i_data].mass_error;                             
                            }
                        }               
                        else{
                            
                            loop_data{
                                
                                Balance[i_phase][i_data].mass_error_prev = Balance[i_phase][i_data].mass_error;                             

                                Balance[i_phase][i_data].mass_error =       
                                
                                    Balance[i_phase][i_data].mass_integral_target - Balance[i_phase][i_data].mass_integral;
                            } 
                        }   
#endif                      
                    }

                    /* B6. Balance_face_swap */
                    {
#if RP_NODE                     
                        if((i_run % Solver_Dict.balance_correction_update) == 0){
                        
                            loop_data{

                                /* smaller balance correction */
                                if((Balance[i_phase][i_data].mass_error_prev * Balance[i_phase][i_data].mass_error) <= 0.){ 
                             
                                    Balance[i_phase][i_data].face_swap /= 2.;
                                    
                                    if(Balance[i_phase][i_data].face_swap < Solver_Dict.face_swap_min){ 
                                        
                                        Balance[i_phase][i_data].face_swap = Solver_Dict.face_swap_min;
                                    }
                                    
                                    if(Balance[i_phase][i_data].face_swap_loops > 1){
                                        
                                        Balance[i_phase][i_data].face_swap_loops /= 2.;
                                    }
                                }
                                /* larger balance correction */
                                else{
                                        
                                    if(Balance[i_phase][i_data].face_swap >= Solver_Dict.face_swap_max_per_loop){
                                        
                                        Balance[i_phase][i_data].face_swap_loops *= 2.;
                                    }
                                    
                                    Balance[i_phase][i_data].face_swap *= 2.;
                                    
                                    if(Balance[i_phase][i_data].face_swap > Solver_Dict.face_swap_max_per_loop){ 
                            
                                        Balance[i_phase][i_data].face_swap = Solver_Dict.face_swap_max_per_loop;
                                    }

                                    if(Balance[i_phase][i_data].face_swap_loops > Solver_Dict.face_swap_max_loops){ 
                            
                                        Balance[i_phase][i_data].face_swap_loops = Solver_Dict.face_swap_max_loops;
                                    }
                                }   
                            }                           
                        }
#endif                      
                    }

                    /* B7. write global Balance (node-0) */
                    {
#if RP_NODE                     
                        if(myid == 0){
                            
                            loop_data{
                                
                                if( (Balance_Dict[i_phase][i_data].write_balance_to_file) &&
                                    ((Solver.global_run_counter % Balance_Dict[i_phase][i_data].write_balance_to_file_interval) == 0)){
                                    
                                    if((Solver.global_run_counter == 0) && (i_phase == 0) && (i_data == 0)){

                                        f_out = fopen(File_Dict.Balance_filename, "w");                                     
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
                }
                
                /* Balance correction by face swaps */
                if(Solver_Dict.balance_correction_on){  
#if RP_NODE 
                    if((i_run % Solver_Dict.balance_correction_update) == 0){

                        /* max_swap_loops */
                        {
                            max_swap_loops = 1.0;   
                            
                            loop_data{
                                
                                if(Balance[i_phase][i_data].face_swap_loops > max_swap_loops){
                                    
                                    max_swap_loops = Balance[i_phase][i_data].face_swap_loops;
                                }
                            }
                        }
                        
                        
                        /* loop_swap (based on max_swap_loops) */
                        loop_max_swap_loops{
                            
                            /* init data_swap */
                            loop_cells{
                                
                                loop_data{
                                    
                                    if(Balance[i_phase][i_data].face_swap_loops >= i_swap){
                                        
                                        C.data_swap[i_phase][i_cell][i_data] = 0.0;
                                    }
                                }           
                            }   
                            
                            /* loop faces and fill data_swap */ 
                            loop_faces{

                                c0 = F.c0[i_face];
                                c1 = F.c1[i_face];
                                
                                loop_data{
                                    
                                    if(Balance[i_phase][i_data].face_swap_loops >= i_swap){
                                        
                                        if(C.data[i_phase][c0][i_data] > C.data[i_phase][c1][i_data]){
                                            
                                            c = c0; c0 = c1; c1 = c;
                                        }
                                      
                                        if(C.data[i_phase][c1][i_data] > C.data[i_phase][c0][i_data]){
                                            
                                            if(C.volume[c0] < C.volume[c1]){
                                                
                                                vol_flip = C.volume[c0] * Balance[i_phase][i_data].face_swap;
                                            }
                                            else{
                                                vol_flip = C.volume[c1] * Balance[i_phase][i_data].face_swap;
                                            }
                                            
                                            if(Balance[i_phase][i_data].mass_integral > Balance[i_phase][i_data].mass_integral_target){
                                                
                                                C.data_swap[i_phase][c1][i_data] -= vol_flip * 
                                            
                                                    (C.data[i_phase][c1][i_data] - C.data[i_phase][c0][i_data]) / 2.;
                                            }
                                            else{

                                                C.data_swap[i_phase][c0][i_data] += vol_flip * 
                                                
                                                    (C.data[i_phase][c1][i_data] - C.data[i_phase][c0][i_data]) / 2.;
                                            }
                                        }
                                    }
                                }   
                            }
                            
                            /* update data by data_swap */
                            loop_cells{
                                
                                loop_data{
                                    
                                    if(Balance[i_phase][i_data].face_swap_loops >= i_swap){
                                    
                                        C.data[i_phase][i_cell][i_data] += C.data_swap[i_phase][i_cell][i_data] / C.volume[i_cell];
                                    }
                                }
                                
                            }                               
                        }
                    }           
#endif                  
                }           

            }
            
        }   /* loop_phases */
        
        Solver.global_run_counter++;    
        
        /* Post processing */
        {
#if RP_NODE                                             
            rCFD_user_post(&Solver_Dict, Phase_Dict, &Cell_Dict, &C);
#endif  
        }
        
    }   /* loop_runs */

    /* D. Message & Transcript */
    {
#if RP_NODE     
        if((myid == 0) && (Solver_Dict.verbal == 1)){
                
            FILE    *f_out = fopen("./Run.trn", "a");
            
            if(f_out){
                
                fprintf(f_out,"\n\nrCFD_run");
            
                fprintf(f_out,"\n\n   %d Runs @ global run counter %d", 
                    
                    Solver_Dict.number_of_runs, Solver.global_run_counter);

                fprintf(f_out,"\n\n   %f seconds real world time took %f seconds compute time", 
                    
                    Solver_Dict.global_time_step * (double)Solver_Dict.number_of_runs,
                    
                    (double)(clock() - Solver.clock)/(double)CLOCKS_PER_SEC);
                    
                fclose(f_out);
            }       
        }
        
        Message0("\n...rCFD_run -> %d Runs in %f seconds @ global run counter %d\n", 
        
            Solver_Dict.number_of_runs, (double)(clock() - Solver.clock)/(double)CLOCKS_PER_SEC, Solver.global_run_counter);
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
    
#if RP_NODE 
    if((myid == 0) && (Solver_Dict.verbal == 1)){
                        
        FILE    *f_out = fopen("./Run.trn", "a");
        
        if(f_out){
            
            fprintf(f_out,"\n\nrCFD_free_all");
        
            time_t      current_time = time(NULL);
            
            char        *c_time_string = ctime(&current_time); 

            fprintf(f_out,"\n\n   Run end @ %s", c_time_string);

            fclose(f_out);
        }   
    }   
#endif
}