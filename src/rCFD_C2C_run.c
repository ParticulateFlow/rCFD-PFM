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

	/* D. set up parallel C2C communication */
	{
#if RP_NODE		
		init_parallel_C2Cs(&Solver_Dict, Phase_Dict, C2Cs);
#endif
	}
	
    /* E. Message & Transcript */
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
    int         i_frame, i_data, i_C2C, i_node, i_node2, i_swap, i_cell, i_face, i_drift, i_dim;

    int         loop_offset0, loop_offset1;
    int         number_of_fill_loops, number_of_unhit_cells;
    int         c, c0, c1;
    int         max_swap_loops;
    int         i_frame_c0, i_frame_c1;
    int         i_tmp;
    
    double      w0, data0, vol_flip;
	double		drift_volume, drift_mass, local_mass_c0, local_mass_c1, local_mass;
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
                    
                        /* C1,2 C2C shifts (local, cross-partitions) */
                        loop_islands{
                            
                            i_frame = Rec.global_frame[i_island];
                            
                            /* C1: local shifts */
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
                        
                            /* C2: cross-partition shifts */
							{
								shift_parallel_C2C_data(&Solver_Dict, Phase_Dict, Balance_Dict, &C2Cs[current_pattern], &C, Balance, i_phase, i_island);
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
                rCFD_user_access_data_after_swap(Balance, Phase_Dict, Data_Dict, &Cell_Dict, &C, i_phase); 
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