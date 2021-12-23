#ifndef RCFD_PARALLEL
#define RCFD_PARALLEL

#include "rCFD_types.h"
#include "rCFD_macros.h"
#include "rCFD_globals.h"

/* (C)  2021 
    Stefan Pirker
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/  

#if 1 /* macros */

#define valid_parallel_face     MPI_faces.principal_face[i_face]

#endif

#if 1 /* typedef */
#if RP_NODE

    typedef struct MPI_cells_struct
    {
        int     number_of_ext_cells;

        int     *number_of_ext_cells_per_node;

        int     *number_of_host_cells_per_node;
                            
        int     *host_of_cell;              /* all nodes, list of host of all cells */

        int     *host_of_ext_cells;         /* node-0, list of all ext. cells */
        
        int     *hosting_cell_index;        /* all nodes, list of cells, which host ext. cells of other nodes */

        int     *host2ext_index;    
        
        double  *data;
        
    } MPI_cells_type;
    
    typedef struct MPI_faces_struct
    {
        int     *principal_face;
        
    } MPI_faces_type;

    typedef struct C2C_MPI_shift_struct
    {
        int         c0,node0;
        int         c1,node1;
        double      w0;

        double      *data;
        
    } C2C_MPI_shift_type;
        
    typedef struct C2C_MPI_struct
    {
        short               format;

        int                 max_number_of_MPI_shifts;
        
        C2C_MPI_shift_type  *shifts_in;
        
        int                 number_of_shifts_in; /* temporally use in C2Cs_prepare */
              
     } C2C_MPI_type;

    
#endif
#endif
    
#if 1 /* global vars */
#if RP_NODE

    static  MPI_cells_type      MPI_cells;
    
    static  MPI_faces_type      MPI_faces;
    
    static  C2C_MPI_type        C2Cs_MPI;
    
#endif
#endif  

#if 1 /* defaults */
#if RP_NODE

    void rCFD_default_MPI_Faces(void)
    {
        Domain  *d=Get_Domain(1);
        Thread  *t;
        int     i_face;
        
        thread_loop_f(t,d){if(THREAD_TYPE(t)==THREAD_F_INTERIOR){begin_f_loop(i_face,t){
                
            if(PRINCIPAL_FACE_P(i_face,t)){
                
                MPI_faces.principal_face[i_face] = 1;
            }
            else{
                
                MPI_faces.principal_face[i_face] = 0;
            }

        }end_f_loop(i_face,t)}}
    }
        
    void rCFD_default_MPI_Cells(void)
    {
        int i_cell;
        
        Domain  *d=Get_Domain(1);
        Thread  *t;
        
        thread_loop_c(t,d){if(FLUID_CELL_THREAD_P(t)){begin_c_loop_all(i_cell,t){
                
                MPI_cells.host_of_cell[i_cell] = C_PART(i_cell,t);

        }end_c_loop_all(i_cell,t)}}     
    }
        
#endif
#endif

#if 1 /* init */
#if RP_NODE

    void init_parallel_grid(const short i_layer)
    {
        int         i_node, i_cell, i_list, i_host;
        int         size;
        int         *list_of_int = NULL;
        double      *list_of_double = NULL;
        
        /* P1: allocate MPI_cells */
        {
            MPI_cells.host_of_cell = (int*)malloc(_Cell_Dict.number_of_cells * sizeof(int));

            MPI_cells.hosting_cell_index = NULL;

            if(myid > 0){
                
                PRF_CSEND_INT(node_zero, &_Cell_Dict.number_of_ext_cells, 1, myid);
                
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
                
                MPI_cells.number_of_ext_cells_per_node[0] = _Cell_Dict.number_of_ext_cells;
                
                /* get number_of_ext_cells (sum up all nodes) */
                MPI_cells.number_of_ext_cells = _Cell_Dict.number_of_ext_cells;
                
                for(i_node = 1; i_node < (node_last+1); i_node++){
                    
                    PRF_CRECV_INT(i_node, &MPI_cells.number_of_ext_cells_per_node[i_node], 1, i_node);
                    
                    MPI_cells.number_of_ext_cells += MPI_cells.number_of_ext_cells_per_node[i_node];
                }
                
                MPI_cells.host_of_ext_cells = (int*)malloc(MPI_cells.number_of_ext_cells * sizeof(int));

                MPI_cells.host2ext_index = (int*)malloc(MPI_cells.number_of_ext_cells * sizeof(int));   
                
                MPI_cells.data = (double*)malloc(MPI_cells.number_of_ext_cells * sizeof(double));
            }
            
            rCFD_default_MPI_Cells();    
        }
            
        /* P2. MPI_cells communication */
        {           
            double *tmp_x;
            
            /* P2.1 fill MPI_cells.host_of_ext_cells, tmp_x */
            {
                if(myid > 0){
                    
                    list_of_int = (int*)malloc( _Cell_Dict.number_of_ext_cells * sizeof(int));
                    
                    list_of_double = (double*)malloc( 3 * _Cell_Dict.number_of_ext_cells * sizeof(double));
                    
                    i_list = 0;
                    
                    loop_ext_cells{
                        
                        list_of_int[i_list] = MPI_cells.host_of_cell[i_cell];
                        
                        list_of_double[(3 * i_list)] =      C.x[i_cell][0];
                        list_of_double[(3 * i_list + 1)] =  C.x[i_cell][1];
                        list_of_double[(3 * i_list + 2)] =  C.x[i_cell][2];
                        
                        i_list++;
                    }
                    
                    PRF_CSEND_INT(node_zero, list_of_int, _Cell_Dict.number_of_ext_cells, myid);

                    PRF_CSEND_REAL(node_zero, list_of_double, (3 * _Cell_Dict.number_of_ext_cells), myid);
                    
                    free(list_of_int);
                    
                    free(list_of_double);
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
                    
                        list_of_double = (double*)malloc( 3 * MPI_cells.number_of_ext_cells_per_node[i_node] * sizeof(double));
                        
                        PRF_CRECV_INT(i_node, list_of_int, MPI_cells.number_of_ext_cells_per_node[i_node], i_node);

                        PRF_CRECV_REAL(i_node, list_of_double, (3 * MPI_cells.number_of_ext_cells_per_node[i_node]), i_node);
                        
                        for(i_cell = 0; i_cell < MPI_cells.number_of_ext_cells_per_node[i_node]; i_cell++){
                            
                            MPI_cells.host_of_ext_cells[i_list] = list_of_int[i_cell];
                            
                            tmp_x[(3 * i_list)] =       list_of_double[(3 * i_cell)];
                            tmp_x[(3 * i_list + 1)] =   list_of_double[(3 * i_cell + 1)];
                            tmp_x[(3 * i_list + 2)] =   list_of_double[(3 * i_cell + 2)];
                                        
                            i_list++;
                        }
                        
                        free(list_of_int);
                        
                        free(list_of_double);
                    }
                }
            }
            
            /* P2.2 analyze MPI_cells.host_of_ext_cells, set MPI_cells.host2ext_index */
            {
                if(myid == 0){
                
                    for(i_host = 0; i_host < (node_last+1); i_host++){
                        
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
            
            /* P2.3 set MPI_cells.hosting_cell_index */
            {
                if(myid == 0){
                    
                    list_of_double = (double*)malloc( 3 * MPI_cells.number_of_host_cells_per_node[0] * sizeof(double));
                    
                    i_list = 0;
                    
                    for(i_cell = 0; i_cell < MPI_cells.number_of_ext_cells; i_cell++){
                        
                        i_host = MPI_cells.host_of_ext_cells[i_cell];
                        
                        if(i_host == 0){
                            
                            list_of_double[(3 * i_list)] =      tmp_x[(3 * i_cell)]; 
                            list_of_double[(3 * i_list + 1)] =  tmp_x[(3 * i_cell + 1)]; 
                            list_of_double[(3 * i_list + 2)] =  tmp_x[(3 * i_cell + 2)]; 
                            
                            i_list++;
                        }
                    }
                    
                    MPI_cells.hosting_cell_index = (int*)malloc(MPI_cells.number_of_host_cells_per_node[0] * sizeof(int));
                    
                    for(i_list = 0; i_list < MPI_cells.number_of_host_cells_per_node[0]; i_list++){
                        
                        loop_cells{
                            
                            if( (C.x[i_cell][0] == list_of_double[(3 * i_list + 0)]) && 
                                (C.x[i_cell][1] == list_of_double[(3 * i_list + 1)]) &&
                                (C.x[i_cell][2] == list_of_double[(3 * i_list + 2)])){
                                   
                                MPI_cells.hosting_cell_index[i_list] = i_cell;
                            }
                        }
                    }               
                    
                    free(list_of_double);
                    
                    /* send tmp_x to Nodes */
                    for(i_node = 1; i_node < (node_last+1); i_node++){
                        
                        PRF_CSEND_INT(i_node, &MPI_cells.number_of_host_cells_per_node[i_node], 1, node_zero);
                        
                        list_of_double = (double*)malloc( 3 * MPI_cells.number_of_host_cells_per_node[i_node] * sizeof(double));
                        
                        i_list = 0;
                        
                        for(i_cell = 0; i_cell < MPI_cells.number_of_ext_cells; i_cell++){
                            
                            i_host = MPI_cells.host_of_ext_cells[i_cell];
                            
                            if(i_host == i_node){
                                
                                list_of_double[(3 * i_list)] =      tmp_x[(3 * i_cell)]; 
                                list_of_double[(3 * i_list + 1)] =  tmp_x[(3 * i_cell + 1)]; 
                                list_of_double[(3 * i_list + 2)] =  tmp_x[(3 * i_cell + 2)]; 
                                
                                i_list++;
                            }
                        }
                        
                        PRF_CSEND_REAL(i_node, list_of_double, (3 * MPI_cells.number_of_host_cells_per_node[i_node]), node_zero);
                        
                        free(list_of_double);
                    }
                
                    free(tmp_x);
                }       

                if(myid > 0){
                    
                    PRF_CRECV_INT(node_zero, &size, 1, node_zero);
                    
                    MPI_cells.hosting_cell_index = (int*)malloc(size * sizeof(int));
                    
                    list_of_double = (double*)malloc(3 * size * sizeof(double));
                    
                    PRF_CRECV_REAL(node_zero, list_of_double, 3 * size, node_zero);
                    
                    for(i_list = 0; i_list < size; i_list++){
                        
                        loop_cells{
                            
                            if( (C.x[i_cell][0] == list_of_double[(3 * i_list + 0)]) && 
                                (C.x[i_cell][1] == list_of_double[(3 * i_list + 1)]) &&
                                (C.x[i_cell][2] == list_of_double[(3 * i_list + 2)])){
                                   
                                MPI_cells.hosting_cell_index[i_list] = i_cell;
                            }
                        }
                    }                   
                    
                    free(list_of_double);
                }
            }
        }   

        /* P3. allocate and set MPI_faces */
        {
            MPI_faces.principal_face = (int*)malloc(_Face_Dict.number_of_faces * sizeof(int));
            
            rCFD_default_MPI_Faces();
        }       
    }
    
    void init_parallel_C2Cs(void)
    {

        int                 i_state, i_phase, i_frame, i_out, i_node, i_C2C;
        int                 size, max_size, index_in_total_list;
        int                 *list_of_int = NULL;
        double              *list_of_double = NULL;
        
        C2C_type            tmp_C2Cs;
        
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
                                    list_of_double=(double*)malloc(size * sizeof(double));
                                    
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
                                    
                                        list_of_double[i_C2C] = C2Cs[current_pattern].shifts_in[i_C2C].w0;
                                    }
                                    
                                    PRF_CSEND_REAL(node_zero, list_of_double, size, myid);
                                    
                                    free(list_of_int);
                                    free(list_of_double);
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
                                        list_of_double=(double*)malloc(size * sizeof(real));
                                
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
                                
                                        PRF_CRECV_REAL(i_node, list_of_double, size, i_node); 
                                        
                                        loop_C2Cs_size{ 
                                        
                                            C2Cs_MPI.shifts_in[index_in_total_list + i_C2C].w0 = list_of_double[i_C2C]; 
                                        }
                                        
                                        free(list_of_int);
                                        free(list_of_double);
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
                                        list_of_double=(double*)malloc(size * sizeof(real));
                                
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
                                            
                                            list_of_double[i_C2C] = tmp_C2Cs.shifts_out[i_C2C].w0;
                                        }
                                        
                                        PRF_CSEND_REAL(i_node, list_of_double, size, node_zero);

                                        free(list_of_int); 
                                        free(list_of_double);
                                    
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
                                    list_of_double=(double*)malloc(size * sizeof(real));
                                    
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
                            
                                    PRF_CRECV_REAL(node_zero, list_of_double, size, node_zero);

                                    loop_C2Cs_size{
                                        
                                        C2Cs[current_pattern].shifts_out[i_C2C].w0=list_of_double[i_C2C];
                                    }
                                    
                                    free(list_of_int);
                                    free(list_of_double); 
                                    
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
        }   

    }

#endif
#endif

#if 1 /* communication */
#if RP_NODE

    void sum_up_parallel_corona_cells(double *cell_data, const short i_layer)
    {
        int     i_cell, i_list, i_node;
        
        int     size;
        
        double  *list_of_double;
                
        /* MPI.1 Node-0 fills MPI_cells.data */
        {
            if(myid > 0){
            
                list_of_double = (double*)malloc(_Cell_Dict.number_of_ext_cells * sizeof(double));
                
                i_list = 0;
                
                loop_ext_cells{
                    
                    list_of_double[i_list] = cell_data[i_cell];
                    
                    i_list++;
                }
                
                PRF_CSEND_INT(node_zero, &_Cell_Dict.number_of_ext_cells, 1, myid);
                
                PRF_CSEND_REAL(node_zero, list_of_double, _Cell_Dict.number_of_ext_cells, myid);
                
                free(list_of_double);
            }
            
            if(myid == 0){
                
                i_list = 0;
                
                loop_ext_cells{
                    
                    MPI_cells.data[i_list] = cell_data[i_cell];
                    
                    i_list++;
                }
                
                for(i_node = 1; i_node < (node_last+1); i_node++){
                    
                    PRF_CRECV_INT(i_node, &size, 1, i_node);
                    
                    list_of_double = (double*)malloc(size * sizeof(double));
                    
                    PRF_CRECV_REAL(i_node, list_of_double, size, i_node);
                    
                    for(i_cell = 0; i_cell < size; i_cell++){
                        
                        MPI_cells.data[i_list] = list_of_double[i_cell];
                        
                        i_list++;
                    }
                    
                    free(list_of_double);
                }
            }
        }
        
        /* MPI.2 Node-0 sends MPI_cells.data to host-Nodes */
        {
            if(myid == 0){
                
                i_list = MPI_cells.number_of_host_cells_per_node[0];
                
                for(i_node = 1; i_node < (node_last+1); i_node++){
                    
                    PRF_CSEND_INT(i_node, &MPI_cells.number_of_host_cells_per_node[i_node], 1, node_zero);
                    
                    list_of_double = (double*)malloc(MPI_cells.number_of_host_cells_per_node[i_node] * sizeof(double));
                    
                    for(i_cell = 0; i_cell < MPI_cells.number_of_host_cells_per_node[i_node]; i_cell++){
                        
                        list_of_double[i_cell] = MPI_cells.data[MPI_cells.host2ext_index[i_list]];
                        
                        i_list++;
                    }
                    
                    PRF_CSEND_REAL(i_node, list_of_double, MPI_cells.number_of_host_cells_per_node[i_node], node_zero);
                    
                    free(list_of_double);
                }           
            }
        }
                
        /* MPI.3 Node-0 collects updated data into MPI_cells.data */
        {
            if(myid == 0){
                
                /* Node-0 fills his own MPI_cells.data */
                
                size = MPI_cells.number_of_host_cells_per_node[0];
                
                list_of_double = (double*)malloc(size * sizeof(double));
                
                for(i_cell = 0; i_cell < size; i_cell++){
                    
                    list_of_double[i_cell] = MPI_cells.data[MPI_cells.host2ext_index[i_cell]];
                }
                
                for(i_list = 0; i_list < size; i_list++){
                    
                    cell_data[MPI_cells.hosting_cell_index[i_list]] += list_of_double[i_list];
                }
                
                for(i_list = 0; i_list < size; i_list++){
                    
                    list_of_double[i_list] = cell_data[MPI_cells.hosting_cell_index[i_list]];
                }
                
                for(i_cell = 0; i_cell < size; i_cell++){
                        
                    MPI_cells.data[MPI_cells.host2ext_index[i_cell]] = list_of_double[i_cell];
                }
                
                free(list_of_double);
                
                /* Node-0 collects MPI_cells.data from other nodes */
                
                i_list = MPI_cells.number_of_host_cells_per_node[0];
                
                for(i_node = 1; i_node < (node_last+1); i_node++){
                
                    PRF_CRECV_INT(i_node, &size, 1, i_node);
                    
                    list_of_double = (double*)malloc(size * sizeof(double));
                    
                    PRF_CRECV_REAL(i_node, list_of_double, size, i_node);
                    
                    for(i_cell = 0; i_cell < size; i_cell++){
                        
                        MPI_cells.data[MPI_cells.host2ext_index[i_list]] = list_of_double[i_cell];
                        
                        i_list++;
                    }
                    
                    free(list_of_double);
                }
            }

            /* other nodes sum up data on their host cells */
            if(myid > 0){

                PRF_CRECV_INT(node_zero, &size, 1, node_zero);
                
                list_of_double = (double*)malloc(size * sizeof(double));
                
                PRF_CRECV_REAL(node_zero, list_of_double, size, node_zero);
                
                for(i_list = 0; i_list < size; i_list++){
                
                    cell_data[MPI_cells.hosting_cell_index[i_list]] += list_of_double[i_list];
                }
            
                for(i_list = 0; i_list < size; i_list++){
                
                    list_of_double[i_list] = cell_data[MPI_cells.hosting_cell_index[i_list]];
                }

                PRF_CSEND_INT(node_zero, &size, 1, myid);
                
                PRF_CSEND_REAL(node_zero, list_of_double, size, myid);

                free(list_of_double);
            }
        }
        
        /* MPI.4 Node-0 sends updated data to nodes */
        {
            if(myid == 0){
                
                i_list = MPI_cells.number_of_ext_cells_per_node[0];
                
                for(i_node = 1; i_node < (node_last+1); i_node++){
                    
                    PRF_CSEND_INT(i_node, &MPI_cells.number_of_ext_cells_per_node[i_node], 1, node_zero);
                    
                    list_of_double = (double*)malloc(MPI_cells.number_of_ext_cells_per_node[i_node] * sizeof(double));
                    
                    for(i_cell = 0; i_cell < MPI_cells.number_of_ext_cells_per_node[i_node]; i_cell++){
                        
                        list_of_double[i_cell] = MPI_cells.data[i_list];
                        
                        i_list++;
                    }
                    
                    PRF_CSEND_REAL(i_node, list_of_double, MPI_cells.number_of_ext_cells_per_node[i_node], node_zero);
                    
                    free(list_of_double);
                }
                
                i_list = 0;
                
                loop_ext_cells{
                    
                    cell_data[i_cell] = MPI_cells.data[i_list];
                    
                    i_list++;
                }
                
            }
            
            if(myid > 0){
                
                PRF_CRECV_INT(node_zero, &size, 1, node_zero);
                
                list_of_double = (double*)malloc(size * sizeof(double));
                
                PRF_CRECV_REAL(node_zero, list_of_double, size, node_zero);
                
                i_list = 0;
                
                loop_ext_cells{
                    
                    cell_data[i_cell] = list_of_double[i_list];
                    
                    i_list++;
                }
                
                free(list_of_double);
            }
        }
        
    }

    void shift_parallel_C2C_data(const int i_state, const int i_phase, const int i_frame, const int i_island)
    {
        int     i_list, i_data, i_C2C, i_node;
        
        int     size, MPI_size, MPI_index, c0, c1, n0, n1, index_in_total_list, loop_offset_in0, loop_offset_in1;
        
        double  w0, data0;

        double  *list_of_double = NULL;

        /* C2.a: fill local C2Cs.shifts_out and send it to Node-0 */
        {
            if (myid > 0){

                size = C2Cs[current_pattern].number_of_shifts_out;
                
                PRF_CSEND_INT(node_zero, &size, 1, myid);
          
                if(size > 0){
          
                    MPI_size = size * Phase_Dict[i_phase].number_of_data;
                    
                    list_of_double = (double*)malloc(MPI_size * sizeof(double));
                    
                    loop_C2Cs_size{

                        c0 = C2Cs[current_pattern].shifts_out[i_C2C].c0;
                        
                        loop_data{
                            
                            i_list =  Phase_Dict[i_phase].number_of_data * i_C2C + i_data;
               
                            list_of_double[i_list] = C.data[i_phase][c0][i_data];
                            
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
                    
                    PRF_CSEND_REAL(node_zero, list_of_double, MPI_size, myid);
              
                    free(list_of_double);
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
                        
                        list_of_double = (double*)malloc(MPI_size * sizeof(double));
                        
                        PRF_CRECV_REAL(i_node, list_of_double, MPI_size, i_node);
                    
                        loop_C2Cs_size{
                       
                            MPI_index = C2Cs[current_pattern].in2out[i_node][i_C2C];
                            
                            for(i_data = 0; i_data < Phase_Dict[i_phase].number_of_data; i_data++){
                                
                                i_list = Phase_Dict[i_phase].number_of_data * i_C2C + i_data;

                                C2Cs_MPI.shifts_in[MPI_index].data[i_data] = list_of_double[i_list]; 
                            }
                        }

                        free(list_of_double);                 
                    }
                }
            }
        }
        
        /* C2.b: Node-0 sends C2Cs.shifts_in to other Nodes; all nodes set C.data_shift */
        {
            if (myid == 0){
                
                size = C2Cs[current_pattern].number_of_shifts_in;
                
                index_in_total_list = size;             

                for(i_node = 1; i_node < (node_last+1); i_node++){
                
                    size = C2Cs[current_pattern].number_of_shifts_to_node_zero[i_node]; 
                    
                    PRF_CSEND_INT(i_node, &size, 1, node_zero);
                    
                    if(size > 0){
                        
                        MPI_size = size * Phase_Dict[i_phase].number_of_data;
                        
                        list_of_double = (double*)malloc(MPI_size * sizeof(double));
                        
                        loop_C2Cs_size{
                        
                            loop_data{
                                
                                i_list = Phase_Dict[i_phase].number_of_data * i_C2C + i_data;
                        
                                list_of_double[i_list] = C2Cs_MPI.shifts_in[(index_in_total_list + i_C2C)].data[i_data];
                            }
                        }

                        PRF_CSEND_REAL(i_node, list_of_double, MPI_size, node_zero);
                        
                        index_in_total_list += size;
                        
                        free(list_of_double);
                    }
                    
                    
                } 
            
                loop_offset_in0 = C2Cs[current_pattern].island_offsets_in[i_island];
                loop_offset_in1 = C2Cs[current_pattern].island_offsets_in[(i_island + 1)];

                for(i_C2C = loop_offset_in0; i_C2C < loop_offset_in1; i_C2C++){
                                                
                    c0 = C2Cs[current_pattern].shifts_in[i_C2C].c0;
                    c1 = C2Cs[current_pattern].shifts_in[i_C2C].c1;
                    w0 = C2Cs[current_pattern].shifts_in[i_C2C].w0;
                    
                    if(w0 > 0.0){
                    
                        loop_data{
             
                            data0 = C2Cs_MPI.shifts_in[i_C2C].data[i_data];
                                            
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
            }

            if (myid > 0){
        
                PRF_CRECV_INT(node_zero, &size, 1, node_zero);
          
                if(size > 0){
              
                    MPI_size = size * Phase_Dict[i_phase].number_of_data;
                    
                    list_of_double=(double*)malloc(MPI_size * sizeof(double));
                    
                    PRF_CRECV_REAL(node_zero, list_of_double, MPI_size, node_zero); 
                    
                    loop_offset_in0 = C2Cs[current_pattern].island_offsets_in[i_island];
                    loop_offset_in1 = C2Cs[current_pattern].island_offsets_in[(i_island + 1)];

                    for(i_C2C = loop_offset_in0; i_C2C < loop_offset_in1; i_C2C++){
                                                    
                        c0 = C2Cs[current_pattern].shifts_in[i_C2C].c0;
                        c1 = C2Cs[current_pattern].shifts_in[i_C2C].c1;
                        w0 = C2Cs[current_pattern].shifts_in[i_C2C].w0;
                        
                        if((w0 > 0.0)){
                        
                            loop_data{
                 
                                i_list = i_C2C * Phase_Dict[i_phase].number_of_data + i_data;
                                
                                data0 = list_of_double[i_list];
                                                                    
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

                    free(list_of_double);
               }
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