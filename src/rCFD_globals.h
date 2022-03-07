#ifndef RCFD_GLOBALS
#define RCFD_GLOBALS

#include "rCFD_types.h"

/* (C)  2021
    Stefan Pirker
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/

#if 1 /* global dicts & vars */

    static  Solver_Dict_type    Solver_Dict;

    static  File_Dict_type      File_Dict;

    static  Rec_Dict_type       Rec_Dict;

    static  Solver_type         Solver;

    static  Rec_type            Rec;

    static  Topo_Dict_type      Topo_Dict;


#if RP_NODE
    static  Tracer_Dict_type    Tracer_Dict;

    static  Norm_Dict_type      Norm_Dict;

    static  Phase_Dict_type     *Phase_Dict = NULL;

    static  Data_Dict_type      **Data_Dict = NULL;

    static  Balance_Dict_type   **Balance_Dict = NULL;

    static  Topo_type           Topo;

    static  Tracer_type         Tracer;

    static  Norm_type           Norms;

    static  C2C_type            ***C2Cs = NULL;             /* [states, phases, frames] */

    static  Balance_type        **Balance = NULL;
#endif

#endif


#endif
