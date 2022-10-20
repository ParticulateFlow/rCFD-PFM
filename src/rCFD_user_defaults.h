#ifndef RCFD_USER_DEFAULT_H
#define RCFD_USER_DEFAULT_H
/* (C)  2022
    Daniel Queteschiner
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/

    void rCFD_user_set_Solver_Dict_default(void){}

    void rCFD_user_set_File_Dict_default(void){}

    void rCFD_user_set_Phase_Dict_default(void){}

    void rCFD_user_set_Tracer_Dict_default(void){}

    void rCFD_user_set_Norm_Dict_default(void){}

    void rCFD_user_set_Rec_Dict_default(void){}

    void rCFD_user_set_Data_Dict_default(void){}

    void rCFD_user_set_Balance_Dict_default(void){}

    void rCFD_user_set_Topo_Dict_default(void){}

    void rCFD_user_set_Cell_Dict_default(const short i_layer){}

    void rCFD_user_set_Face_Dict_default(const short i_layer){}

    void rCFD_user_pre_proc_default(void){}

    void rCFD_user_init_Data_default(const short i_layer){}

    short rCFD_user_set_layer_default(const short current_layer){return 0;}

    void rCFD_user_access_data_before_shift_default(const short i_phase, const short i_layer){}

    void rCFD_user_access_data_after_swap_default(const short i_phase, const short i_layer){}

    void rCFD_user_post_default(void){}

    double rCFD_user_set_random_walk_velocity_default(void){return 0.0;}

    void rCFD_user_set_Norm_default(void){}

    void rCFD_user_register_default_functions()
    {
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_set_Solver_Dict );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_set_File_Dict );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_set_Phase_Dict );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_set_Tracer_Dict );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_set_Norm_Dict );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_set_Rec_Dict );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_set_Data_Dict );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_set_Balance_Dict );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_set_Topo_Dict );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_set_Cell_Dict );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_set_Face_Dict );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_pre_proc );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_init_Data );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_set_layer );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_access_data_before_shift );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_access_data_after_swap );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_post );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_set_random_walk_velocity );
        REGISTER_RCFD_UDF_DEFAULT( rCFD_user_set_Norm );
    }


#endif

