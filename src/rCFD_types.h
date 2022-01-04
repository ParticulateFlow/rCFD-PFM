#ifndef RCFD_TYPES
#define RCFD_TYPES

#include "udf.h"

/* (C)  2021
    Stefan Pirker
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/

    /* A. Global Dicts */

    typedef struct Solver_Dict_struct
    {
        short   version_year, version_month;

        short   verbal;

        char    run_transcript_filename[80];

        /* main characteristics */

        short   number_of_frames;
        short   number_of_states;
        short   number_of_phases;
        short   number_of_islands;
        short   number_of_layers;
        short   number_of_runs;

        /* sub-model flags */

        short   recurrence_process_on;
        short   data_convection_on;
        short   face_diffusion_on;
        short   data_binarization_on;
        short   data_drifting_on;
        short   balance_correction_on;
        short   on_the_fly_post_on;

        /* sub-model params (order of occurrence) */

        int     analyse_CFD_count;

        short   max_number_of_cells_per_time_step;

        int     time_steps_per_monitoring_interval;
        double  start_time_for_monitoring;

        short   C2C_loading_reduction;

        double  global_time_step;

        short   max_fill_loops;

        double  face_swap_max_per_loop;
        double  face_swap_min;
        short   face_swap_max_loops;

        short   number_of_drift_loops;

        short   balance_correction_update;
        short   control_conc_sum_on;

    } Solver_Dict_type;

    typedef struct File_Dict_struct
    {
        char    *tracer_start_position_filename;

        char    *C2C_filename;

        char    *Norm_filename;

        char    *vof_filename;

        char    *Jump_filename;

        char    *Matrix_filename;

        char    *Balance_filename;

        char    *Dict_filename;

    } File_Dict_type;

    typedef struct Phase_Dict_struct
    {
        int         number_of_data;

        double      time_step;

        double      shift_probability;

        double      density;

        double      heat_capacity;

        double      reference_temperature;

        int         number_of_user_vars;

        double      *user;

    } Phase_Dict_type;

    typedef struct Tracer_Dict_struct
    {
        int     number_of_Tracers_per_cell;

        short   region_of_interest_exists;

        double  ROI_x_min;
        double  ROI_x_max;
        double  ROI_y_min;
        double  ROI_y_max;
        double  ROI_z_min;
        double  ROI_z_max;

        int     coarse_graining;        /* 1: take every Tracer, 2: take every second Tracer, .... */

        short   *random_walk;           /* [i_phase] */

        short   C2C_format;

    } Tracer_Dict_type;

    enum{ /* norm_types */
        standard,
        number_of_norm_type_names
    };

    typedef struct Norm_Dict_struct
    {
        short   format;

        int     coarse_graining;        /* if ((i_cell % coarse_graining) == 0) */
        
        short   larger_than_mean_norm_only;

    } Norm_Dict_type;

    enum{ /* rec methods */
        quarter_jumps,
        number_of_rec_methods
    };

    typedef struct Rec_Dict_struct
    {
        short   method;

        short   min_seq_length;
        short   max_seq_length;

    } Rec_Dict_type;

    enum{ /* data_types */
        generic_data,
        concentration_data,
        binary_data,
        temperature_data,
        number_of_data_type_names
    };

    typedef struct Data_Dict_struct
    {
        short       type;

        double      physical_diff;

        double      binarization_art_diff;

        short       drifting_data;

        double*     drift_velocity;

        int         number_of_user_vars;

        double*     user;

    } Data_Dict_type;

    enum{ /* balancing_types */
        no_balancing,
        global_balancing,
        per_node_balancing,
        number_of_balancing_modes
    };

    typedef struct Balance_Dict_struct
    {
        short       type;
        
        short       max_correction_loops;
        
        double      accuracy_level;

        short       write_balance_to_file;

        int         write_balance_to_file_interval;

    }Balance_Dict_type;

    typedef struct Cell_Dict_struct
    {
        int     number_of_cells;

        int     number_of_int_cells;

        int     number_of_ext_cells;

        int     number_of_user_vars;

    } Cell_Dict_type;

    typedef struct Face_Dict_struct
    {
        int     number_of_int_faces;

        int     number_of_ext_faces;

        int     number_of_faces;

        } Face_Dict_type;

    typedef struct Topo_Dict_struct
    {
        Cell_Dict_type  *Cell_Dict;

        Face_Dict_type  *Face_Dict;

    } Topo_Dict_type;

    /* B. Global Vars */

    typedef struct Solver_struct
    {
        int         current_state;

        int         global_run_counter;
        
        double      global_time;
        
        double      *timestep_width_per_layer;

        clock_t     clock;

    }Solver_type;

    typedef struct Cell_struct
    {
        double      **x;
        double      *volume;

        double      **average_velocity; /* [i_phase][i_cell] */
        double      **crossing_time;

        short       *hit_by_other_cell;
        short       *island_id;

        double      *weight_after_shift;
        double      *weight_after_swap;

        double      ***vof;             /* [i_frame][i_phase][i_cell] */

        double      ***data;            /* [i_phase][i_cell][i_data] */
        double      ***data_shift;
        double      ***data_swap;

        double      *drift_exchange;    /* [i_cell] */

        double      **user;             /* [i_cell][i_user] */
        double      ***rec_user;        /* [i_frame][i_cell][i_rec_user] */

    } Cell_type;

    typedef struct Face_struct
    {
        int     *c0, *c1;

        double  **area;

    } Face_type;

    enum{   /* C2C shift formats */

        c0_n0_c1_n1_w0,

        number_of_C2C_formats
    };

    typedef struct C2C_shift_struct
    {
      int           c0,node0;
      int           c1,node1;
      double        w0;

    } C2C_shift_type;

    typedef struct C2C_struct
    {
        short               format;

        int                 number_of_shifts;
        int                 number_of_shifts_in;
        int                 number_of_shifts_out;

        C2C_shift_type      *shifts;
        C2C_shift_type      *shifts_in;
        C2C_shift_type      *shifts_out;

        int                 *island_offsets;
        int                 *island_offsets_in;

        /* MPI variables, only used by Node-0 */

        int                 number_of_MPI_shifts;
        int                 *number_of_shifts_to_node_zero;
        int                 *number_of_shifts_from_node_zero;
        int                 **in2out;

     } C2C_type;

    typedef struct Balance_struct
    {
        double      mass_integral, mass_integral_target;

        double      mass_integral_global, mass_integral_target_global;

        double      mass_source, mass_source_global;

        double      **node2node_flux, **node2node_data_flux;

        double      mass_error, mass_error_global;

    }Balance_type;

    typedef struct Tracer_struct
    {
        short   allocated;
        short   initialized;
        short   ready2write;
        short   monitoring_started;

        int     *monitor_counter;       /* Tracers per Frame [i_phase] */

        int     frame_counter;

        int     *number_of_shifts;      /* [i_phase] */

        C2C_shift_type  **shifts;

    } Tracer_type;

    typedef struct Norm_struct
    {
        int     number_of_norms;

        int     frame_counter;

        double  *norm;

    } Norm_type;

    typedef struct Rec_struct
    {
        int     *global_frame;              /* [i_island] */
        int     frame_in_sequence;
        int     sequence_length;

        int     ****jumps;                  /* [state, state2, island, frame] */

    } Rec_type;

#endif