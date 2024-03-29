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

    enum{ /* solver modes */
        preparation_mode,
        run_mode
    };

    typedef struct Solver_Dict_struct
    {
        short   version_year, version_month;

        short   verbose;

        char    run_transcript_filename[80];

        /* main characteristics */

        short   mode;

        short   number_of_frames;
        short   number_of_states;
        short   number_of_phases;
        short   number_of_islands;
        short   number_of_layers;
        short   number_of_runs;

        /* sub-model flags */

        short   recurrence_process_on;
        short   data_convection_on;
        short   face_swap_diffusion_on;
        short   data_binarization_on;
        short   data_drifting_on;
        short   balance_correction_on;
        short   on_the_fly_post_on;

        /* sub-model params (order of occurrence) */

        int     analyse_CFD_count;

        short   max_number_of_cells_per_time_step;

        int     time_steps_per_monitoring_interval;
        double  start_time_for_monitoring;

        int     min_number_of_cells_per_layer;
        short   max_number_of_faces_per_cell;
        short   max_number_of_children;
        short   number_redist_loops_per_layer;

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
        const char    *tracer_start_position_filename;

        const char    *Prep_Transscript_filename;

        const char    *Run_Transscript_filename;

        const char    *C2C_filename;

        const char    *Norm_filename;

        const char    *vof_filename;

        const char    *Jump_filename;

        const char    *Matrix_filename;

        const char    *Balance_filename;

        const char    *Rec_Frames_filename;

        const char    *Dict_filename;

    } File_Dict_type;

    typedef struct Phase_Dict_struct
    {
        int         number_of_data;

        int         number_of_concentration_data;

        double      time_step;

        double      shift_probability;

        double      density;

        double      heat_capacity;

        double      reference_temperature;

        double      vof_max;

        short       hindered_drift_on;

        int         number_of_user_vars;

        double      *user;

    } Phase_Dict_type;

    enum{ /* tracer guiding format */
        guide_by_force_format,
        guide_by_value_format
    };

    typedef struct Tracer_Dict_struct
    {
        short   format;

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

        double  *random_walk_velocity_ratio;

        short   C2C_format;

        short   avoid_information_lock_cells_on;

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

    enum{ /* rec format */
        quarter_jumps_format,
        off_diagonal_band_format,
        replay_format,
        number_of_rec_formats
    };

    typedef struct Rec_Dict_struct
    {
        short   format;

        short   monitor_rec_frames_on;
        short   adapt_vof_stitching_on;

        short   min_seq_length;
        short   max_seq_length;

        short   off_diagonal_band_width;

        short   number_of_adapt_vof_loops;

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

        double      face_swap_diffusion;

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

        int         current_layer;

        int         global_run_counter;

        short       balance_file_opened;

        short       rec_frames_monitor_file_opened;

        double      global_time;

        double      *timestep_width;

        clock_t     clock;

    }Solver_type;

    typedef struct Cell_struct
    {
        double      **x;
        double      *volume;
        double      *grid_spacing;

        double      **average_velocity; /* [i_phase][i_cell] */
        double      **crossing_time;

        short       *hit_by_other_cell;
        short       *island_id;

        double      *weight_after_shift;
        double      *weight_after_swap;

        double      ***vof;             /* [i_frame][i_cell][i_phase] */

        double      **vof_changed;      /* [i_phase][i_cell] */

        double      ***data;            /* [i_phase][i_cell][i_data] */
        double      ***data_shift;
        double      ***data_swap;

        double      *drift_exchange;    /* [i_cell] */

        double      **user;             /* [i_cell][i_user] */
        double      ***rec_user;        /* [i_frame][i_cell][i_rec_user] */

        short       *marked;            /* [i_cell] */
        int         *parent_cell;
        short       *number_of_children;
        int         **child_index;

    } Cell_type;

    typedef struct Face_struct
    {
        int     *c0, *c1;

        double  **area;

    } Face_type;

    typedef struct Topo_struct
    {
        Cell_type   *Cell;

        Face_type   *Face;

    } Topo_type;

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

        int                 *number_of_shifts;      /* [i_layer] */
        int                 *number_of_shifts_in;
        int                 *number_of_shifts_out;

        C2C_shift_type      **shifts;               /* [i_layer][i_shift] */
        C2C_shift_type      **shifts_in;
        C2C_shift_type      **shifts_out;

        int                 **island_offsets;       /* [i_layer][i_island] */
        int                 **island_offsets_in;

        /* MPI variables, only used by Node-0 */

        int                 **number_of_shifts_to_node_zero;        /* [i_layer][i_node] */
        int                 **number_of_shifts_from_node_zero;
        int                 ***in2out;                              /* [i_layer][i_node][i_node] */

     } C2C_type;

    typedef struct tmp_C2C_struct
    {
        /* same as C2C_type, but only for one layer */
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

        int                 *number_of_shifts_to_node_zero;
        int                 *number_of_shifts_from_node_zero;
        int                 **in2out;

    } tmp_C2C_type;


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

        int     *prev_global_frame;

        int     frame_in_sequence;

        short   jumped_at_last_frame;

        int     sequence_length;

        int     ****jumps;                  /* [state, state2, island, frame] */

    } Rec_type;

    typedef struct rCFD_UDF_struct
    {
        void (*_rCFD_user_set_Solver_Dict)(void);

        void (*_rCFD_user_set_File_Dict)(void);

        void (*_rCFD_user_set_Phase_Dict)(void);

        void (*_rCFD_user_set_Tracer_Dict)(void);

        void (*_rCFD_user_set_Norm_Dict)(void);

        void (*_rCFD_user_set_Rec_Dict)(void);

        void (*_rCFD_user_set_Data_Dict)(void);

        void (*_rCFD_user_set_Balance_Dict)(void);

        void (*_rCFD_user_set_Topo_Dict)(void);

        void (*_rCFD_user_set_Cell_Dict)(const short i_layer);

        void (*_rCFD_user_set_Face_Dict)(const short i_layer);

        void (*_rCFD_user_pre_proc)(void);

        void (*_rCFD_user_init_Data)(const short i_layer);

        short (*_rCFD_user_set_layer)(const short current_layer);

        void (*_rCFD_user_access_data_before_shift)(const short i_phase, const short i_layer);

        void (*_rCFD_user_access_data_after_swap)(const short i_phase, const short i_layer);

        void (*_rCFD_user_post)(void);

        double (*_rCFD_user_set_random_walk_velocity)(void);

        void (*_rCFD_user_set_Norm)(void);

    } rCFD_UDF_type;

#endif
