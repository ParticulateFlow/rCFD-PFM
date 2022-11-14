;; (C)  2022
;;  Daniel Queteschiner
;;  Particulate Flow Modelling
;;  Johannes Kepler University, Linz, Austria
;;  www.particulate-flow.at
;;
(define number_of_processors 4)
(define rCFD_src_dir  "../../src")
(define rCFD_user_src_dir  "./user_src")
(define ANSYS_Fluent_case_dir "./ansys_fluent")
(define ANSYS_Fluent_case_file  "quad_80k")
(define number_of_Timesteps_for_rCFD_analyse_CFD 5)
(define ANSYS_Fluent_number_of_CFD_episodes 101)
(define ANSYS_Fluent_simulation_timestep_width 0.005)
(define ANSYS_Fluent_number_of_timesteps_per_episode 1)
(define i 0)
;;
(ti-menu-load-string (format  #f "!cp ~a/~a.cas.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))
(ti-menu-load-string (format  #f "!cp ~a/~a.dat.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))
(ti-menu-load-string (format  #f "!cp ~a/CFD_user.c ." rCFD_user_src_dir))
(ti-menu-load-string "!rm -r libudf_cfd")
(ti-menu-load-string "!rm -r data/ip")
;;
(ti-menu-load-string "/define/user-defined/compiled-functions compile
    libudf_cfd
    yes
    CFD_user.c
    \"\"
    \"\" ")
;;
(ti-menu-load-string (format #f "/file/read-case-data ./~a.cas.h5 OK" ANSYS_Fluent_case_file))
(ti-menu-load-string "/define/user-defined/compiled-functions load \"libudf_cfd\"")
(ti-menu-load-string "/define/user-defined/execute-on-demand \"CFD_convert_csv2ip::libudf_cfd\"")
(ti-menu-load-string "!cp ./data/ip/0000.ip ./data.ip")
(ti-menu-load-string "/file/interpolate/read-data yes data.ip")
;;
;;
(ti-menu-load-string (format  #f "!cp ~a/*.c ." rCFD_src_dir))
(ti-menu-load-string (format  #f "!cp ~a/*.h ." rCFD_src_dir))
(ti-menu-load-string (format  #f "!cp ~a/*.h ." rCFD_user_src_dir))
(ti-menu-load-string "!rm -r libudf_rcfd_prep")
(ti-menu-load-string "!rm -r data/c2c")
(ti-menu-load-string "!rm -r data/tmp")
(ti-menu-load-string "!rm -r data/vof")
(ti-menu-load-string "!rm -r rec")
(ti-menu-load-string "!rm -r post")
(ti-menu-load-string "!rm rCFD_prep.trn")
;;
(ti-menu-load-string "/define/user-defined/compiled-functions compile
    libudf_rcfd_prep
    yes
    yes
    rCFD_C2C_prep.c
    \"\"
    rCFD_user.h
    rCFD_types.h
    rCFD_globals.h
    rCFD_parallel.h
    rCFD_defaults.h
    rCFD_macros.h
    rCFD_init.h
    rCFD_layer.h
    rCFD_memory.h
    rCFD_free.h
    rCFD_user_defaults.h
    \"\" "
)
;;
(ti-menu-load-string "/define/user-defined/user-defined-memory 30 q" )
(ti-menu-load-string "/define/user-defined/compiled-functions load \"libudf_rcfd_prep\"")
(ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_init_all::libudf_rcfd_prep\"")
;;
(ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_write_Tracer_Positions::libudf_rcfd_prep\"")
(ti-menu-load-string "/define/models/dpm/unsteady-tracking yes yes")
(ti-menu-load-string "/define/models/dpm/injections/create-injection
    rcfd_tracer
    no yes file no
    \"./data/tmp/tracer_start_pos.inj\"
    no no yes
    \"rcfd_init_tracers::libudf_rcfd_prep\"
    no no
    0 9999"
)
(ti-menu-load-string "/define/materials/change-create
    anthracite
    tracer
    yes
    constant 1.
    no"
)
(ti-menu-load-string "/define/models/dpm/injections/set-injection-properties
    rcfd_tracer rcfd_tracer
    no no no \"./data/tmp/tracer_start_pos.inj\"
    yes no 1 0.15
    no no no no 0 100000"
)
(ti-menu-load-string "/define/injection/injection-properties/set/pick-injections-to-set no rcfd_tracer")
(ti-menu-load-string "/define/injection/injection-properties/set/physical-models/drag-parameters rcfd_no_standard_drag::libudf_rcfd_prep")
(ti-menu-load-string "/define/models/dpm/user-defined
    \"rCFD_guide_Tracers::libudf_rcfd_prep\"
    \"none\"
    \"none\"
    \"rCFD_update_Tracers::libudf_rcfd_prep\"
    \"none\"
    15"
)
(ti-menu-load-string "define/boundary-conditions/set
    wall w_top ()
    mixture
    dpm-bc-type yes escape q"
)
(ti-menu-load-string "/define/user-defined/function-hooks/execute-at-end
    \"rCFD_write_fields::libudf_rcfd_prep\"
    \"rCFD_write_C2Cs::libudf_rcfd_prep\"
    \"\""
)
;;
(do ((i  0 (+ i  1)))
    ((= i  ANSYS_Fluent_number_of_CFD_episodes))
    (ti-menu-load-string (format #f "!cp ./data/ip/~04d.ip ./data.ip" i))
    (ti-menu-load-string "/file/interpolate/read-data yes data.ip")
    (ti-menu-load-string (format #f "/solve/set/time-step ~d" ANSYS_Fluent_simulation_timestep_width))
    (ti-menu-load-string (format #f "/solve/dual-time-iterate ~d 20" ANSYS_Fluent_number_of_timesteps_per_episode))
)
;;
(ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_write_Rec::libudf_rcfd_prep\"")
;;
(ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_free_all::libudf_rcfd_prep\"")
(ti-menu-load-string "/define/injections/delete-injection rcfd_tracer")
;;
(ti-menu-load-string "!rm  ./*.c")
(ti-menu-load-string "!rm  ./*.h")
(ti-menu-load-string "!rm  ./*.h5")
(ti-menu-load-string "!rm  ./fluent*")
(ti-menu-load-string "!rm  ./log")
;;
(ti-menu-load-string "/exit OK")

