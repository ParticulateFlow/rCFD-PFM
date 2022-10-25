;; (C)  2022
;;  Daniel Queteschiner
;;  Particulate Flow Modelling
;;  Johannes Kepler University, Linz, Austria
;;  www.particulate-flow.at
;;
(define rCFD_src_dir "../../src")
(define rCFD_user_src_dir "./user_src")
(define ANSYS_Fluent_case_dir "./ansys_fluent")
(define ANSYS_Fluent_case_file "LabScale_FB")
(define number_of_rCFD_episodes 15000)  ;; (5 * 15 s * 200 = 15000)
(define i 0)
;;
(ti-menu-load-string (format  #f "!cp ~a/~a.cas.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))
(ti-menu-load-string (format  #f "!cp ~a/~a.dat.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))
(ti-menu-load-string (format  #f "!cp ~a/*.c ." rCFD_src_dir))
(ti-menu-load-string (format  #f "!cp ~a/*.h ." rCFD_src_dir))
(ti-menu-load-string (format  #f "!cp ~a/*.h ." rCFD_user_src_dir))
(ti-menu-load-string "!rm -r libudf_rcfd_run")
;;
(ti-menu-load-string "/define/user-defined/compiled-functions compile
    libudf_rcfd_run
    yes
    rCFD_C2C_run.c
    \"\"
    rCFD_user.h
    rCFD_types.h
    rCFD_globals.h
    rCFD_parallel.h
    rCFD_defaults.h
    rCFD_user_defaults.h
    rCFD_macros.h
    rCFD_init.h
    rCFD_memory.h
    rCFD_layer.h
    rCFD_free.h
    \"\" "
)
;;
(ti-menu-load-string (format #f "/file/read-case-data ./~a.cas.h5 OK" ANSYS_Fluent_case_file))
(ti-menu-load-string "/define/user-defined/user-defined-memory 30 q" )
(ti-menu-load-string "/define/user-defined/compiled-functions load \"libudf_rcfd_run\"")
(ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_init_all::libudf_rcfd_run\"")
;;
(ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_read_C2Cs::libudf_rcfd_run\"")
;;
(do  ((i  0 (+ i  1)))
    ((= i  number_of_rCFD_episodes))

    (ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_run::libudf_rcfd_run\"")
)
;;
(ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_free_all::libudf_rcfd_run\"")
(ti-menu-load-string "/define/injections/delete-injection rcfd_tracer")
;;
(ti-menu-load-string "!mv ./*.jpg ./post")
(ti-menu-load-string "!mv ./Run.trn ./post")
(ti-menu-load-string "!rm  ./*.c")
(ti-menu-load-string "!rm  ./*.h")
(ti-menu-load-string "!rm  ./*.h5")
(ti-menu-load-string "!rm  ./*.tiff")
(ti-menu-load-string "!rm  ./*.tif")
(ti-menu-load-string "!rm  ./fluent*")
(ti-menu-load-string "!rm  ./log")
;;
(ti-menu-load-string "/exit OK")

