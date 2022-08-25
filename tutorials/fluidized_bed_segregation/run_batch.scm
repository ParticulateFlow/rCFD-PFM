;; (C)  2022
;;  Stefan Pirker
;;  Particulate Flow Modelling
;;  Johannes Kepler University, Linz, Austria
;;  www.particulate-flow.at

(ti-menu-load-string (format  #f "!cp ~a/~a.cas.h5 ." "./ansys_fluent" "quad_80k"))
(ti-menu-load-string (format  #f "!cp ~a/~a.dat.h5 ." "./ansys_fluent" "quad_80k"))

(ti-menu-load-string (format  #f "!cp ~a/*.c ." "../../src"))
(ti-menu-load-string (format  #f "!cp ~a/*.h ." "../../src"))
(ti-menu-load-string (format  #f "!cp ~a/*.h ." "./user_src"))

(ti-menu-load-string (format  #f "!mv ./rCFD_user_batch.h ./rCFD_user.h"))

(ti-menu-load-string "!rm -r libudf_rcfd_run")

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
    rCFD_macros.h
    rCFD_init.h
    rCFD_layer.h
    rCFD_utilities.h
    rCFD_free.h
    \"\" ")

(ti-menu-load-string (format #f "/file/read-case-data ./~a.cas.h5 OK" "quad_80k"))

(ti-menu-load-string "/define/user-defined/user-defined-memory 30 q" )

(ti-menu-load-string "/define/user-defined/compiled-functions load \"libudf_rcfd_run\"")

(ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_init_all::libudf_rcfd_run\"")

(ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_read_C2Cs::libudf_rcfd_run\"")

(ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_run::libudf_rcfd_run\"")

(ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_free_all::libudf_rcfd_run\"")

(ti-menu-load-string "!mv ./run_batch.trn ./post")

(ti-menu-load-string "!mv ./Run.trn ./post")

(ti-menu-load-string "!rm -r ./*.c")
(ti-menu-load-string "!rm -r ./*.h")
(ti-menu-load-string "!rm -r ./*.h5")

(ti-menu-load-string "!rm -r ./log")
(ti-menu-load-string "!rm -r ./*.log")

(ti-menu-load-string "!rm -r ./*.trn")

(ti-menu-load-string "!rm -r ./*.sh")
