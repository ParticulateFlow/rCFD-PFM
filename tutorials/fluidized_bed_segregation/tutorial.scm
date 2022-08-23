;; (C)  2021-22 
;;  Stefan Pirker
;;  Particulate Flow Modelling
;;  Johannes Kepler University, Linz, Austria
;;  www.particulate-flow.at 


(define rCFD_src_dir  "../../src")
(define rCFD_user_src_dir  "./user_src")
(define ANSYS_Fluent_case_dir "./ansys_fluent")
(define ANSYS_Fluent_case_file  "quad_80k")

(define number_of_Timesteps_for_fullCFD_analyse_CFD 2)
(define number_of_Timesteps_for_rCFD_analyse_CFD 5)

(define ANSYS_Fluent_number_of_CFD_episodes 1480)
(define ANSYS_Fluent_simulation_timestep_width 0.005)
(define ANSYS_Fluent_number_of_timesteps_per_episode 1)

(define number_of_rCFD_episodes 735)

(define number_of_cfd_episodes 1050) ;;1000
(define number_of_timesteps_for_CFD_episode 10)

(define i 0) 

(define img_rCFD_1       "vof_small_large_d32_")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (CFD_conv1)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string (format  #f "!cp ~a/~a.cas.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))    
    (ti-menu-load-string (format  #f "!cp ~a/~a.dat.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))    
    
    (ti-menu-load-string (format  #f "!cp ~a/CFD_user.c ."        rCFD_user_src_dir))
	
    (ti-menu-load-string "!rm -r libudf_cfd")
)   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (CFD_conv2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string "/define/user-defined/compiled-functions compile 
        libudf_cfd
        yes
        CFD_user.c
        \"\"
        \"\" ")
)   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (CFD_conv3)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string (format #f "/file/read-case-data ./~a.cas.h5 OK" ANSYS_Fluent_case_file))
	
	(ti-menu-load-string "/define/user-defined/compiled-functions load \"libudf_cfd\"")
	
    (ti-menu-load-string "/define/user-defined/execute-on-demand \"CFD_convert_csv2ip::libudf_cfd\"")
	
	;; check if ip file can be loaded 
	
	(ti-menu-load-string "!cp ./data/ip/0000.ip ./data.ip")
	
	(ti-menu-load-string "/file/interpolate/read-data yes data.ip")
)   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (CFD_conv)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (CFD_conv1)
    (CFD_conv2)
    (CFD_conv3)	
)   


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_P1)       ;; prep.1   load all you need in tutorial's main folder
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string (format  #f "!cp ~a/~a.cas.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))    
    (ti-menu-load-string (format  #f "!cp ~a/~a.dat.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))    
    		
    (ti-menu-load-string (format  #f "!cp ~a/*.h ." rCFD_user_src_dir))
	
    (ti-menu-load-string (format  #f "!cp ~a/rCFD_C2C_prep.c ."    rCFD_src_dir))
    (ti-menu-load-string (format  #f "!cp ~a/*.h ."                rCFD_src_dir))

    (ti-menu-load-string "!rm -r libudf_rcfd_prep")
)   


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_P2)       ;; prep.2    compile libudf_preparation needed for simulation (UDF 4)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string "/define/user-defined/compiled-functions compile 
        libudf_rcfd_prep 
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
        rCFD_free.h
        \"\" ")
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_P3)       ;; prep.3:  load start file, allocate UDMI's and init rCFD
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string (format #f "/file/read-case-data ./~a.cas.h5 OK" ANSYS_Fluent_case_file))
    
    (ti-menu-load-string "/define/user-defined/user-defined-memory 30 q" )
    
    (ti-menu-load-string "/define/user-defined/compiled-functions load \"libudf_rcfd_prep\"")
    
    (ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_init_all::libudf_rcfd_prep\"")
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_P4)       ;; prep.4:  run fluent simulations for analyzing time-scales
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
      (ti-menu-load-string "/define/user-defined/function-hooks/execute-at-end 
        \"rCFD_analyse_CFD::libudf_rcfd_prep\"
		\"\"")

    (do  ((i  0 (+ i  1)))
        ((= i  number_of_Timesteps_for_rCFD_analyse_CFD))

		(ti-menu-load-string (format #f "!cp ./data/ip/~04d.ip ./data.ip" i))
		(ti-menu-load-string "/file/interpolate/read-data yes data.ip")

        (ti-menu-load-string (format #f "/solve/set/time-step ~d" ANSYS_Fluent_simulation_timestep_width)) 
        (ti-menu-load-string "/solve/dual-time-iterate 1 1")
		
    )      
	
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_P5)       ;; prep.5:  write intermediate fluent case
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    (ti-menu-load-string (format #f "/file/write-case-data ./data/tmp/~a_After_Analyzing.cas.h5 OK" ANSYS_Fluent_case_file))
)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_P6)       ;; prep.6:  write tracer_start_pos
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_write_Tracer_Positions::libudf_rcfd_prep\"")
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_P7)       ;; prep.7:  define anysy_fluent model settings for discrete phase,
                        ;;          write intermediate fluent file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ;; Note: injection name "rcfd_tracer" must not be already defined
    
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
        
    (ti-menu-load-string "/define/user-defined/function-hooks/execute-at-end 
        \"rCFD_write_Norms::libudf_rcfd_prep\" 
        \"rCFD_write_C2Cs::libudf_rcfd_prep\"
        \"\""
    )
)
    

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_P8)       ;; prep.8:  run fluent simulations for monitoring C2C paths and norms,
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


    (do  ((i  0 (+ i  1)))
        ((= i  ANSYS_Fluent_number_of_CFD_episodes))
		
		(ti-menu-load-string (format #f "!cp ./data/ip/~04d.ip ./data.ip" i))
		(ti-menu-load-string "/file/interpolate/read-data yes data.ip")

        (ti-menu-load-string (format #f "/solve/set/time-step ~d" ANSYS_Fluent_simulation_timestep_width)) 
        (ti-menu-load-string (format #f "/solve/dual-time-iterate ~d 20" ANSYS_Fluent_number_of_timesteps_per_episode))
    )
)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_P9)       ;; prep.9:  analyze and write recurrence path,
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    (ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_write_Rec::libudf_rcfd_prep\"")
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_P10)      ;; prep.10: free memory, delete injection
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_free_all::libudf_rcfd_prep\"")
    
    (ti-menu-load-string "/define/injections/delete-injection rcfd_tracer")
)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_P11)      ;; prep.11: clean-up folder
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


    (ti-menu-load-string "!rm  ./*.c")
    (ti-menu-load-string "!rm  ./*.h")
    (ti-menu-load-string "!rm  ./*.h5") 
    (ti-menu-load-string "!rm  ./*.tiff")
    (ti-menu-load-string "!rm  ./*.tif")
    (ti-menu-load-string "!rm  ./fluent*")
    (ti-menu-load-string "!rm  ./log")
)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_prep)     	;; prep.1-11
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    (rCFD_P1)
    (rCFD_P2) 
    (rCFD_P3)       
    ;;(rCFD_P4)
    ;;(rCFD_P5) 
    (rCFD_P6)       
    (rCFD_P7)
    (rCFD_P8) 
    (rCFD_P9) 
    (rCFD_P10)      
    (rCFD_P11)      
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_R1)       ;; run.1    clean-up folder and load source files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	(ti-menu-load-string "!rm *.h5")
	(ti-menu-load-string "!rm *.c")
	(ti-menu-load-string "!rm *.h")

    (ti-menu-load-string "!rm -r libudf_rcfd_run")	

    (ti-menu-load-string (format  #f "!cp ~a/~a.cas.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))    
    (ti-menu-load-string (format  #f "!cp ~a/~a.dat.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))    
	
    (ti-menu-load-string (format  #f "!cp ~a/*.c ." rCFD_src_dir))
    (ti-menu-load-string (format  #f "!cp ~a/*.h ." rCFD_src_dir))
    (ti-menu-load-string (format  #f "!cp ~a/*.h ." rCFD_user_src_dir))
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_R2)       ;; run.2    compile everything needed for simulation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
        rCFD_free.h
        \"\" ")
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_R3)       ;; run.3:   load start file, allocate UDMI's and init rCFD
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string (format #f "/file/read-case-data ./~a.cas.h5 OK" ANSYS_Fluent_case_file))
    
    (ti-menu-load-string "/define/user-defined/user-defined-memory 20 q" )
    
    (ti-menu-load-string "/define/user-defined/compiled-functions load \"libudf_rcfd_run\"")

    (ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_init_all::libudf_rcfd_run\"")
)   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_R4)       ;; run.4:   read c2c's
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_read_C2Cs::libudf_rcfd_run\"")
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_R5)       ;; run.5:   run simulation w/o post-processing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (do  ((i  0 (+ i  1)))
        ((= i  number_of_rCFD_episodes))              

        (ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_run::libudf_rcfd_run\"")
    )
)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_R5_post)  ;; run.5:   run simulation with post-processing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ;; set graphics setting 
    (ti-menu-load-string "/display/open-window 1")
    (ti-menu-load-string "/display/set-window 1")
    (ti-menu-load-string "/display/set/picture/driver tiff")
    (ti-menu-load-string "/display/set/picture/color-mode color")
    (ti-menu-load-string "/display/set/picture/invert-background? no")
    (ti-menu-load-string "/display/set/picture/invert-background? no")
    (ti-menu-load-string "/display/set/picture/use-window-resolution? no")
    (ti-menu-load-string "/display/set/picture/x-resolution 1189")
    (ti-menu-load-string "/display/set/picture/y-resolution 642")
    (ti-menu-load-string "/display/set/overlays no")        
        
    (do  ((i  0 (+ i  1)))
        ((= i  number_of_rCFD_episodes))              

        (ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_run::libudf_rcfd_run\"")
        
        (ti-menu-load-string "/display/set/overlays no") 

        (ti-menu-load-string "/display/objects/display contour-1")  
        (ti-menu-load-string "/display/set/overlays yes")
        (ti-menu-load-string "/display/views/restore-view view_4x")
        (ti-menu-load-string "/display/set/lights/lights-on no")        
        (ti-menu-load-string "/display/update-scene/select-geometry yes contour-1-3 no no") 
        (ti-menu-load-string "/display/update-scene/display yes no yes yes no no no no no 0 0 0 0")  
        (ti-menu-load-string "/display/update-scene/transform no yes 0.0 0.0 0.0 no no")    
        (ti-menu-load-string "/display/update-scene/select-geometry no yes contour-1-3")    
		
        (ti-menu-load-string "/display/objects/display contour-4")  
        (ti-menu-load-string "/display/set/overlays yes")
        (ti-menu-load-string "/display/views/restore-view view_4x")
        (ti-menu-load-string "/display/set/lights/lights-on no")       
        (ti-menu-load-string "/display/update-scene/select-geometry yes contour-4-3 no no") 
        (ti-menu-load-string "/display/update-scene/display yes no yes yes no no no no no 0 0 0 0")  
        (ti-menu-load-string "/display/update-scene/transform no yes -0.25 0.0 0.0 no no")    
        (ti-menu-load-string "/display/update-scene/select-geometry no yes contour-4-3")    
		
        (ti-menu-load-string "/display/objects/display contour-2")  
        (ti-menu-load-string "/display/set/overlays yes")
        (ti-menu-load-string "/display/views/restore-view view_4x")
        (ti-menu-load-string "/display/set/lights/lights-on no")       
        (ti-menu-load-string "/display/update-scene/select-geometry yes contour-2-3 no no") 
        (ti-menu-load-string "/display/update-scene/display yes no yes yes no no no no no 0 0 0 0")  
        (ti-menu-load-string "/display/update-scene/transform no yes 0.25 0.0 0.0 no no")    
        (ti-menu-load-string "/display/update-scene/select-geometry no yes contour-2-3")    
		
        (ti-menu-load-string "/display/objects/display contour-5")  
        (ti-menu-load-string "/display/set/overlays yes")
        (ti-menu-load-string "/display/views/restore-view view_4x")
        (ti-menu-load-string "/display/set/lights/lights-on no")       
        (ti-menu-load-string "/display/update-scene/select-geometry yes contour-5-3 no no") 
        (ti-menu-load-string "/display/update-scene/display yes no yes yes no no no no no 0 0 0 0")  
        (ti-menu-load-string "/display/update-scene/transform no yes 0.57 0.0 0.0 no no")    
        (ti-menu-load-string "/display/update-scene/select-geometry no yes contour-5-3")    

        (ti-menu-load-string "!rm temp1.tif")
        (ti-menu-load-string "display/hardcopy temp1.tif")
        (ti-menu-load-string (format #f "! convert temp1.tif ~a~04d.jpg &" img_rCFD_1 i))                

	)
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_R6)       ;; run.6: free memory
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_free_all::libudf_rcfd_run\"")
)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_R7)       ;; run.7: clean-up folder
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string "!mv ./*.jpg ./post")      
    (ti-menu-load-string "!mv ./*.out ./post")      
    (ti-menu-load-string "!mv ./Run.trn ./post")
    
    (ti-menu-load-string "!rm  ./*.c")
    (ti-menu-load-string "!rm  ./*.h")
    (ti-menu-load-string "!rm  ./*.h5") 
    (ti-menu-load-string "!rm  ./*.tiff")
    (ti-menu-load-string "!rm  ./*.tif")
    (ti-menu-load-string "!rm  ./fluent*")
    (ti-menu-load-string "!rm  ./log")
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_run)      ;; run.1-7
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (rCFD_R1)
    (rCFD_R2)
    (rCFD_R3)
    (rCFD_R4)
	(rCFD_R5)
    (rCFD_R6)
    (rCFD_R7)   
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_run_post)      ;; run.1-7 with post-processing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (rCFD_R1)
    (rCFD_R2)
    (rCFD_R3)
    (rCFD_R4)
    (rCFD_R5_post)
    (rCFD_R6)
    (rCFD_R7)   
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
