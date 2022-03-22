;; (C)  2021-22 
;;  Stefan Pirker
;;  Particulate Flow Modelling
;;  Johannes Kepler University, Linz, Austria
;;  www.particulate-flow.at 


(define rCFD_src_dir  "../../src")
(define rCFD_user_src_dir  "./user_src")
(define ANSYS_Fluent_case_dir "./ansys_fluent")
(define ANSYS_Fluent_case_file  "LabScale_FB")

(define number_of_Timesteps_for_rCFD_analyse_CFD 5)

(define ANSYS_Fluent_simulation_timestep 0.0005)
(define ANSYS_Fluent_number_of_simulation_timesteps 3050)

(define number_of_rCFD_episodes 100)

(define number_of_cfd_episodes 1000)
(define number_of_timesteps_for_CFD_episode 10)

(define i 0) 

(define img-ref "gas_solid_mixing_")
(define img-1 "gas_phase_")
(define img-2 "solid_phase_")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (CFD_r1)       	;; ref.1   load all you need in tutorial's main folder, create /ref folder
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string (format  #f "!cp ~a/~a.cas.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))    
    (ti-menu-load-string (format  #f "!cp ~a/~a.dat.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))    
    
    (ti-menu-load-string (format  #f "!cp ~a/CFD_ref.c ." rCFD_user_src_dir))

    (ti-menu-load-string "!rm -r libudf_ref")
	
	(ti-menu-load-string "!mkdir ref")
)   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (CFD_r2)       ;; ref.2   compile everything needed for CFD reference case
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string "/define/user-defined/compiled-functions compile 
        libudf_ref 
        yes
        CFD_ref.c
        \"\"
        \"\" ")
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (CFD_r3)       ;; ref.3:  load start file, allocate UDMI's and hook udf's for reference case
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string (format #f "/file/read-case-data ./~a.cas.h5 OK" ANSYS_Fluent_case_file))
    
    (ti-menu-load-string "/define/user-defined/user-defined-memory 30 q" )
    
    (ti-menu-load-string "/define/user-defined/compiled-functions load \"libudf_ref\"")
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (CFD_r4)       ;; ref.4:  hook monitors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
 	(ti-menu-load-string "/define/user-defined/function-hooks/execute-at-end 
		\"CFD_ref_monitors::libudf_ref\" 
		\"\"")

)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (CFD_r5)       ;; ref. 5:   run simulation
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
    
    (do  ((i  0 (+ i  1)))
        ((= i  number_of_CFD_episodes))              

		(ti-menu-load-string (format #f "/solve/set/time-step ~d" ANSYS_Fluent_simulation_timestep)) 
		(ti-menu-load-string (format #f "/solve/dual-time-iterate ~d 30" number_of_timesteps_for_CFD_episode))

        (ti-menu-load-string "/display/objects/display contour-2")  
        (ti-menu-load-string "/display/set/overlays yes")
        (ti-menu-load-string "/display/views/restore-view view_x=0")
        (ti-menu-load-string "/display/set/lights/lights-on no")        
        (ti-menu-load-string "/display/update-scene/select-geometry yes contour-2-5 no no") 
        (ti-menu-load-string "/display/update-scene/display yes no yes no no no no no no 0 0 0 0")  
        (ti-menu-load-string "/display/update-scene/transform no yes 0 0.0 0 no no")    
        (ti-menu-load-string "/display/update-scene/select-geometry no yes contour-2-5")    
		
        (ti-menu-load-string "/display/objects/display contour-3")  
        (ti-menu-load-string "/display/set/overlays yes")
        (ti-menu-load-string "/display/views/restore-view view_x=0")
        (ti-menu-load-string "/display/set/lights/lights-on no")        
        (ti-menu-load-string "/display/update-scene/select-geometry yes contour-3-5 no no") 
        (ti-menu-load-string "/display/update-scene/display yes no yes no no no no no no 0 0 0 0")  
        (ti-menu-load-string "/display/update-scene/transform no yes 0 0.2 0 no no")    
        (ti-menu-load-string "/display/update-scene/select-geometry no yes contour-3-5")    

        (ti-menu-load-string "!rm temp1.tif")
        (ti-menu-load-string "display/hardcopy temp1.tif")
        (ti-menu-load-string (format #f "! convert temp1.tif ~a~04d.jpg &" img-ref i))                
		
    )
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (CFD_r6)       ;; run.7: clean-up folder
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string "!mv ./*.out ./ref")
    (ti-menu-load-string "!mv ./*.jpg ./ref")	
    
    (ti-menu-load-string "!rm  ./*.c")
    (ti-menu-load-string "!rm  ./*.h")
    (ti-menu-load-string "!rm  ./*.h5") 
    (ti-menu-load-string "!rm  ./*.tiff")
    (ti-menu-load-string "!rm  ./*.tif")
    (ti-menu-load-string "!rm  ./fluent*")
    (ti-menu-load-string "!rm  ./log")
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (CFD_ref)       ;; ref   	set up and run reference cfd simulation for validation 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	(CFD_1)
	(CFD_2)
	(CFD_3)
	(CFD_4)
	(CFD_5)
	(CFD_6)
)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_P1)       ;; prep.1   load all you need in tutorial's main folder
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string (format  #f "!cp ~a/~a.cas.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))    
    (ti-menu-load-string (format  #f "!cp ~a/~a.dat.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))    
    
    (ti-menu-load-string (format  #f "!cp ~a/*.c ." rCFD_src_dir))
    (ti-menu-load-string (format  #f "!cp ~a/*.h ." rCFD_src_dir))
    (ti-menu-load-string (format  #f "!cp ~a/*.h ." rCFD_user_src_dir))

    (ti-menu-load-string "!rm -r libudf_rcfd_prep")
    
)   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_P2)       ;; prep.2   compile everything needed for monitoring
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
        \"\""
    )

    (ti-menu-load-string (format #f "/solve/set/time-step ~d" ANSYS_Fluent_simulation_timestep)) 
    (ti-menu-load-string (format #f "/solve/dual-time-iterate ~d 30" number_of_Timesteps_for_rCFD_analyse_CFD))
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
        tracer 
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
    
    ;; At the moment, we do not consider other case-specific (non rCFD) execute-at-end udf's
    
    (ti-menu-load-string "/define/user-defined/function-hooks/execute-at-end 
        \"rCFD_write_Norms::libudf_rcfd_prep\" 
        \"rCFD_write_C2Cs::libudf_rcfd_prep\"
        \"\""
    )

    (ti-menu-load-string (format #f "/file/write-case-data ./data/tmp/~a_Before_Monitoring.cas.h5 OK" ANSYS_Fluent_case_file))
)
    

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_P8)       ;; prep.8:  run fluent simulations for monitoring C2C paths and norms,
                        ;;          write intermediate fluent file
                        ;;          in case of restart, uncomment first lines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ;;(ti-menu-load-string (format #f "/file/read-case-data ./~a_Before_Monitoring.cas.h5 OK" ANSYS_Fluent_case_file))
    
    ;;(ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_init_all::libudf_rcfd_prep\"")

    (ti-menu-load-string (format #f "/solve/set/time-step ~d" ANSYS_Fluent_simulation_timestep)) 
    (ti-menu-load-string (format #f "/solve/dual-time-iterate ~d 30" ANSYS_Fluent_number_of_simulation_timesteps))

    (ti-menu-load-string (format #f "/file/write-case-data ./data/tmp/~a_After_Monitoring.cas.h5 OK" ANSYS_Fluent_case_file))   
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

    (ti-menu-load-string "!rm -r ./*.c")
    (ti-menu-load-string "!rm -r ./*.h")
    (ti-menu-load-string "!rm -r ./*.h5")   
)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_prep)     ;; prep.1-11
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    (rCFD_P1)
    (rCFD_P2) 
    (rCFD_P3)       
    (rCFD_P4)
    (rCFD_P5) 
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

    (ti-menu-load-string (format  #f "!cp ~a/~a.cas.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))    
    (ti-menu-load-string (format  #f "!cp ~a/~a.dat.h5 ." ANSYS_Fluent_case_dir ANSYS_Fluent_case_file))    
    
    (ti-menu-load-string (format  #f "!cp ~a/*.c ." rCFD_src_dir))
    (ti-menu-load-string (format  #f "!cp ~a/*.h ." rCFD_src_dir))
    (ti-menu-load-string (format  #f "!cp ~a/*.h ." rCFD_user_src_dir))

    (ti-menu-load-string "!rm -r libudf_rcfd_run")
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
    
    (ti-menu-load-string "/define/user-defined/user-defined-memory 30 q" )
    
    (ti-menu-load-string "/define/user-defined/compiled-functions load \"libudf_rcfd_run\"")

    (ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_init_all::libudf_rcfd_run\"")
)   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_R4)       ;; run.4:   read c2c's
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_read_C2Cs::libudf_rcfd_run\"")
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_R5)       ;; run.5:   run simulation with post-processing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (do  ((i  0 (+ i  1)))
        ((= i  number_of_rCFD_episodes))              

        (ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_run::libudf_rcfd_run\"")
    )
)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_R5_post)       ;; run.5:   run simulation with post-processing
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
        
        ;; gas phase and solid fraction
        (ti-menu-load-string "/display/set/overlays no")            
        
        (ti-menu-load-string "/display/objects/display contour-2")  
        (ti-menu-load-string "/display/set/overlays yes")
        (ti-menu-load-string "/display/views/restore-view view_x=0")
        (ti-menu-load-string "/display/set/lights/lights-on no")        
        (ti-menu-load-string "/display/update-scene/select-geometry yes contour-2-5 no no") 
        (ti-menu-load-string "/display/update-scene/display yes no yes no no no no no no 0 0 0 0")  
        (ti-menu-load-string "/display/update-scene/transform no yes 0 0.0 0 no no")    
        (ti-menu-load-string "/display/update-scene/select-geometry no yes contour-2-5")    

        (ti-menu-load-string "/display/objects/display contour-3")  
        (ti-menu-load-string "/display/set/overlays yes")
        (ti-menu-load-string "/display/views/restore-view view_x=0")
        (ti-menu-load-string "/display/set/lights/lights-on no")        
        (ti-menu-load-string "/display/update-scene/select-geometry yes contour-3-5 no no") 
        (ti-menu-load-string "/display/update-scene/display yes no yes no no no no no no 0 0 0 0")  
        (ti-menu-load-string "/display/update-scene/transform no yes 0 0.2 0 no no")    
        (ti-menu-load-string "/display/update-scene/select-geometry no yes contour-3-5")    
        
        (ti-menu-load-string "/display/objects/display contour-1")  
        (ti-menu-load-string "/display/set/overlays yes")
        (ti-menu-load-string "/display/views/restore-view view_x=0")
        (ti-menu-load-string "/display/set/lights/lights-on no")        
        (ti-menu-load-string "/display/update-scene/select-geometry yes contour-1-5 no no") 
        (ti-menu-load-string "/display/update-scene/display yes no yes no no no no no no 0 0 0 0")  
        (ti-menu-load-string "/display/update-scene/transform no yes 0 0.4 0 no no")    
        (ti-menu-load-string "/display/update-scene/select-geometry no yes contour-1-5") 
        
        (ti-menu-load-string "!rm temp1.tif")
        (ti-menu-load-string "display/hardcopy temp1.tif")
        (ti-menu-load-string (format #f "! convert temp1.tif ~a~04d.jpg &" img-1 i))                

        ;; solid phase
        (ti-menu-load-string "/display/set/overlays no")        
        
        (ti-menu-load-string "/display/objects/display contour-4")  
        (ti-menu-load-string "/display/set/overlays yes")
        (ti-menu-load-string "/display/views/restore-view view_x=0")
        (ti-menu-load-string "/display/set/lights/lights-on no")        
        (ti-menu-load-string "/display/update-scene/select-geometry yes contour-4-5 no no") 
        (ti-menu-load-string "/display/update-scene/display yes no yes no no no no no no 0 0 0 0")  
        (ti-menu-load-string "/display/update-scene/transform no yes 0 0.0 0 no no")    
        (ti-menu-load-string "/display/update-scene/select-geometry no yes contour-4-5")    

        (ti-menu-load-string "/display/objects/display contour-5")  
        (ti-menu-load-string "/display/set/overlays yes")
        (ti-menu-load-string "/display/views/restore-view view_x=0")
        (ti-menu-load-string "/display/set/lights/lights-on no")        
        (ti-menu-load-string "/display/update-scene/select-geometry yes contour-5-5 no no") 
        (ti-menu-load-string "/display/update-scene/display yes no yes no no no no no no 0 0 0 0")  
        (ti-menu-load-string "/display/update-scene/transform no yes 0 0.2 0 no no")    
        (ti-menu-load-string "/display/update-scene/select-geometry no yes contour-5-5")    
        
        (ti-menu-load-string "/display/objects/display contour-6")  
        (ti-menu-load-string "/display/set/overlays yes")
        (ti-menu-load-string "/display/views/restore-view view_x=0")
        (ti-menu-load-string "/display/set/lights/lights-on no")        
        (ti-menu-load-string "/display/update-scene/select-geometry yes contour-6-5 no no") 
        (ti-menu-load-string "/display/update-scene/display yes no yes no no no no no no 0 0 0 0")  
        (ti-menu-load-string "/display/update-scene/transform no yes 0 0.4 0 no no")    
        (ti-menu-load-string "/display/update-scene/select-geometry no yes contour-6-5") 
                
        (ti-menu-load-string "!rm temp2.tif")
        (ti-menu-load-string "display/hardcopy temp2.tif")
        (ti-menu-load-string (format #f "! convert temp2.tif ~a~04d.jpg &" img-2 i))        
    )
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_R6)       ;; run.6: free memory
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string "/define/user-defined/execute-on-demand \"rCFD_free_all::libudf_rcfd_run\"")
    
    (ti-menu-load-string "/define/injections/delete-injection rcfd_tracer")
)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (rCFD_R7)       ;; run.7: clean-up folder
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    (ti-menu-load-string "!mv ./*.jpg ./post")      
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
