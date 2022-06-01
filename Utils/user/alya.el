;; $Id: alya.el,v 1.2 2007/01/23
;; 
;; Simple mode for Alya input files
;;
;; Author: G. Houzeaux
;;         based on the tochnog template done by Martin Lüthi, tn@tnoo.net
;;         and Miguel A. Pasenau de Riera, miguel@cimne.upc.es
;;
;;
;; Version 0.1, 2007.01.23
;;
;; 
;; Add to your .emacs
;;
;; (load-library "alya")
;; (load-alya-modes)
;;
;; or
;; 
;; (autoload 'alya-dat-mode "alya" "Alya mode" t)
;; (setq auto-mode-alist (append '(("\\.dat$" . alya-dat-mode)) auto-mode-alist))
;;
;; (autoload 'alya-fix-mode "alya" "Alya .fix mode" t)
;; (setq auto-mode-alist (append '(("\\.fix$" . alya-fix-mode)) auto-mode-alist))
;; 
;;
;; (-*- alya -*-)
;;

(defun fortran-fontify-string (limit)
  (let ((match (match-string 1)))
    (cond ((string= "'" match)
	   (re-search-forward "\\([^'\n]*'?\\)" limit))
	  ((string= "\"" match)
	   (re-search-forward "\\([^\"\n]*\"?\\)" limit)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; For Alya fixity file (boundary conditions)
;;

(defvar alya-fix-mode-hook nil)

;; all Alya keywords
(defvar alya-fix-font-lock-keywords
  (list

   '("END_VALUES"                 . font-lock-function-name-face)
   '("END_FUNCTIONS"              . font-lock-function-name-face)
   '("END_FUNCTION"               . font-lock-function-name-face)
   '("END_ON_NODES"               . font-lock-function-name-face)
   '("END_ON_BOUNDARIES"          . font-lock-function-name-face)
   '("VALUES"                     . font-lock-function-name-face)
   '("FUNCTIONS"                  . font-lock-function-name-face)
   '("FUNCTION"                   . font-lock-function-name-face)
   '("ON_NODES"                   . font-lock-function-name-face)
   '("ON_BOUNDARIES"              . font-lock-function-name-face)
   
   )
  )

(defvar alya-comment-prefix "(** "
  "*The comment `alya-insert-comment' inserts."
  )


;;;###autoload
(defun alya-fix-mode ()
  "Major mode for editing Alya input files."
  (interactive)

  ;; set up local variables
  (kill-all-local-variables)
  (make-local-variable 'font-lock-defaults)
  (make-local-variable 'paragraph-separate)
  (make-local-variable 'paragraph-start)
  (make-local-variable 'require-final-newline)
  (make-local-variable 'comment-start)
  (make-local-variable 'comment-end)
  (make-local-variable 'comment-start-skip)
  (make-local-variable 'comment-column)
  (make-local-variable 'comment-indent-function)
  (make-local-variable 'indent-region-function)
  (make-local-variable 'indent-line-function)
  (make-local-variable 'add-log-current-defun-function)

  (setq font-lock-defaults '(alya-fix-font-lock-keywords nil t))
  (setq major-mode          'alya-fix-mode
	mode-name               "Alya .fix"
	font-lock-defaults      '(alya-fix-font-lock-keywords nil t)
	paragraph-separate      "^[ \t]*$"
	paragraph-start         "^[ \t]*$"
	require-final-newline   t
	comment-start           ""
	comment-end             ""
	comment-start-skip      ""
	comment-column          40
	comment-indent-function 'alya-comment-indent-function
	indent-region-function  'alya-indent-region
	indent-line-function    'alya-indent-line
	)

  ;; Run the mode hook.
  (if alya-fix-mode-hook
      (run-hooks 'alya-fix-mode-hook))
  )

(provide 'alya-fix-mode)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Data and fixity
;;

(defvar alya-dat-mode-hook nil)
(defvar alya-dat-font-lock-keywords
  (list

;;
;; All fonts:
;; ----------
;; font-lock-string-face         rose
;; font-lock-comment-face        dark red
;; font-lock-variable-name-face  yellow
;; font-lock-type-face           light green
;; font-lock-function-name-face  light blue
;; font-lock-defaults            white
;;

;;
;; POSTALYA.DAT FILE
;;

   '("MARK_ELEMENTS"             . font-lock-type-face)
   '("FORMAT"                    . font-lock-type-face)
   '("ELIMINATE_BOUNDARY_NODES"  . font-lock-type-face)
   '("MULTIPLE_FILE"             . font-lock-type-face)
   '("BOUNDARY "                 . font-lock-type-face)
   '("BOUNDARY:"                 . font-lock-type-face)
 
;;
;; GENERAL
;;
   '("$-"                        . font-lock-comment-face)
   '("-"                         . font-lock-comment-face)
   '("ALYA:"                     . font-lock-comment-face)
   '("INCLUDE"                   . font-lock-variable-name-face)
   '("BAR02"                     . font-lock-variable-name-face)
   '("BAR03"                     . font-lock-variable-name-face)
   '("BAR04"                     . font-lock-variable-name-face) 
   '("TRI03"                     . font-lock-variable-name-face) 
   '("TRI06"                     . font-lock-variable-name-face) 
   '("QUA04"                     . font-lock-variable-name-face) 
   '("QUA08"                     . font-lock-variable-name-face) 
   '("QUA09"                     . font-lock-variable-name-face) 
   '("QUA16"                     . font-lock-variable-name-face) 
   '("TET04"                     . font-lock-variable-name-face) 
   '("TET10"                     . font-lock-variable-name-face) 
   '("PYR05"                     . font-lock-variable-name-face) 
   '("PYR14"                     . font-lock-variable-name-face) 
   '("PEN06"                     . font-lock-variable-name-face) 
   '("PEN15"                     . font-lock-variable-name-face) 
   '("PEN18"                     . font-lock-variable-name-face) 
   '("HEX08"                     . font-lock-variable-name-face) 
   '("HEX20"                     . font-lock-variable-name-face) 
   '("HEX27"                     . font-lock-variable-name-face) 
   '("HEX64"                     . font-lock-variable-name-face) 
;;
;; IMPORTANT WORDS
;;
   '("BINARY"                    . font-lock-variable-name-face)
   '("BIN"                       . font-lock-variable-name-face)
   '("ASCII"                     . font-lock-variable-name-face)
   '("FINE"                      . font-lock-variable-name-face)
   '("COARSE"                    . font-lock-variable-name-face)
;;
;; MUST BE BEFORE BECAUSE THEY INCLUDE COLORED WORD
;;
   '("END_DISCRETE_FUNCTIONS"     . font-lock-type-face)
   '("END_SPACE_&_TIME_FUNCTIONS" . font-lock-type-face)
   '("END_TIME_FUNCTIONS"         . font-lock-type-face)
   '("END_SPACE_FUNCTIONS"        . font-lock-type-face)
   '("DISCRETE_FUNCTIONS"         . font-lock-type-face)
   '("SPACE_FUNCTIONS"            . font-lock-type-face)
   '("UNKNOWN_BOUNDARIES"         . font-lock-type-face)
   '("END_FUNCTIONS"              . font-lock-type-face)
   '("END_NODES"                  . font-lock-type-face)
   '("END_ON_NODES"               . font-lock-type-face)
   '("END_ON_BOUNDARIES"          . font-lock-type-face)
   '("END_GEOMETRICAL_CONDITIONS" . font-lock-type-face)
   '("END_SKEW_SYSTEMS"           . font-lock-type-face)
   '("END_MATERIALS"              . font-lock-comment-face)
   '("END_MATERIAL"               . font-lock-comment-face)
   '("SPACE_&_TIME_FUNCTIONS"     . font-lock-type-face)
   '("TIME_FUNCTIONS"             . font-lock-type-face)
   '("FUNCTIONS"                  . font-lock-type-face)
   '("ON_NODES"                   . font-lock-type-face)
   '("ON_BOUNDARIES"              . font-lock-type-face)
   '("GEOMETRICAL_CONDITIONS"     . font-lock-type-face)
   '("PERIODIC_NODES"             . font-lock-type-face)
   '("END_NODES_PER_ELEMENT"      . font-lock-type-face)
   '("NODES_PER_ELEMENT"          . font-lock-type-face)
   '("END_BOUNDARIES"             . font-lock-type-face)
;;
;; NASTIN
;;
   '("TAU_STRATEGY"              . font-lock-defaults)
   '("SCHUR_COMPLEMENT"          . font-lock-variable-name-face)
   '("FRACTIONAL_STEP"           . font-lock-variable-name-face)
   '("SEMI_IMPLICIT"             . font-lock-variable-name-face)
   '("MOMENTUM_PRESERVING"       . font-lock-variable-name-face)
   '("CONTINUITY_PRESERVING"     . font-lock-variable-name-face)
   '("END_MOMENTUM"              . font-lock-type-face)
   '("END_VISCOUS_TERM"          . font-lock-type-face)
   '("END_HYDROSTATIC_PRESSURE"  . font-lock-type-face)
   '("END_CONTINUITY"            . font-lock-type-face)
   '("END_ALGORITHM"             . font-lock-type-face)
   '("MOMENTUM"                  . font-lock-type-face)
   '("END_VISCOUS_TERM"          . font-lock-type-face)
   '("HYDROSTATIC_PRESSURE"      . font-lock-type-face)
   '("CONTINUITY"                . font-lock-type-face)
   '("ALGORITHM"                 . font-lock-type-face)

;   TEMPORAL_DERIVATIVES
;   CONVECTIVE_TERM:          On
;   VISCOUS_TERM:             LAPLACIAN   
;   '("ELEMENT_LENGTH"               . font-lock-variable-name-face)
;   '("STABILIZATION"                . font-lock-variable-name-face)  
;   '("TRACKING"                     . font-lock-variable-name-face)
;   '("TIME_INTEGRATION"             . font-lock-variable-name-face)  
;   '("SAFETY_FACTOR"                . font-lock-variable-name-face)   
;   '("STEADY_STATE_TOLER"           . font-lock-variable-name-face)   
;   '("STEADY_STATE_TOLERANCE"       . font-lock-variable-name-face)   
;   '("VECTORIZED_ASSEMBLY"          . font-lock-variable-name-face)   
;   '("NORM_OF_CONVERGENCE"          . font-lock-variable-name-face)   
;   '("MAXIMUM_NUMBER_OF_IT"         . font-lock-variable-name-face)
;   '("MAXIMUM_NUMBER_OF_ITERATIONS" . font-lock-variable-name-face)
;   '("CONVERGENCE_TOLERANCE"        . font-lock-variable-name-face)  

;;
;; DOMAIN
;;
   '("NODAL_POINTS"              . font-lock-type-face)
   '("SPACE_DIMENSIONS"          . font-lock-type-face)
   '("TYPES_OF"                  . font-lock-type-face)
   '("TYPES_OF_ELEMENTS"         . font-lock-type-face)
   '("NODES"                     . font-lock-type-face) 
   '("EDGES"                     . font-lock-type-face) 
   '("BOUNDARIES"                . font-lock-type-face)
   '("SKEW_SYSTEMS"              . font-lock-type-face)
   '("MATERIALS="                . font-lock-type-face)
   '("FIELDS="                   . font-lock-type-face)
   '("SUBDOMAINS="               . font-lock-type-face)
   '("GROUPS="                   . font-lock-type-face)
   '("MATERIALS"                 . font-lock-comment-face)
   '("MATERIAL"                  . font-lock-comment-face)
   '("INTEGRATION_RULE"          . font-lock-type-face)
   '("END_IMMERSED_BOUNDARY"     . font-lock-type-face)
   '("IMMERSED_BOUNDARY"         . font-lock-type-face)
   '("DOMAIN_INTEGRATION_POINTS" . font-lock-type-face)
   '("SCALE"                     . font-lock-type-face)
   '("OUTPUT_MESH_DATA"          . font-lock-type-face)
   '("OUTPUT:"                   . font-lock-type-face)
   '("END_ELSEST"                . font-lock-type-face)
   '("ELSEST"                    . font-lock-type-face)
   '("PERIODICITY_STRATEGY"      . font-lock-type-face)
   '("END_DIMENSIONS"            . font-lock-function-name-face)
   '("END_STRATEGY"              . font-lock-function-name-face)
   '("DIMENSIONS"                . font-lock-function-name-face)
   '("STRATEGY"                  . font-lock-function-name-face)
   '("END_SETS"                  . font-lock-function-name-face)
   '("SETS"                      . font-lock-function-name-face)
   '("BOUNDARY_ELEMENT"          . font-lock-type-face)
   '("EXTRAPOLATE_BOUNDARY_CONDITIONS" . font-lock-type-face)
;;
;; GEOMETRY
;;
   '("END_GEOMETRY"              . font-lock-function-name-face)
   '("END_TYPES_OF_ELEMENTS"     . font-lock-type-face)
   '("END_PERIODICITY"           . font-lock-type-face)
   '("END_ELEMENTS"              . font-lock-type-face)
   '("END_COORDINATES"           . font-lock-type-face)
   '("END_TYPES"                 . font-lock-type-face)
   '("END_MATERIALS"             . font-lock-type-face)
   '("END_SUBDOMAINS"            . font-lock-type-face)
   '("END_LELBO"                 . font-lock-type-face)
   '("GEOMETRY"                  . font-lock-function-name-face)
   '("ELEMENTS"                  . font-lock-type-face)
   '("TYPES_OF_ELEMENTS"         . font-lock-type-face)
   '("BOUNDARIES"                . font-lock-type-face)
   '("COORDINATES"               . font-lock-type-face)
   '("LELBO"                     . font-lock-type-face)
   '("TYPES"                     . font-lock-type-face)
   '("MATERIALS"                 . font-lock-type-face)
   '("SUBDOMAINS"                . font-lock-type-face)
   '("ELEMENTS"                  . font-lock-type-face)
   '("PERIODICITY"               . font-lock-type-face)
;;
;; FIELDS
;;
   '("END_FIELDS"                . font-lock-function-name-face)
   '("FIELDS"                    . font-lock-function-name-face)
   '("END_FIELD"                 . font-lock-variable-name-face) 
   '("FIELD"                     . font-lock-variable-name-face)   
;;
;; PROPERTIES
;;
   '("DENSITY"                   . font-lock-variable-name-face)
   '("VISCOSITY"                 . font-lock-variable-name-face)
   '("SPECIFIC_HEAT"             . font-lock-variable-name-face)
   '("CONDUCTIVITY"              . font-lock-variable-name-face)   
   '("DIAMETER"                  . font-lock-variable-name-face)   
   '("MODEL"                     . font-lock-variable-name-face)   
   '("NEWMARK"                   . font-lock-variable-name-face)   
;;
;; SOLVERS
;;
   '("DEFLATED_CG"               . font-lock-string-face)
   '("CG"                        . font-lock-string-face)
   '("GMRES"                     . font-lock-string-face)
   '("BICGSTAB"                  . font-lock-string-face)
   '("GAUSS_SEIDEL"              . font-lock-string-face)
   '("LINELET"                   . font-lock-string-face)
   '("RAS"                       . font-lock-string-face)
   '("DIAGONAL"                  . font-lock-string-face)
   '("ORTHOMIN"                  . font-lock-string-face)
   '("RICHARDSON"                . font-lock-string-face)
   '("BLOCK_DIAGONAL"            . font-lock-string-face)
   '("AII"                       . font-lock-string-face)
   '("CONJUGATE_GRADIENT"        . font-lock-string-face)
   '("AGMG"                      . font-lock-string-face)
   '("SPARSE"                    . font-lock-string-face)
   '("DIRECT"                    . font-lock-string-face)
   '("CSR"                       . font-lock-string-face)
   '("SKYLINE"                   . font-lock-string-face)
   '("CHOLESKY"                  . font-lock-string-face)
   '("KRYLOV"                    . font-lock-string-face)
   '("STATIC"                    . font-lock-string-face)
   '("DYNAMIC"                   . font-lock-string-face)
   '("CHUNK_SIZE"                . font-lock-variable-name-face)
   '("AUTO_TUNING"               . font-lock-variable-name-face)
;;
;; RUN DATA
;;
   '("ON_LAST_MESH"              . font-lock-variable-name-face)
   '("ON_ORIGINAL_MESH"          . font-lock-variable-name-face)
   '("GRAVITY"                   . font-lock-variable-name-face)
   '("ABSORPTION"                . font-lock-variable-name-face)
   '("SCATTERING"                . font-lock-variable-name-face)
   '("DIVISION"                  . font-lock-variable-name-face)
   '("END_MULTIPLICATION"        . font-lock-variable-name-face)
   '("MULTIPLICATION"            . font-lock-variable-name-face)
   '("ROUGHNESS"                 . font-lock-variable-name-face)
   '("U_REFERENCE"               . font-lock-variable-name-face)
   '("H_REFERENCE"               . font-lock-variable-name-face)
   '("K_REFERENCE"               . font-lock-variable-name-face)
   '("US_REFERENCE"              . font-lock-variable-name-face)
   '("WALL_LAW"                  . font-lock-variable-name-face)

   '("END_ALGEBRAIC_SOLVER"      . font-lock-type-face)
   '("ALGEBRAIC_SOLVER"          . font-lock-type-face)

   '("PRECONDITIONING"           . font-lock-variable-name-face)
   '("PRECONDITIONER"            . font-lock-variable-name-face)
   '("EXTENDED_GRAPH"            . font-lock-variable-name-face)

   '("RUN_TYPE"                  . font-lock-type-face)
   '("LATEX_INFO_FILE"           . font-lock-type-face)
   '("CUSTOMER"                  . font-lock-type-face)
   '("POSTPROCESS"               . font-lock-type-face)
   '("OUTPUT_FORMAT"             . font-lock-type-face)
   '("LIVE_INFORMATION"          . font-lock-type-face)
   '("LOG_FILE"                  . font-lock-type-face)
   '("LIVE_INFO"                 . font-lock-type-face)
   '("LOT_OF_MEMORY"             . font-lock-type-face)
   '("MAXIMUM_MEMORY"            . font-lock-variable-name-face)
   '("MEMORY"                    . font-lock-type-face)
   '("TIME_COUPLING"             . font-lock-type-face)
   '("TIME_INTERVAL"             . font-lock-type-face)
   '("TIME_STEP_SIZE"            . font-lock-type-face)
   '("NUMBER_OF_STEPS"           . font-lock-type-face)
   '("MAXIMUM_NUMBER_GLOBAL"     . font-lock-type-face)
   '("END_BLOCK_ITERATION"       . font-lock-function-name-face)
   '("END_RUN_DATA"              . font-lock-function-name-face)
   '("END_PROBLEM_DATA"          . font-lock-function-name-face)
   '("END_MPI_IO"                . font-lock-function-name-face)
   '("BLOCK_ITERATION"           . font-lock-function-name-face)
   '("RUN_DATA"                  . font-lock-function-name-face)
   '("PROBLEM_DATA"              . font-lock-function-name-face)
   '("MPI_IO"                    . font-lock-function-name-face)
   '("END_LAGRANGIAN_PARTICLES"  . font-lock-type-face)
   '("LAGRANGIAN_PARTICLES"      . font-lock-type-face)

   '("EXECUTION_MODE"            . font-lock-type-face)
   '("END_PARTITIONING"          . font-lock-type-face)
   '("PARTITIONING"              . font-lock-type-face)
   '("END_REPARTITIONING"        . font-lock-type-face)
   '("REPARTITIONING"            . font-lock-type-face)
   '("WEIGHT"                    . font-lock-variable-name-face)
   '("SFC"                       . font-lock-variable-name-face)
   '("GID"                       . font-lock-variable-name-face)
   '("METIS"                     . font-lock-variable-name-face)
   '("METIS4"                    . font-lock-variable-name-face)
   '("ONLY_PREPROCESS"           . font-lock-variable-name-face)
   '("READ_PREPROCESS"           . font-lock-variable-name-face)
   '("PARTITION_AND_RUN"         . font-lock-variable-name-face)

;;
;; KERMOD
;;

   '("END_MESH"                  . font-lock-type-face)
   '("MESH"                      . font-lock-type-face)
   '("END_WITNESS_POINTS"        . font-lock-type-face)
   '("WITNESS_POINTS"            . font-lock-type-face)
   '("END_WITNESS"               . font-lock-type-face)
   '("WITNESS"                   . font-lock-type-face)
   '("END_FILTERS"               . font-lock-type-face)
   '("FILTERS"                   . font-lock-type-face)
   
;;
;; RUN DATA: MODULES AND SERVICES
;;
   '("GLOBAL"                    . font-lock-comment-face)
   '("LOCAL"                     . font-lock-comment-face)
   '("PRESCRIBED"                . font-lock-comment-face)
   '("FROM_CRITICAL"             . font-lock-comment-face)

   '("END_NASTIN_MODULE"         . font-lock-comment-face)
   '("END_TEMPER_MODULE"         . font-lock-comment-face)
   '("END_CODIRE_MODULE"         . font-lock-comment-face)
   '("END_TURBUL_MODULE"         . font-lock-comment-face)
   '("END_EXMEDI_MODULE"         . font-lock-comment-face)
   '("END_NASTAL_MODULE"         . font-lock-comment-face)
   '("END_ALEFOR_MODULE"         . font-lock-comment-face)
   '("END_LATBOL_MODULE"         . font-lock-comment-face)
   '("END_APELME_MODULE"         . font-lock-comment-face)
   '("END_SOLIDZ_MODULE"         . font-lock-comment-face)
   '("END_GOTITA_MODULE"         . font-lock-comment-face)
   '("END_WAVEQU_MODULE"         . font-lock-comment-face)
   '("END_LEVELS_MODULE"         . font-lock-comment-face)
   '("END_QUANTY_MODULE"         . font-lock-comment-face)
   '("END_PARTIS_MODULE"         . font-lock-comment-face)
   '("END_CHEMIC_MODULE"         . font-lock-comment-face)
   '("END_HELMOZ_MODULE"         . font-lock-comment-face)
   '("END_IMMBOU_MODULE"         . font-lock-comment-face)
   '("END_RADIAT_MODULE"         . font-lock-comment-face)
   '("END_POROUS_MODULE"         . font-lock-comment-face)
   '("END_NEUTRO_MODULE"         . font-lock-comment-face)

   '("NASTIN_MODULE"             . font-lock-comment-face)
   '("TEMPER_MODULE"             . font-lock-comment-face)
   '("CODIRE_MODULE"             . font-lock-comment-face)
   '("TURBUL_MODULE"             . font-lock-comment-face)
   '("EXMEDI_MODULE"             . font-lock-comment-face)
   '("NASTAL_MODULE"             . font-lock-comment-face)
   '("ALEFOR_MODULE"             . font-lock-comment-face)
   '("LATBOL_MODULE"             . font-lock-comment-face)
   '("APELME_MODULE"             . font-lock-comment-face)
   '("SOLIDZ_MODULE"             . font-lock-comment-face)
   '("GOTITA_MODULE"             . font-lock-comment-face)
   '("WAVEQU_MODULE"             . font-lock-comment-face)
   '("LEVELS_MODULE"             . font-lock-comment-face)
   '("QUANTY_MODULE"             . font-lock-comment-face)
   '("PARTIS_MODULE"             . font-lock-comment-face)
   '("CHEMIC_MODULE"             . font-lock-comment-face)
   '("HELMOZ_MODULE"             . font-lock-comment-face)
   '("IMMBOU_MODULE"             . font-lock-comment-face)
   '("RADIAT_MODULE"             . font-lock-comment-face)
   '("POROUS_MODULE"             . font-lock-comment-face)
   '("NEUTRO_MODULE"             . font-lock-comment-face)
  
   '("NASTIN"                    . font-lock-variable-name-face)
   '("TEMPER"                    . font-lock-variable-name-face)
   '("CODIRE"                    . font-lock-variable-name-face)
   '("TURBUL"                    . font-lock-variable-name-face)
   '("EXMEDI"                    . font-lock-variable-name-face)
   '("NASTAL"                    . font-lock-variable-name-face)
   '("ALEFOR"                    . font-lock-variable-name-face)
   '("LATBOL"                    . font-lock-variable-name-face)
   '("APELME"                    . font-lock-variable-name-face)
   '("SOLIDZ"                    . font-lock-variable-name-face)
   '("GOTITA"                    . font-lock-variable-name-face)
   '("WAVEQU"                    . font-lock-variable-name-face)
   '("LEVELS"                    . font-lock-variable-name-face)
   '("QUANTY"                    . font-lock-variable-name-face)
   '("PARTIS"                    . font-lock-variable-name-face)
   '("CHEMIC"                    . font-lock-variable-name-face)
   '("HELMOZ"                    . font-lock-variable-name-face)
   '("IMMBOU"                    . font-lock-variable-name-face)
   '("RADIAT"                    . font-lock-variable-name-face)
   '("POROUS"                    . font-lock-variable-name-face)

   '("END_HDFPOS_SERVICE"        . font-lock-comment-face)
   '("END_PARALL_SERVICE"        . font-lock-comment-face)
   '("END_DODEME_SERVICE"        . font-lock-comment-face)
   '("HDFPOS_SERVICE"            . font-lock-comment-face)
   '("PARALL_SERVICE"            . font-lock-comment-face)
   '("DODEME_SERVICE"            . font-lock-comment-face)

   '("POSTPROCESS"               . font-lock-type-face)
   '("OUTPUT_FILE"               . font-lock-type-face)
   '("TASK"                      . font-lock-type-face)
   '("PARTITION_TYPE"            . font-lock-type-face)
   '("COMMUNICATIONS"            . font-lock-type-face)
   '("FILE_HIERARCHY"            . font-lock-type-face)
   '("VIRTUAL_FILE"              . font-lock-type-face)

;;
;; MODULES
;;
   '("END_PHYSICAL_PROBLEM"      . font-lock-function-name-face)
   '("END_PROBLEM_DEFINITION"    . font-lock-type-face)
   '("END_PROPERTIES"            . font-lock-type-face)
   '("END_NUMERICAL_TREATMENT"   . font-lock-function-name-face)
   '("END_SUBGRID_SCALE_MODELING". font-lock-type-face)
   '("END_TIME_PROBLEM"          . font-lock-type-face)
   '("END_INNER_ITERATIONS"      . font-lock-type-face)
   '("END_OTHERS"                . font-lock-type-face)

   '("END_OUTPUT_&_POST_PROCESS" . font-lock-function-name-face)
   '("END_OUTPUT_&_POSTPROCESS"  . font-lock-function-name-face)
   '("END_BOUNDARY_CONDITIONS"   . font-lock-function-name-face)

   '("END_ELEMENT_SET"           . font-lock-type-face)
   '("END_BOUNDARY_SET"          . font-lock-type-face)
   '("END_NODE_SET"              . font-lock-type-face)
   '("END_CODES"                 . font-lock-type-face)
   '("END_PARAMETERS"            . font-lock-type-face)

   '("PHYSICAL_PROBLEM"          . font-lock-function-name-face)
   '("PROBLEM_DEFINITION"        . font-lock-type-face)
   '("PROPERTIES"                . font-lock-type-face)
   '("NUMERICAL_TREATMENT"       . font-lock-function-name-face)
   '("SUBGRID_SCALE_MODELING"    . font-lock-type-face)
   '("TIME_PROBLEM"              . font-lock-type-face)
   '("INNER_ITERATIONS"          . font-lock-type-face)
   '("OTHERS"                    . font-lock-type-face)
   '("EDGE_ELEMENTS"             . font-lock-type-face)

   '("OUTPUT_&_POST_PROCESS"     . font-lock-function-name-face)
   '("OUTPUT_&_POSTPROCESS"      . font-lock-function-name-face)
   '("BOUNDARY_CONDITIONS"       . font-lock-function-name-face)

   '("ELEMENT_SET"               . font-lock-type-face)
   '("BOUNDARY_SET"              . font-lock-type-face)
   '("NODE_SET"                  . font-lock-type-face)
   '("CODES"                     . font-lock-type-face)
   '("PARAMETERS"                . font-lock-type-face)

;;
;; CHEMIC
;;
   '("END_REACTIONS"             . font-lock-variable-name-face)
   '("END_CLASSES"               . font-lock-variable-name-face)
   '("END_ODES"                  . font-lock-variable-name-face)
   '("ODES"                      . font-lock-variable-name-face)
   '("CLASSES"                   . font-lock-variable-name-face)
   '("REACTIONS"                 . font-lock-variable-name-face)

   '("END_DATA"                  . font-lock-comment-face)
   '("DATA"                      . font-lock-comment-face)
   
;;
;; VERY SMALL WORDS
;;
   '(" ON"                        . font-lock-keyword-face)
   '(" YES"                       . font-lock-keyword-face)
   '(" OFF"                       . font-lock-keyword-face)
   '(" NO"                        . font-lock-keyword-face)

   )
  )

(defvar alya-comment-prefix "(** "
  "*The comment `alya-insert-comment' inserts."
  )


;;;###autoload
(defun alya-dat-mode ()
  "Major mode for editing Alya input files."
  (interactive)
  ;; set up local variables
  (kill-all-local-variables)
  (make-local-variable 'font-lock-defaults)
  (make-local-variable 'paragraph-separate)
  (make-local-variable 'paragraph-start)
  (make-local-variable 'require-final-newline)
  (make-local-variable 'comment-start)
  (make-local-variable 'comment-end)
  (make-local-variable 'comment-start-skip)
  (make-local-variable 'comment-column)
  (make-local-variable 'comment-indent-function)
  (make-local-variable 'indent-region-function)
  (make-local-variable 'indent-line-function)
  (make-local-variable 'add-log-current-defun-function)
  ;;
  (setq font-lock-defaults '(alya-dat-font-lock-keywords nil t))
  ;;(set-syntax-table alya-dat-mode-syntax-table)
  (setq major-mode          'alya-dat-mode
	mode-name               "Alya .dat"
	font-lock-defaults      '(alya-dat-font-lock-keywords nil t)
	paragraph-separate      "^[ \t]*$"
	paragraph-start         "^[ \t]*$"
	require-final-newline   t
	comment-start           ""
	comment-end             ""
	comment-start-skip      ""
	comment-column          40
	comment-indent-function 'alya-comment-indent-function
	indent-region-function  'alya-indent-region
	indent-line-function    'alya-indent-line
	)

  ;; Run the mode hook.
  (if alya-dat-mode-hook
      (run-hooks 'alya-dat-mode-hook))
  )

(provide 'alya-dat-mode)

(defun load-alya-modes ()
  ;; para que al cargar un .dat/.fix se cargue el modo Alya.  
  (autoload 'alya-fix-mode "alya" "Alya .fix mode" t)
  (setq auto-mode-alist (append '(("\\.fix$" . alya-fix-mode)) auto-mode-alist))
  
  (autoload 'alya-dat-mode "alya" "Alya .dat mode" t)
  (setq auto-mode-alist (append '(("\\.dat$" . alya-dat-mode)) auto-mode-alist))

  (autoload 'alya-dat-mode "alya" "Alya .alyadat mode" t)
  (setq auto-mode-alist (append '(("\\.alyadat$" . alya-dat-mode)) auto-mode-alist))
)



(defun alya-indent-line ()
  )
(setq alya-indent-line 2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 
;; 
;; MENU BAR
;; 
;; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defcustom alya-comment-region "$"
  "String inserted by \\[alya-comment-region] at start of each line in region."
  :type  'string
  :safe  'stringp
  :group 'alya-indent)

(defun alya-separator ()
  "Insert separator at cursor point."
  (interactive)
  (insert "$-------------------------------------------------------------------" ?\n)
  (backward-char 69))

(defun alya-markup-region (start end)
  "Wrap a the region with separators."
  (interactive "r")
  (save-excursion
    (goto-char end)   (insert "$-------------------------------------------------------------------" ?\n)
    (goto-char start) (insert "$-------------------------------------------------------------------" ?\n)
    ))

(defun alya-comment-region (beg-region end-region)
  "Comment/uncomment every line in the region.
Insert the variable `alya-comment-region' at the start of every line
in the region, or, if already present, remove it."
  (interactive "*r")
  (let ((end (copy-marker end-region)))
    (goto-char beg-region)
    (beginning-of-line)
    (if (looking-at (regexp-quote alya-comment-region))
        (delete-region (point) (match-end 0))
      (insert alya-comment-region))
    (while (and (zerop (forward-line 1))
                (< (point) end))
      (if (looking-at (regexp-quote alya-comment-region))
          (delete-region (point) (match-end 0))
        (insert alya-comment-region)))
    (set-marker end nil)))

(defun alya-downcase ()
  (interactive)
  (when (region-active-p)
    (let ((i 0)
	  (return-string "")
	  (input (buffer-substring-no-properties (region-beginning) (region-end))))
      (while (< i (- (region-end) (region-beginning)))
	(let ((current-char (substring input i (+ i 1))))
	  (if (string= (substring input i (+ i 1)) (downcase (substring input i (+ i 1))))
	      (setq return-string
		    (concat return-string (downcase (substring input i (+ i 1)))))
	    (setq return-string
		  (concat return-string (downcase (substring input i (+ i 1)))))))
	(setq i (+ i 1)))
      (delete-region (region-beginning) (region-end))
      (insert return-string))))

(defun alya-upcase ()
  (interactive)
  (when (region-active-p)
    (let ((i 0)
	  (return-string "")
	  (input (buffer-substring-no-properties (region-beginning) (region-end))))
      (while (< i (- (region-end) (region-beginning)))
	(let ((current-char (substring input i (+ i 1))))
	  (if (string= (substring input i (+ i 1)) (upcase (substring input i (+ i 1))))
	      (setq return-string
		    (concat return-string (upcase (substring input i (+ i 1)))))
	    (setq return-string
		  (concat return-string (upcase (substring input i (+ i 1)))))))
	(setq i (+ i 1)))
      (delete-region (region-beginning) (region-end))
      (insert return-string))))

(defun alya-parall ()
  "Insert Parall service options."
  (interactive)
  (insert "PARALL_SERVICE: ON" ?\n)
  (insert "  OUTPUT:          OFF" ?\n)
  (insert "  OPENMP" ?\n)
  (insert "    $ ELEMENT_CHUNK=  500" ?\n)
  (insert "    $ BOUNDARY_CHUNK= 500" ?\n)
  (insert "    $ NODE_CHUNK=     500" ?\n)
  (insert "    COLORING:         MINCOLOR" ?\n)
  (insert "  END_OPENMP" ?\n)
  (insert "  COMMUNICATIONS:  ASYNCHRONOUS" ?\n)
  (insert "  $ TASK:            READ_PREPROCESS, BINARY, SUBDOMAINS = 10 $ ONLY_PREPROCESS" ?\n)
  (insert "  $ FILE_HIERARCHY:  ON $ OFF" ?\n)
  (insert "  $ FILE_OPEN_CLOSE: ON $ OFF" ?\n)
  (insert "  $ VIRTUAL_MEMORY:  ON, MAXIMUM= X" ?\n)
  (insert "  $ PARTITIONING" ?\n)
  (insert "      $ METHOD:         METIS                  $ SFC / FIELD=int / ORIENTED_BIN" ?\n)
  (insert "      $ ELEMENT_GRAPH:  FACES                  $ NODES (ONLY METIS)" ?\n)
  (insert "      $ BOXES:          int,int,int  COARSE    " ?\n)
  (insert "      $ BOXES:          int,int,int  FINE      " ?\n)
  (insert "      $ DIRECTION:      rea,rea,rea            " ?\n)
  (insert "      $ EXECUTION_MODE: PARALLEL               $ SEQUENTIAL (PARALLEL ONLY SFC AND ORUENTED BIN)" ?\n)
  (insert "      $ WEIGHT:         GAUSS                  $ OFF / FIELD / ELEMENT, TRI03=1, QU04=2.4..." ?\n)
  (insert "  $ END_PARTITIONING" ?\n)
  (insert "END_PARALL_SERVICE" ?\n)
  )

(defun alya-mpiio ()
  "Insert MPIO/IO options."
  (interactive)
  (insert "MPI_IO: ON  " ?\n)
  (insert "  GEOMETRY:     On $ OFF / EXPORT " ?\n)
  (insert "  RESTART:      On $ OFF " ?\n)
  (insert "  POSTPPROCESS: On $ OFF " ?\n)
  (insert "  MERGE:        On $ OFF " ?\n)
  (insert "  END_MERGE              " ?\n)
  (insert "END_MPI_IO" ?\n)
  )
  
(defun alya-doxygen-subroutine ()
  "Insert Doxygen header for a subroutine."
  (interactive)
  (insert "!-----------------------------------------------------------------------" ?\n)
  (insert "!> @addtogroup ???" ?\n)
  (insert "!> @{" ?\n)
  (insert "!> @file    "(file-relative-name (buffer-file-name)) ?\n)
  (insert "!> @author  "(shell-command-to-string "echo -n $(whoami)") ?\n)
  (insert "!> @date    "(shell-command-to-string "echo -n $(date +%Y-%m-%d)") ?\n)
  (insert "!> @brief   ???" ?\n)
  (insert "!> @details ???" ?\n)
  (insert "!>          ???" ?\n)
  (insert "!> @} " ?\n)
  (insert "!-----------------------------------------------------------------------" ?\n)
  )

(defun alya-doxygen-defgroup ()
  "Insert Doxygen header for a module."
  (interactive)
  (insert "!-----------------------------------------------------------------------" ?\n)
  (insert "!> @defgroup ???" ?\n)
  (insert "!> ... description of the group... " ?\n)
  (insert "!-----------------------------------------------------------------------" ?\n)
  )
  
(defun alya-doxygen-module ()
  "Insert Doxygen header for a module."
  (interactive)
  (insert "!-----------------------------------------------------------------------" ?\n)
  (insert "!> @addtogroup ???" ?\n)
  (insert "!> @{" ?\n)
  (insert "!> @file    "(file-relative-name (buffer-file-name)) ?\n)
  (insert "!> @author  "(shell-command-to-string "echo -n $(whoami)") ?\n)
  (insert "!> @date    "(shell-command-to-string "echo -n $(date +%Y-%m-%d)") ?\n)
  (insert "!> @brief   ???" ?\n)
  (insert "!> @details ???" ?\n)
  (insert "!>          ???" ?\n)
  (insert "!-----------------------------------------------------------------------" ?\n)
  )

(defun alya-doxygen-module-close ()
  "Insert Doxygen end of a module."
  (interactive)
  (insert "!> @}" ?\n)
  )

(defun alya-doxygen-subroutine-module ()
  "Insert Doxygen header for a module subroutine."
  (interactive)
  (insert "!-----------------------------------------------------------------------" ?\n)
  (insert "!> " ?\n)
  (insert "!> @author  "(shell-command-to-string "echo -n $(whoami)") ?\n)
  (insert "!> @date    "(shell-command-to-string "echo -n $(date +%Y-%m-%d)") ?\n)
  (insert "!> @brief   ???" ?\n)
  (insert "!> @details ???" ?\n)
  (insert "!> " ?\n)
  (insert "!-----------------------------------------------------------------------" ?\n)
  )

(defun alya-solver ()
  "Insert Elsest field."
  (interactive)
  (insert "ALGEBRAIC_SOLVER" ?\n)
  (insert "  SOLVER:           CG/GMRES/DEFLATED_CG, KRYLOV=10, COARSE: SPARSE/CHOLESKY" ?\n)
  (insert "  CONVERGENCE:      ITERATIONS=1000, TOLERANCE=1.0e-10 $ ADAPTIVE, RATIO=1.0e-2 " ?\n)
  (insert "  OUTPUT:           OFF                                $ CONVERGENCE/MARKET/NON_PRECONDITIONED_RESIDUAL" ?\n)
  (insert "  PRECONDITIONER:   DIAGONAL, RIGHT                    $ RAS/AII/BLOCK_DIAGONAL/LINELET, NEVER_CHANGE, THRESHOLD=1.0e-3" ?\n)
  (insert "  OPTIONS:          OFF                                $ FIXITY, REACTION_FORCE, FULL_ROW, SAVE_KRYLOV, SIZE=10, " ?\n)
  (insert "  $ OPENMP:         AUTO_TUNING/STATIC/GUIDED/DYNAMIC, CHUNK_SIZE=200, INTERFACE=0" ?\n)
  (insert "  $ RESIDUAL:       RHS                                $ USER/CONSTANT, VALUE=0.01" ?\n)
  (insert "  $ MATRIX_FORMAT:  CSR/ELL/COO" ?\n)
  (insert "  $ RHS_MINIMUM=    1.0e-14" ?\n)
  (insert "  $ COARSE_SOLVER:  ON" ?\n)
  (insert "  $ CONFIGURATION:  maphys-config.tx" ?\n)
  (insert "END_ALGEBRAIC_SOLVER" ?\n)
  )

(defun alya-elsest ()
  "Insert solver field."
  (interactive)
  (insert "ELSEST" ?\n)
  (insert "  STRATEGY:           BIN" ?\n)
  (insert "  NUMBER_BOXES=       50,50,50" ?\n)
  (insert "  $ STRATEGY=         OCT_TREE" ?\n)
  (insert "  $ MAXIMUM_ELEMENTS= 20" ?\n)
  (insert "  TOLERANCE=          1.0e-6" ?\n)
  (insert "END_ELSEST" ?\n)
  )

(defun alya-domain ()
  "Insert domain field."
  (interactive)
  (insert "$-------------------------------------------------------------" ?\n)
  (insert "DIMENSIONS" ?\n)
  (insert "  NODAL_POINTS=          xxx" ?\n)
  (insert "  ELEMENTS=              xxx" ?\n)
  (insert "  SPACE_DIMENSIONS=        2" ?\n)
  (insert "  TYPES_OF_ELEMENTS=     xxx" ?\n)
  (insert "  BOUNDARIES=            xxx" ?\n)
  (insert "  MATERIALS=               1" ?\n)
  (insert "  SUBDOMAINS=              1" ?\n)
  (insert "$  FIELDS:                  2" ?\n)
  (insert "$    FIELD = 1, DIMENSIONS=2, ELEMENTS" ?\n)
  (insert "$    FIELD = 2, DIMENSIONS=1, NODES" ?\n)
  (insert "$  END_FIELDS" ?\n)
  (insert "END_DIMENSIONS" ?\n)
  (insert "$-------------------------------------------------------------" ?\n)
  (insert "STRATEGY" ?\n)
  (insert "  INTEGRATION_RULE:                 Open " ?\n) 
  (insert "  DOMAIN_INTEGRATION_POINTS:        0" ?\n)
  (insert "  BOUNDARY_ELEMENT:                 ON" ?\n)
  (insert "  EXTRAPOLATE_BOUNDARY_CONDITIONS:  ON" ?\n)
  (insert "  GROUPS=                          100, SEQUENTIAL_FRONTAL $ PARALLEL_FRONTAL / SFC, FINE=xxx, COARSE=xxx / PARTITION / FIELD" ?\n)
  (insert "END_STRATEGY" ?\n)
  (insert "$-------------------------------------------------------------" ?\n)
  (insert "GEOMETRY" ?\n)
  (insert "  INCLUDE ./example.geo.dat" ?\n)
  (insert "END_GEOMETRY" ?\n)
  (insert "$-------------------------------------------------------------" ?\n)
  (insert "SETS" ?\n)
  (insert "  INCLUDE ./example.set.dat" ?\n)
  (insert "END_SETS" ?\n)
  (insert "$-------------------------------------------------------------" ?\n)
  (insert "BOUNDARY_CONDITIONS" ?\n)
  (insert "  INCLUDE ./example.fix.dat" ?\n)
  (insert "END_BOUNDARY_CONDITIONS" ?\n)
  (insert "$-------------------------------------------------------------" ?\n)
  (insert "FIELDS" ?\n)
  (insert "  INCLUDE ./example.fie.dat" ?\n)
  (insert "END_FIELDS" ?\n)
  (insert "$-------------------------------------------------------------" ?\n)
  )

(require 'easymenu)

(easy-menu-define my-menu global-map "My own menu"
  '("Alya"
    ["(Un)Comment Region" alya-comment-region t]
    "--"
    ["Insert Separator"   alya-separator t]
    ["Insert Separator Around Region"  alya-markup-region t]
    "--"
    ["Insert group definition"                       alya-doxygen-defgroup t]
    ["Insert Doxygen header for a subroutine"        alya-doxygen-subroutine t]
    ["Insert Doxygen header for a module"            alya-doxygen-module t]
    ["Insert Doxygen end of a module"                alya-doxygen-module-close t]
    ["Insert Doxygen header for a module subroutine" alya-doxygen-subroutine-module t]
    "--"
    ("Change Case"
     :help "Change the case of buffer"
     ["Upcase Region"   alya-upcase   t]
     ["Downcase Region" alya-downcase t]
     "--"
     ["Upcase Word"     upcase-word   t]
     ["Downcase Word"   downcase-word t]
     )
    "--" 
    ("Data"
     :help "Write a field"
     ["Parall"     alya-parall   t]
     ["MPI/IO"     alya-mpiio    t]
     ["Solver"     alya-solver   t]
     ["Elsest"     alya-elsest   t]
     ["Domain"     alya-domain   t]
     )    
    )
  )

(define-key-after (lookup-key global-map [menu-bar])
[Alya] ; shortcut for our menu
(cons "Alya" my-menu) 'tools) ; Our menu's name in cons.

;(easy-menu-add my-menu global-map)
