#########################################################################
#
#  tcl main file for the problem type Alya
#
########################################################################
proc InitGIDProject  { dir } {
    global Alya
    global Text
    global textfactor textfactors
    global FDPveloc FDPlengt FDPtempe FDPra FDPre FDPar FDPri
    global TurmuK TurmuO TurmuE TurmuMut TurmuReT TurmuMutMu

    AlyaGetPaths
    set Alya(ptypepath) $dir
    source  [file join $Alya(ptypepath) tcls      Graphs.tcl]
    source  [file join $Alya(ptypepath) tcls      Text_en.tcl]
    source  [file join $Alya(ptypepath) tcls      Colors.tcl]
    source  [file join $Alya(ptypepath) styles    Styles.tcl]
    AlyaText
    AlyaColor
    SetAlyaMenus

    # List of all windows
    set Alya(ListWindow) { .gid.postwin }
    
    # Constants for the Dimensionless Parameters
    set FDPveloc    1.0
    set FDPlengt    1.0
    set FDPtempe    1.0
    set textfactor  0
    set textfactors 0
    set TurmuK      1.0
    set TurmuO      1.0
    set TurmuE      1.0
    set TurmuMut    1.0
    set TurmuReT    1.0
    set TurmuMutMu  1.0
}

proc EndGIDProject { } {
    AlyaDestroyWindows
}

proc LoadGIDProject { filespd } {
    AlyaCheckModuleService
    AlyaDestroyWindows
}

proc InitGIDPostProcess { } {
    global Alya
    global NASTIN TEMPER CODIRE TURBUL EXMEDI NASTAL 
    global ALEFOR SOLIDZ GOTITA WAVEQU LEVELS
    global DODEME PARALL
    
    AlyaCheckModuleService
    if { $NASTIN == "On" } {AlyaNASTINPostprocess}
    if { $TEMPER == "On" } {AlyaTEMPERPostprocess}
    

}

proc BeforeMeshGeneration { meshsize } {
    AlyaBeforeMeshGeneration
}

proc AlyaCheckModuleService {  } {
    global NASTIN TEMPER CODIRE TURBUL EXMEDI NASTAL 
    global ALEFOR SOLIDZ GOTITA WAVEQU LEVELS
    global DODEME
    # Check if module or service is On 
    set Prbdata   [GiD_Info gendata]
    set NASTIN    [lindex $Prbdata [expr [lsearch $Prbdata "NASTIN_Module:#CB#(On,Off)"]+1 ] ]
    set TEMPER    [lindex $Prbdata [expr [lsearch $Prbdata "TEMPER_Module:#CB#(On,Off)"]+1 ] ]
    set CODIRE    [lindex $Prbdata [expr [lsearch $Prbdata "CODIRE_Module:#CB#(On,Off)"]+1 ] ]
    set TURBUL    [lindex $Prbdata [expr [lsearch $Prbdata "TURBUL_Module:#CB#(On,Off)"]+1 ] ]
    set EXMEDI    [lindex $Prbdata [expr [lsearch $Prbdata "EXMEDI_Module:#CB#(On,Off)"]+1 ] ]
    set NASTAL    [lindex $Prbdata [expr [lsearch $Prbdata "NASTAL_Module:#CB#(On,Off)"]+1 ] ]
    set ALEFOR    [lindex $Prbdata [expr [lsearch $Prbdata "ALEFOR_Module:#CB#(On,Off)"]+1 ] ]
    set SOLIDZ    [lindex $Prbdata [expr [lsearch $Prbdata "SOLIDZ_Module:#CB#(On,Off)"]+1 ] ]
    set GOTITA    [lindex $Prbdata [expr [lsearch $Prbdata "GOTITA_Module:#CB#(On,Off)"]+1 ] ]
    set WAVEQU    [lindex $Prbdata [expr [lsearch $Prbdata "WAVEQU_Module:#CB#(On,Off)"]+1 ] ]
    set LEVELS    [lindex $Prbdata [expr [lsearch $Prbdata "LEVELS_Module:#CB#(On,Off)"]+1 ] ]
    set DODEME    [lindex $Prbdata [expr [lsearch $Prbdata "DODEME_Service:#CB#(On,Off)"]+1 ] ]
}

proc AfterWriteCalcFileGIDProject { file error } {
    global Alya    
    global NASTIN TEMPER CODIRE TURBUL EXMEDI NASTAL 
    global ALEFOR SOLIDZ GOTITA WAVEQU LEVELS
    global DODEME 

    AlyaGetPaths
    set Basdirectory $Alya(projectpath).gid

    # Check if module or service is On 
    AlyaCheckModuleService

    # If problem is not solved, delete bas files
    if { $NASTIN != "On" } {
        set dat [file join $Basdirectory $Alya(projectname)-7.dat]
        set error [file delete $dat]
    }
    if { $TEMPER != "On" } {
        set dat [file join $Basdirectory $Alya(projectname)-8.dat]
        set error [file delete $dat]
    }
    if { $CODIRE != "On" } {
        set dat [file join $Basdirectory $Alya(projectname)-9.dat]
        set error [file delete $dat]
    }
    if { $TURBUL != "On" } {
        set dat [file join $Basdirectory $Alya(projectname)-10.dat]
        set error [file delete $dat]
    }
    if { $EXMEDI != "On" } {
        set dat [file join $Basdirectory $Alya(projectname)-11.dat]
        set error [file delete $dat]
    }
    if { $NASTAL != "On" } {
        set dat [file join $Basdirectory $Alya(projectname)-12.dat]
        set error [file delete $dat]
    }
    if { $ALEFOR != "On" } { 
        set dat [file join $Basdirectory $Alya(projectname)-13.dat]
        set error [file delete $dat]
    }
    if { $SOLIDZ != "On" } {
        set dat [file join $Basdirectory $Alya(projectname)-14.dat]
        set error [file delete $dat]
    }
    if { $GOTITA != "On" } {
        set dat [file join $Basdirectory $Alya(projectname)-15.dat]
        set error [file delete $dat]
    }
    if { $DODEME != "On" } {
        set dat [file join $Basdirectory $Alya(projectname)-16.dat]
        set error [file delete $dat]
    }
    if { $WAVEQU != "On" } {
        set dat [file join $Basdirectory $Alya(projectname)-17.dat]
        set error [file delete $dat]
    }
    if { $LEVELS != "On" } {
        set dat [file join $Basdirectory $Alya(projectname)-18.dat]
        set error [file delete $dat]
    }
}


proc EndGIDPostProcess { } {  
    global Alya
    #AlyaDestroyWindows
}

proc SetAlyaMenus { } { 
    CreateMenu       "Alya" "PREPOST"
    InsertMenuOption "Alya" "Update version"                                  0 "AlyaUpdate"            "PREPOST"
    InsertMenuOption "Alya" "---"                                             1 ""                      "PREPOST"
    InsertMenuOption "Alya" "Utilities"                                       2 ""                      "PRE"
    InsertMenuOption "Alya" "Utilities>Problem Data"                          3 ""                      "PRE"
    InsertMenuOption "Alya" "Utilities>Calculator"                            3 ""                      "PRE"
    InsertMenuOption "Alya" "Utilities>Problem Data>Air properties"           0 "AlyaAirProperties"     "PRE"
    InsertMenuOption "Alya" "Utilities>Problem Data>Blood properties"         1 "AlyaBloodProperties"   "PRE"
    InsertMenuOption "Alya" "Utilities>Problem Data>Water properties"         2 "AlyaWaterProperties"   "PRE"
    InsertMenuOption "Alya" "Utilities>Problem Data>All transient "           3 "AlyaTransient"         "PRE"
    InsertMenuOption "Alya" "Utilities>Problem Data>Transient to stationary " 4 "AlyaTransStat"         "PRE"
    InsertMenuOption "Alya" "Utilities>Problem Data>Stationary "              5 "AlyaStationary"        "PRE"
    InsertMenuOption "Alya" "Utilities>Calculator>Dimensionless parameters"   0 "AlyaComputeParameters" "PRE"
    InsertMenuOption "Alya" "Utilities>Calculator>Turbulent parameters"       1 "AlyaComputeTurmu"      "PRE"
    InsertMenuOption "Alya" "---"                                            10 ""                      "PRE"
    InsertMenuOption "Alya" "Output Files"                                   11 "AlyaFindOutputFiles"   "PREPOST"
    InsertMenuOption "Alya" "Live results"                                   12 ""                      "PREPOST"
    InsertMenuOption "Alya" "Live results>Global convergence"                13 "AlyaConvergence 0"     "PREPOST"
    InsertMenuOption "Alya" "Live results>NASTIN convergence"                14 "AlyaConvergence 1"     "PREPOST"
    InsertMenuOption "Alya" "Live results>TEMPER convergence"                15 "AlyaConvergence 2"     "PREPOST"
    InsertMenuOption "Alya" "Live results>CODIRE convergence"                16 "AlyaConvergence 3"     "PREPOST"
    InsertMenuOption "Alya" "Live results>TURBUL convergence"                17 "AlyaConvergence 4"     "PREPOST"
    InsertMenuOption "Alya" "Live results>EXMEDI convergence"                18 "AlyaConvergence 5"     "PREPOST"
    InsertMenuOption "Alya" "Live results>NASTAL convergence"                19 "AlyaConvergence 6"     "PREPOST"
    InsertMenuOption "Alya" "Live results>ALEFOR convergence"                20 "AlyaConvergence 7"     "PREPOST"
    InsertMenuOption "Alya" "Live results>SOLIDZ convergence"                20 "AlyaConvergence 10"    "PREPOST"
    InsertMenuOption "Alya" "Live results>GOTITA convergence"                15 "AlyaConvergence 11"    "PREPOST"
    InsertMenuOption "Alya" "Live results>WAVEQU convergence"                15 "AlyaConvergence 12"    "PREPOST"
    InsertMenuOption "Alya" "Live results>LEVELS convergence"                15 "AlyaConvergence 13"    "PREPOST"
    InsertMenuOption "Alya" "---"                                            21 ""                      "PREPOST"
    InsertMenuOption "Alya" "Alya About"                                     22 "AlyaAbout"             "PREPOST"
    UpdateMenus
}

proc AlyaGlobalConvergence { } {
    global AlyaAutomaticUpdate
    global whichplot
    set whichplot 0
    set AlyaAutomaticUpdate 1
    AlyaPostprocesWindow
}

proc AlyaConvergence { igraph } {
    global AlyaAutomaticUpdate
    global whichplot
    set whichplot $igraph
    set AlyaAutomaticUpdate 1
    AlyaPostprocesWindow
}

proc AlyaGetPaths { } {
    global Alya
    set Alya(projectpath) [GiD_Info Project]
    set Alya(projectpath) [lindex $Alya(projectpath) 1]
    set lname             [split  $Alya(projectpath) "\\ /"] 
    set Alya(projectname) [lindex $lname [expr {[llength $lname]-1}]]
    set Alya(projectpath) [join   $lname /]
    set Alya(resultspath) $Alya(projectpath).gid/results
    #set Alya(ptypepath)   [GiD_Info problemtypepath]
}

proc AlyaFindOutputFiles { } {
    global Alya
    global whichfile
    AlyaGetPaths
    set lfiles [glob -nocomplain "$Alya(resultspath)/*"]
    set nfiles [llength $lfiles]
    if { $nfiles == 0 } { 
        tk_messageBox -icon error  -type ok   \
            -message "There are no output files for this project"
        return
    } 
    for {set i 0} {$i < $nfiles} {incr i 1} {
        set ifile  [lindex $lfiles $i]
        set ifile  [string tolower $ifile]
        set lfiles [lreplace $lfiles $i $i $ifile]
    } 
    set lfiles [lsort -increasing $lfiles] 
    set W .selectfile
    toplevel $W 
    wm title $W " List of output files "
    wm geom  $W +200+200
    wm focusmodel $W active
    frame $W.top -borderwidth 2 -relief groove
    label $W.top.t -text " Select an output file "
    pack  $W.top.t -side top -anchor w
    pack  $W.top -side top -fill both -padx 2m -pady 1m
    set whichfile [lindex $lfiles 0]
    for {set i 0} {$i < $nfiles} {incr i 1} {
        set ifile [lindex $lfiles $i] 
        frame $W.top.$i
        radiobutton $W.top.$i.1 -text $ifile -variable whichfile -value $ifile 
        pack  $W.top.$i.1 -side top    -anchor w 
        pack  $W.top.$i   -side top -anchor w
    }
    frame  $W.buttons 
    button $W.buttons.1 -text " Accept " -width 14 -relief raised \
        -command "AlyaShowFile $W"
    pack   $W.buttons.1 -side left 
    button $W.buttons.2 -text " Cancel " -width 14 -relief raised \
        -command "destroy $W"
    pack   $W.buttons.2 -side left 
    pack   $W.buttons   -side top -fill both -padx 2m -pady 2m

}

proc AlyaShowFile { W } {
    global whichfile
    destroy $W
    set W .outputwindow
    toplevel $W
    wm geometry $W +100+100 
    wm focusmodel $W active
    wm title      $W " Output files window "
    frame $W.top -bd 2 -relief groove
    label $W.top.1 -text " File: $whichfile"
    pack  $W.top.1 -side top -anchor w
    pack  $W.top   -side top -fill both -padx 2m -pady 2m
    frame $W.screen -bd 5 -relief groove
    scrollbar $W.screen.sy -command "$W.screen.text yview" -orient vertical 
    text  $W.screen.text -bg white -height 15 -width 45   \
        -yscrollcommand [list $W.screen.sy set]  -wrap none -font "courier 8"
    pack $W.screen.sy -side right  -fill y
    pack $W.screen.text -side top -expand yes -fill both
    pack $W.screen -side top -expand yes -fill both 
    frame  $W.orders -bd 5 
    button $W.orders.c -text "  Close   " -height 1 -width 14 -relief raised \
        -command "destroy $W"
    pack  $W.orders.c -anchor s -side right
    button $W.orders.r -text "  Refresh   " -height 1 -width 14 -relief raised \
        -command "AlyaReadFile $W $whichfile"
    pack  $W.orders.r -anchor s -side right
    button $W.orders.o -text "  Other File   " -height 1 -width 14 -relief raised \
        -command "AlyaSelectOtherFile $W"
    pack  $W.orders.o -anchor s -side right
    pack  $W.orders   -anchor s -side bottom -ipadx 5 -ipady 5 

    AlyaReadFile $W $whichfile 
}

proc AlyaSelectOtherFile { W } { 
    destroy $W
    AlyaFindOutputFiles 
}

proc AlyaReadFile { W whichfile } {
    set luout [open $whichfile r 0600]
    $W.screen.text yview moveto 0.0
    $W.screen.text delete @0,0 end
    while { ![eof $luout] } {
        set linea [gets $luout]
        $W.screen.text insert end "$linea\n"
    }
    close $luout
}

proc AlyaAbout { } { 
    set w .about
    toplevel $w
    wm focusmodel $w active
    wm title    $w " Alya about "
    frame $w.top
    label $w.top.title -text " Alya: Computational Mechanics and Design" -font "courier 8"
    pack  $w.top.title -side top 
    pack  $w.top -padx 2m -pady 2m
    frame $w.information 
    label $w.information.version -text " Version 1.0 "
    pack  $w.information.version -side top -anchor w
    pack  $w.information  -side top -padx 2m -pady 2m
    frame $w.authors
    label $w.authors.authors -text " Program written by Guillaume Houzeaux"
    label $w.authors.adress  -text " Please, send any comment or suggestion to \n  \
                                 guillaume.houzeaux@upc.edu"
    pack $w.authors.authors -side top -anchor w
    pack $w.authors.adress  -side top -anchor w
    pack $w.authors -padx 2m -pady 2m
    frame  $w.bottom
    button $w.bottom.start -text "Accept" -height 1 -width 14 -command "destroy $w"
    pack   $w.bottom.start   -side left -anchor w
    pack $w.bottom -side bottom -padx 2m -pady 2m
}

proc AlyaEvalDirichletCondition { value x y z } {
    set result 0.0
    if { [catch { set result [format "%12.10e" [expr $value]] }]!=0} {
        WarnWinText "Error in Dirichlet boundary conditions!" 
    }
    return $result
}

proc AlyaEvalNeumannCondition { value x y z } {
    set result 0.0
    if { [catch { set result [format "%12.10e" [expr $value]] }]!=0} {
        WarnWinText "Error in Neumann boundary conditions!" 
    }
    return $result
}

proc AlyaProblemName { } {
    global Alya
    AlyaGetPaths
    set Alya(projectpath) [GiD_Info Project]
    set Alya(projectpath) [lindex $Alya(projectpath) 1]
    set lname             [split  $Alya(projectpath) "\\ /"] 
    set Alya(projectname) [lindex $lname [expr {[llength $lname]-1}]]
    return [file join $Alya(projectpath).gid data $Alya(projectname)]
}

proc AlyaProblemSolo { } {
    global Alya
    AlyaGetPaths
    set Alya(projectpath) [GiD_Info Project]
    set Alya(projectpath) [lindex $Alya(projectpath) 1]
    set lname             [split  $Alya(projectpath) "\\ /"] 
    set Alya(projectname) [lindex $lname [expr {[llength $lname]-1}]]
    return [file join $Alya(projectname)]
}

proc AlyaUpdate { } {
    set Alya(ptypepath)   [GiD_Info problemtypepath]
    GiD_Process Mescape Data Defaults TransfProblem $Alya(ptypepath)   
}

proc AlyaCreateGeneralWindow { w title x y wbg } {
    global Color 
    global Alya

    AlyaGetPaths

    catch         {destroy $w }
    toplevel      $w
    wm geometry   $w   +$x+$y
    wm resizable  $w   0 0 
    wm transient  $w   $wbg
    wm title      $w   $title
    #wm iconbitmap $w   [file join $Alya(ptypepath) pics Favent.ico]
    focus         $w
    #background of the window    
    set  Color(background) [$w cget -background]
    #border of the window
    set  Color(border)     [AlyaCCColorActivo $Color(background)]
    #boton importante
    set  Color(boton1)     $Color(background)
    set  Color(boton1fg)   $Color(background)
    #boton general           
    set  Color(boton2)     $Color(background)

    bind $w   <Control-w> "ShowWizard"
    #bind $w   <Alt-w>     "FaventHideWindows"
    #bind .gid <Alt-w>     "FaventHideWindows"
}

proc AlyaNumberDecode { linea } {
    global number

    # Initializations
    set ndime 0

    # Look for comments
    set first [string index $linea 0]
    if {$first == "$"} {
        set ndime -1
        return $ndime
    } elseif {$first == "#"} {
        set ndime -1
        return $ndime
    }
    # Decodes linea
    set llinea [split   $linea " "]
    set len    [llength $llinea]
    if {$len != "0"} {   
        set ndime 0
        for {set j 0} {$j < $len} {incr j 1} {
            if {[lindex $llinea $j] != "" } {
                incr ndime 1
                set number($ndime) [lindex $llinea $j]
                if { $number($ndime)== "NaN" } { set number($ndime) 0 }
            }
        }
    } else {
        set ndime -1
    }
    return $ndime
}

proc AlyaDestroyWindows { } {
    global Alya
    foreach namewindow $Alya(ListWindow) {
        if [winfo exists $namewindow] {
            destroy $namewindow
        }
    }      
}

proc AlyaAirProperties { } {
    # Put air properties at 20 degrees

    # NASTIN
    #kike: it is possible to modify general data question values without process, with the GiD_AccessValue command
    #e.g.
    #GiD_Process Mescape Data ProblemData -SingleField- NSI_Viscosity                 1.8205e-05    
    GiD_AccessValue set gendata NSI_Viscosity 1.8205e-05
    
    GiD_Process Mescape Data ProblemData -SingleField- NSI_Density                   1.2047
    GiD_Process Mescape Data ProblemData -SingleField- NSI_Beta                      3.4112e-03
    GiD_Process Mescape Data ProblemData -SingleField- NSI_Modulus_of_G              9.81
    # TEMPER
    GiD_Process Mescape Data ProblemData -SingleField- TEM_Viscosity                 1.8205e-05
    GiD_Process Mescape Data ProblemData -SingleField- TEM_Density                   1.2047
    GiD_Process Mescape Data ProblemData -SingleField- TEM_Specific_heat             1.0061e+03
    GiD_Process Mescape Data ProblemData -SingleField- TEM_Thermal_conductivity      2.5596e-02
    GiD_Process Mescape Data ProblemData -SingleField- TEM_Turbulent_Prandtl_number  0.9
    # TURBUL
    GiD_Process Mescape Data ProblemData -SingleField- TUR_Viscosity                 1.8205e-05
    GiD_Process Mescape Data ProblemData -SingleField- TUR_Density                   1.2047
    GiD_Process Mescape Data ProblemData -SingleField- TUR_Beta                      3.4112e-03
    GiD_Process Mescape Data ProblemData -SingleField- TUR_Norm_gravity              9.81
    GiD_Process Mescape Data ProblemData -SingleField- TUR_Turbulent_Prandtl_number  0.9
}

proc AlyaBloodProperties { } {
    # Put blood properties

    # NASTIN
    GiD_Process Mescape Data ProblemData -SingleField- NSI_Viscosity                 0.0036
    GiD_Process Mescape Data ProblemData -SingleField- NSI_Density                   994.0
}

proc AlyaWaterProperties { } {
    # Put air properties at 20 degrees

    # NASTIN
    GiD_Process Mescape Data ProblemData -SingleField- NSI_Viscosity                 9.7720e-4
    GiD_Process Mescape Data ProblemData -SingleField- NSI_Density                   9.9778e+2
    GiD_Process Mescape Data ProblemData -SingleField- NSI_Beta                      3.4112E-3
    GiD_Process Mescape Data ProblemData -SingleField- NSI_Modulus_of_G              9.81
    # TEMPER
    GiD_Process Mescape Data ProblemData -SingleField- TEM_Viscosity                 9.7720e-4
    GiD_Process Mescape Data ProblemData -SingleField- TEM_Density                   9.9778e+2
    GiD_Process Mescape Data ProblemData -SingleField- TEM_Specific_heat             4.0764e+3
    GiD_Process Mescape Data ProblemData -SingleField- TEM_Thermal_conductivity      0.60475
    GiD_Process Mescape Data ProblemData -SingleField- TEM_Turbulent_Prandtl_number  0.9
    # TURBUL
    GiD_Process Mescape Data ProblemData -SingleField- TUR_Viscosity                 9.7720e-4
    GiD_Process Mescape Data ProblemData -SingleField- TUR_Density                   9.9778e+2 
    GiD_Process Mescape Data ProblemData -SingleField- TUR_Beta                      3.4112E-3
    GiD_Process Mescape Data ProblemData -SingleField- TUR_Norm_gravity              9.81
    GiD_Process Mescape Data ProblemData -SingleField- TUR_Turbulent_Prandtl_number  0.9
}


proc AlyaTransient { } {
    global tstep tfina tinit mitgl
    global AskValueName  AskValueTitle  AskValueQuestion AskValueNValues AskValueTextValue \
           AskValueValue AskValueButton AskValueCommand 

    set Prbdata   [GiD_Info gendata]
    set tinit     [lindex $Prbdata [expr [lsearch $Prbdata "Time_initial="]+1 ] ]
    set tfina     [lindex $Prbdata [expr [lsearch $Prbdata "Time_final="]+1 ] ]
    set tstep     [lindex $Prbdata [expr [lsearch $Prbdata "Time_step_size="]+1 ] ]
    set mitgl     [lindex $Prbdata [expr [lsearch $Prbdata "Maximum_number_of_steps="]+1 ] ]       
    
    set AskValueNValues      4
    set AskValueName         transient
    set AskValueTitle        "Enter values"
    set AskValueQuestion     "Time integration parameters"
    set AskValueTextValue(1) "Initial time"
    set AskValueTextValue(2) "Final time"
    set AskValueTextValue(3) "Time step"
    set AskValueTextValue(4) "Maximum number of time steps"

    set AskValueValue(1)     $tinit
    set AskValueValue(2)     $tfina
    set AskValueValue(3)     $tstep
    set AskValueValue(4)     $mitgl
    set AskValueButton(1)    "Ok"
    set AskValueButton(2)    "Cancel"
    set AskValueCommand(1)   "AlyaTransientValue"
    set AskValueCommand(2)   ""

    AlyaAskValues

}

proc AlyaTransientValue { } {
    global tinit tfina tstep mitgl
    global AskValueValue 

    set tinit $AskValueValue(1)
    set tfina $AskValueValue(2)
    set tstep $AskValueValue(3)
    set mitgl $AskValueValue(4)

    # Put all transient
    GiD_Process Mescape Data ProblemData -SingleField- Time_coupling                 Global_prescribed
    GiD_Process Mescape Data ProblemData -SingleField- Time_step_size                $tstep
    GiD_Process Mescape Data ProblemData -SingleField- Time_initial                  $tinit
    GiD_Process Mescape Data ProblemData -SingleField- Time_final                    $tfina
    GiD_Process Mescape Data ProblemData -SingleField- Maximum_number_of_steps       $mitgl

    GiD_Process Mescape Data ProblemData -SingleField- NSI_Temporal_derivatives      On
    GiD_Process Mescape Data ProblemData -SingleField- NSI_Temporal_weighting        All

    GiD_Process Mescape Data ProblemData -SingleField- TEM_Temporal_derivatives      On
    GiD_Process Mescape Data ProblemData -SingleField- TEM_Temporal_weighting        All

    GiD_Process Mescape Data ProblemData -SingleField- TUR_Temporal_derivatives      On
    GiD_Process Mescape Data ProblemData -SingleField- TUR_Temporal_weighting        All
}

proc AlyaTransStat { } {
    global safet mitgl
    global AskValueName  AskValueTitle  AskValueQuestion AskValueNValues AskValueTextValue \
           AskValueValue AskValueButton AskValueCommand 

    set Prbdata   [GiD_Info gendata]
    set safet     [lindex $Prbdata [expr [lsearch $Prbdata "NSI_Safety_factor="]+1 ] ]
    set mitgl     [lindex $Prbdata [expr [lsearch $Prbdata "Maximum_number_of_steps="]+1 ] ]       


    set AskValueNValues      2
    set AskValueName         transient
    set AskValueTitle        "Enter values"
    set AskValueQuestion     "Transient-to-stationary method parameters"
    set AskValueTextValue(1) "Safety factor"
    set AskValueTextValue(2) "Maximum number of time steps"

    set AskValueValue(1)     $safet
    set AskValueValue(2)     $mitgl
    set AskValueButton(1)    "Ok"
    set AskValueButton(2)    "Cancel"
    set AskValueCommand(1)   "AlyaTransStatValue"
    set AskValueCommand(2)   ""

    AlyaAskValues

}

proc AlyaTransStatValue { } {
    global safet mitgl
    global AskValueValue 

    set safet $AskValueValue(1)
    set mitgl $AskValueValue(2)

    # Put all transient
    GiD_Process Mescape Data ProblemData -SingleField- Time_coupling                    Global_from_critical
    GiD_Process Mescape Data ProblemData -SingleField- Time_initial                     0.0
    GiD_Process Mescape Data ProblemData -SingleField- Time_final                       1.0e8
    GiD_Process Mescape Data ProblemData -SingleField- Maximum_number_of_steps          $mitgl
    GiD_Process Mescape Data ProblemData -SingleField- Maximum_coupling_iterations      1
    GiD_Process Mescape Data ProblemData -SingleField- Maximum_block_global_iterations  1,1,1


    GiD_Process Mescape Data ProblemData -SingleField- NSI_Temporal_derivatives         On
    GiD_Process Mescape Data ProblemData -SingleField- NSI_Safety_factor                $safet
    GiD_Process Mescape Data ProblemData -SingleField- NSI_Temporal_weighting           Galerkin

    GiD_Process Mescape Data ProblemData -SingleField- TEM_Temporal_derivatives         On
    GiD_Process Mescape Data ProblemData -SingleField- TEM_Safety_factor                $safet
    GiD_Process Mescape Data ProblemData -SingleField- TEM_Temporal_weighting           Galerkin

    GiD_Process Mescape Data ProblemData -SingleField- TUR_Temporal_derivatives         On
    GiD_Process Mescape Data ProblemData -SingleField- TUR_Safety_factor                $safet
    GiD_Process Mescape Data ProblemData -SingleField- TUR_Temporal_weighting          Galerkin

}


proc AlyaStationary { } {
    global tstep tfina tinit mitgl
    global AskValueName  AskValueTitle  AskValueQuestion AskValueNValues AskValueTextValue \
           AskValueValue AskValueButton AskValueCommand 

    set Prbdata   [GiD_Info gendata]
    set micou     [lindex $Prbdata [expr [lsearch $Prbdata "Maximum_coupling_iterations="]+1 ] ]    
    
    set AskValueNValues      1
    set AskValueName         stationary
    set AskValueTitle        "Enter values"
    set AskValueQuestion     "Stationary parameters"
    set AskValueTextValue(1) "Maximum number of global iterations"

    set AskValueValue(1)     $micou
    set AskValueButton(1)    "Ok"
    set AskValueButton(2)    "Cancel"
    set AskValueCommand(1)   "AlyaStationaryValue"
    set AskValueCommand(2)   ""

    AlyaAskValues 

}

proc AlyaStationaryValue { } {
    global tinit tfina tstep mitgl
    global AskValueValue 

    set micou $AskValueValue(1)

    # Put all stationary
    GiD_Process Mescape Data ProblemData -SingleField- Maximum_coupling_iterations   $micou
    GiD_Process Mescape Data ProblemData -SingleField- Time_coupling                 Off
    GiD_Process Mescape Data ProblemData -SingleField- NSI_Temporal_derivatives      Off
    GiD_Process Mescape Data ProblemData -SingleField- TEM_Temporal_derivatives      Off
    GiD_Process Mescape Data ProblemData -SingleField- TUR_Temporal_derivatives      Off
}




proc AlyaAskValues { } {
    global AskValueName  AskValueTitle  AskValueQuestion AskValueNValues AskValueTextValue \
           AskValueValue AskValueButton AskValueCommand 
    global Text Color


    # Widget geometry
    set      w .gid.askvalues$AskValueName
    AlyaCreateGeneralWindow $w $AskValueTitle 400 140 .gid

    frame  $w.tout      -bg $Color(border)
    frame  $w.tout.all  -bg $Color(border)

    frame  $w.tout.all.top   -bg $Color(border)
    # Question
    label  $w.tout.all.top.info  -text $AskValueQuestion -bg $Color(border) -font {-weight bold -size 8} 
    pack   $w.tout.all.top.info  -side top -anchor w -pady 2m

    frame $w.tout.all.top.center -bg $Color(background) -borderwidth 3 -relief groove
    # Loop over values
    for {set ivalu 1} {$ivalu <= $AskValueNValues } {incr ivalu 1} {
        set            jvalu [expr $ivalu-1]
        label          $w.tout.all.top.center.vname$ivalu  -text $AskValueTextValue($ivalu)  -bg $Color(background)
        grid   config  $w.tout.all.top.center.vname$ivalu  -row $jvalu -column 0 -sticky w 
        entry          $w.tout.all.top.center.entry$ivalu  -width 10 -relief sunken -textvariable AskValueValue($ivalu) -bg white
        grid   config  $w.tout.all.top.center.entry$ivalu  -row $jvalu -column 1 -sticky w 
    }
    pack   $w.tout.all.top.center   -side top  -anchor center -ipady 4m -fill x

    pack   $w.tout.all.top -padx 1m -pady 1m

    pack   $w.tout.all       -side top  

    # Buttons
    frame  $w.tout.bot -bg $Color(border)
    button $w.tout.bot.apply -text $AskValueButton(1) -height 1 -width 10  \
           -command "AlyaAskValueOk $AskValueCommand(1)"
    pack   $w.tout.bot.apply -side left -padx 1m 
    button $w.tout.bot.close -text $AskValueButton(2) -height 1 -width 10  -command "destroy $w"
    pack   $w.tout.bot.close -side left -padx 1m 
    pack   $w.tout.bot       -side top  -pady 2m
   
    pack   $w.tout 

    bind   $w <Return>     "$w.tout.bot.apply invoke"
    bind   $w <Escape>     "$w.tout.bot.close invoke"
}

proc AlyaAskValueOk { MyCommand1 } {
    # Check if entries are numbers
    global AskValueName  AskValueTitle  AskValueQuestion AskValueNValues AskValueTextValue \
           AskValueValue AskValueButton AskValueCommand 
    global Text
    for {set ivalu 1} {$ivalu <= $AskValueNValues } {incr ivalu 1} {

        if { [catch { set result [expr $AskValueValue($ivalu)] }]!=0} {
            AlyaError [concat "A value is not a number: " $AskValueValue($ivalu)]
            return -1
        }
    }
    catch {destroy .gid.askvalues$AskValueName}
    $MyCommand1
}

proc AlyaGetNdime { } {
    set Ndime 2
    foreach layer [GiD_Info layers] { 
        set bbox [lindex [eval GiD_Info layers -bbox -use geometry  [GiD_Info layers]] 0]
        set zmin [lindex $bbox 2]
        set zmax [lindex $bbox 5]
        if { [expr (abs($zmin)+abs($zmax))>1e-6] } {
            set Ndime 3
            break
        }
    }
    return $Ndime
}


proc AlyaComputeParameters { } {
    AlyaCalculateDimensionlessParameters
    AlyaAskParameters
}

proc AlyaComputeTurmu { } {
    AlyaCalculateTurmu
    AlyaAskParametersTurmu
}

proc AlyaAskParameters { } {
    global FDPveloc FDPlengt FDPtempe FDPra FDPre FDPar FDPri FDPgr
    global Text Color Alya

    set   w .gid.faventparameters
    AlyaCreateGeneralWindow $w $Text(tit013) 400 140 .gid

    frame         $w.tout                     -bg $Color(border)     
    frame         $w.tout.wall                -bg $Color(border)      

    label         $w.tout.wall.title          -text $Text(lab278)\n$Text(lab279) -bg $Color(border) -font {-weight bold -size 8} -justify left
    pack          $w.tout.wall.title          -side top -anchor w  -pady 1m 

    frame         $w.tout.wall.up             -borderwidth 3 -relief groove     
    # Flow scales
    frame         $w.tout.wall.up.left 
    label         $w.tout.wall.up.left.veloca -text [concat $Text(lab280) " U="] 
    grid   config $w.tout.wall.up.left.veloca -row 0 -column 0 -sticky w  
    entry         $w.tout.wall.up.left.velocb -relief sunken -width 10 -textvariable FDPveloc
    grid   config $w.tout.wall.up.left.velocb -row 0 -column 1 -sticky w  
    label         $w.tout.wall.up.left.lengta -text [concat $Text(lab281) " L="]
    grid   config $w.tout.wall.up.left.lengta -row 1 -column 0 -sticky w  
    entry         $w.tout.wall.up.left.lengtb -relief sunken -width 10 -textvariable FDPlengt
    grid   config $w.tout.wall.up.left.lengtb -row 1 -column 1 -sticky w   
    label         $w.tout.wall.up.left.tempea -text [concat $Text(lab282) " T="] 
    grid   config $w.tout.wall.up.left.tempea -row 2 -column 0 -sticky w  
    entry         $w.tout.wall.up.left.tempeb -relief sunken -width 10 -textvariable FDPtempe 
    grid   config $w.tout.wall.up.left.tempeb -row 2 -column 1 -sticky w  
    pack          $w.tout.wall.up.left        -side left -anchor w -padx 2m
    # Calculate button
    frame         $w.tout.wall.up.center      

    image  create photo imarrow  -file [file join $Alya(ptypepath)/pics/favent-arrow.gif] 
    #canvas $w.tout.wall.up.center.canvas -width 33 -height 33 -bd 0 
    #$w.tout.wall.up.center.canvas create image 2 2 -image imarrow -anchor nw
    #pack   $w.tout.wall.up.center.canvas -side top -anchor w -padx 2m -pady 2m
 
    button        $w.tout.wall.up.center.text  -image imarrow \
                  -height 21 -width 21 -relief raised -command "AlyaUpdateParameters $w"
    pack          $w.tout.wall.up.center.text  -side top

    pack          $w.tout.wall.up.center       -side left -anchor w -padx 2m
    # Dimensionless parameters
    frame         $w.tout.wall.up.righ        
    label         $w.tout.wall.up.righ.re1    -text  [concat $Text(lab283) " Re="] -font {-weight bold -size 8}
    grid   config $w.tout.wall.up.righ.re1    -row 0 -column 0 -sticky w  
    label         $w.tout.wall.up.righ.re2    -text  $FDPre 
    grid   config $w.tout.wall.up.righ.re2    -row 0 -column 1 -sticky w
    label         $w.tout.wall.up.righ.ra1    -text  [concat $Text(lab284) " Ra="] -font {-weight bold -size 8}
    grid   config $w.tout.wall.up.righ.ra1    -row 1 -column 0  -sticky w
    label         $w.tout.wall.up.righ.ra2    -text  $FDPra 
    grid   config $w.tout.wall.up.righ.ra2    -row 1 -column 1 -sticky w
    label         $w.tout.wall.up.righ.ar1    -text  [concat $Text(lab285) " Gr="] -font {-weight bold -size 8} 
    grid   config $w.tout.wall.up.righ.ar1    -row 2 -column 0 -sticky w
    label         $w.tout.wall.up.righ.ar2    -text  $FDPgr 
    grid   config $w.tout.wall.up.righ.ar2    -row 2 -column 1 -sticky w
    label         $w.tout.wall.up.righ.ri1    -text  [concat $Text(lab286) " Ri="] -font {-weight bold -size 8} 
    grid   config $w.tout.wall.up.righ.ri1    -row 3 -column 0 -sticky w
    label         $w.tout.wall.up.righ.ri2    -text  $FDPri 
    grid   config $w.tout.wall.up.righ.ri2    -row 3 -column 1 -sticky w
    pack          $w.tout.wall.up.righ        -side left -anchor w  -padx 2m

    pack          $w.tout.wall.up             -side top -anchor w -ipady 2m 
    pack          $w.tout.wall                -padx 2m -pady 2m 
    pack          $w.tout 

    # Close button
    frame         $w.buttons        -bd 5 -bg [AlyaCCColorActivo [$w cget -background]]
    button        $w.buttons.close  -text $Text(bot002) -height 1 -width 10 -relief raised -command "destroy $w" 
    pack          $w.buttons.close  -pady 1m
    pack          $w.buttons        -side top -fill both -expand yes


    bind         $w <Return> "tkButtonInvoke $w.tout.wall.up.center.text"
    bind         $w <Escape> "destroy $w"

}

proc AlyaCheckNumber { a } {
    global Text
    set patron {^((-)?[0-9]+(.[0-9]+)?)$}
    if { [regexp $patron $a] != 1 } {
        AlyaError $Text(lab111)
        return -1
    }
    return 0
}

proc AlyaUpdateParameters { w } {
    global FDPre FDPra FDPgr FDPri
    AlyaCalculateDimensionlessParameters
    $w.tout.wall.up.righ.re2 conf -text $FDPre 
    $w.tout.wall.up.righ.ra2 conf -text $FDPra 
    $w.tout.wall.up.righ.ar2 conf -text $FDPgr 
    $w.tout.wall.up.righ.ri2 conf -text $FDPri 

}

proc AlyaCalculateDimensionlessParameters { } {
    global FDPveloc FDPlengt FDPtempe FDPra FDPre FDPar FDPri FDPgr
    global Text

    if { [catch {set results [expr $FDPveloc]}] != 0 } { 
        AlyaError [= "The velocity scale must be a real number"] 
        return
    }
    if { $FDPveloc < 0.0 } {
        AlyaError $Text(lab116)
        return
    }

    if { [catch {set results [expr $FDPlengt]}] != 0 } { 
        AlyaError [= "The length scale must be a real number"] 
        return
    } 
    if { $FDPlengt < 0.0 } {
        AlyaError $Text(lab117)
        return
    }

    if { [catch {set results [expr $FDPtempe]}] != 0 } { 
        AlyaError [= "The temperature scale must be a real number"] 
        return
    } 
    if { $FDPtempe < 0.0 } {
        AlyaError $Text(lab118)
        return
    }

    set Prbdata             [GiD_Info gendata]
    set PrbdataDensity      [lindex $Prbdata [expr [lsearch $Prbdata "NSI_Density="]+1 ] ]
    set PrbdataViscosity    [lindex $Prbdata [expr [lsearch $Prbdata "NSI_Viscosity="]+1 ] ]
    set PrbdataSpecificHeat [lindex $Prbdata [expr [lsearch $Prbdata "TEM_Specific_heat="]+1 ] ]
    set PrbdataConductivity [lindex $Prbdata [expr [lsearch $Prbdata "TEM_Thermal_conductivity="]+1 ] ]
    set PrbdataVolExpansion [lindex $Prbdata [expr [lsearch $Prbdata "NSI_Beta="]+1 ] ]
    set PrbdataBousGravity  [lindex $Prbdata [expr [lsearch $Prbdata "NSI_Modulus_of_G="]+1 ] ]

    set FDPgr [expr $PrbdataVolExpansion*$PrbdataBousGravity*double($FDPtempe)*double($FDPlengt)*double($FDPlengt)*double($FDPlengt)*$PrbdataDensity*$PrbdataDensity/($PrbdataViscosity*$PrbdataViscosity)]
    if { $FDPveloc > 0.0000001 } {
        set FDPre [expr double($FDPveloc)*double($FDPlengt)*$PrbdataDensity/double($PrbdataViscosity)]
        set FDPri [expr $PrbdataVolExpansion*$PrbdataBousGravity*double($FDPlengt)*double($FDPtempe)/(double($FDPveloc)*double($FDPveloc))]
    } else {
        set Uscal [expr sqrt($PrbdataVolExpansion*$PrbdataBousGravity*double($FDPlengt)*double($FDPtempe))]
        set FDPre [expr $PrbdataDensity*double($Uscal)*double($FDPlengt)/double($PrbdataViscosity)]
        set FDPri [expr 1.0]
    }
    set FDPre [format "%5.3E" $FDPre] 
    set FDPgr [format "%5.3E" $FDPgr] 
    set FDPri [format "%5.3E" $FDPri] 

    set FDPra [expr $PrbdataVolExpansion*$PrbdataSpecificHeat*$PrbdataBousGravity*$PrbdataDensity*$PrbdataDensity*double($FDPlengt)*double($FDPlengt)*double($FDPlengt)*$FDPtempe/double($PrbdataViscosity*$PrbdataConductivity)]
    set FDPra [format "%5.3E" $FDPra] 
    #AlyaAskParameters
}

proc AlyaError { message } {
    global Color Text Alya
    set w .gid.faventerror
    AlyaCreateGeneralWindow $w "Alya - Error" 400 140 .gid
    frame  $w.all             -bg   $Color(border)
    # Question
    frame  $w.all.info        -bg   $Color(background)
    image  create photo imerror  -file [file join $Alya(ptypepath)/pics/alya-error.gif] 
    frame  $w.all.info.titl0        
    canvas $w.all.info.titl0.canvas -width 33 -height 33 -bd 0 
    $w.all.info.titl0.canvas create image 2 2 -image imerror -anchor nw
    pack   $w.all.info.titl0.canvas -side top -anchor w -padx 2m -pady 2m
    pack   $w.all.info.titl0        -side left
    label  $w.all.info.titl1  -text $message  -font {-weight bold -size 8} -justify left
    pack   $w.all.info.titl1  -side left -anchor w -pady 1m -padx 2m
    pack   $w.all.info        -side top -anchor w 
    # Buttons
    frame  $w.all.button   -bg   $Color(border)
    button $w.all.button.1 -text $Text(bot003) -height 1 -width 10 -command "destroy $w"
    pack   $w.all.button.1 -side left    -padx 1m 
    pack   $w.all.button   -side top     -pady 2m
    pack   $w.all
    bind   $w <Return> "tkButtonInvoke $w.all.button.1"
}


proc AlyaPostFillVELOC { } {
    global Text
    catch {destroy .gid.faventfactor} 
    set rlist [GiD_Info post get cur_results_list Contour_Fill] 
    set Iresu [lsearch $rlist "VELOC"] 
    if { $Iresu == "-1" } {
        AlyaError $Text(lab100)
        return
    }
    GiD_Process Mescape utilities Comments Automatic No escape
    GiD_Process Mescape utilities render Normal
    GiD_Process Mescape Results ContOptions NumberOfColors 20
    GiD_Process Mescape Results ContOptions SetMinOptions ResetValue
    GiD_Process Mescape Results ContOptions SetMaxOptions ResetValue
    GiD_Process Mescape DisplayStyle BodyBound
    GiD_Process Mescape Results Contourfill VELOC |VELOC|
    GiD_Process Mescape Results ContOptions ChangeTitleLeg $Text(lab034) escape
    GiD_Process Mescape
}
proc AlyaPostVectVELOC { } {
    global Text  textfactor

    set rlist [GiD_Info post get cur_results_list Display_Vectors] 
    set Iresu [lsearch $rlist "VELOC"] 
    if { $Iresu == "-1" } {
        AlyaError $Text(lab100)
        return
    }
    # Find factor
    set res_anal [GiD_Info post get cur_analisis]
    set res_step [GiD_Info post get cur_step]
    set res_view "Display_Vectors"
    set res_res  "VELOC"
    set textfactor [GiD_Info post get cur_vector_factor $res_view $res_res $res_anal $res_step]
    AlyaPostFacto
    AlyaPostVeloV
}

proc AlyaPostVeloV { } {
    global textfactor Text

    set factor [expr double($textfactor)]
    GiD_Process Mescape utilities Comments Automatic No escape
    GiD_Process Mescape utilities render FlatRender
    GiD_Process Mescape results options vectorcolour colourmodules
    GiD_Process Mescape Results DisplayVectors VELOC |VELOC| $factor
    GiD_Process $factor
    GiD_Process Mescape Results ContOptions ChangeTitleLeg $Text(lab034) escape   
    GiD_Process Mescape DisplayStyle Boundaries escape
    GiD_Process Mescape
} 

proc AlyaPostFillTEMPE { } {
    global Text

    catch {destroy .gid.faventfactor}    
    set rlist [GiD_Info post get cur_results_list Contour_Fill] 
    set Iresu [lsearch $rlist "TEMPE"] 
    if { $Iresu == "-1" } {
        AlyaError $Text(lab100)
        return
    }
    GiD_Process Mescape utilities Comments Automatic No escape
    GiD_Process Mescape utilities render Normal
    GiD_Process Mescape Results ContOptions NumberOfColors 20
    GiD_Process Mescape Results ContOptions SetMinOptions ResetValue
    GiD_Process Mescape Results ContOptions SetMaxOptions ResetValue
    GiD_Process Mescape DisplayStyle BodyBound
    GiD_Process Mescape Results Contourfill TEMPE
    GiD_Process Mescape Results ContOptions ChangeTitleLeg $Text(lab036) escape
    GiD_Process Mescape
}

proc AlyaPostLineSTREA { } {
    global Text
    catch {destroy .gid.faventfactor}    
    set rlist [GiD_Info post get cur_results_list Contour_Fill] 
    set Iresu [lsearch $rlist "STREA"] 
    if { $Iresu == "-1" } {
        AlyaError $Text(lab100)
        return
    }
    GiD_Process Mescape utilities Comments Automatic No escape
    GiD_Process Mescape utilities render Normal
    GiD_Process Mescape Results ContOptions NumberOfColors 20
    GiD_Process Mescape Results ContOptions SetMinOptions ResetValue
    GiD_Process Mescape Results ContOptions SetMaxOptions ResetValue
    GiD_Process Mescape DisplayStyle BodyBound
    GiD_Process Mescape Results Contourline STREA
    GiD_Process Mescape Results ContOptions ChangeTitleLeg [= "STREA"] escape
    GiD_Process Mescape
}

proc AlyaPostFillPRESS { } {
    global Text
    catch {destroy .gid.faventfactor}    
    set rlist [GiD_Info post get cur_results_list Contour_Fill] 
    set Iresu [lsearch $rlist "PRESS"] 
    if { $Iresu == "-1" } {
        AlyaError $Text(lab100)
        return
    }
    GiD_Process Mescape utilities Comments Automatic No escape
    GiD_Process Mescape utilities render Normal
    GiD_Process Mescape Results ContOptions NumberOfColors 20
    GiD_Process Mescape Results ContOptions SetMinOptions ResetValue
    GiD_Process Mescape Results ContOptions SetMaxOptions ResetValue
    GiD_Process Mescape DisplayStyle BodyBound
    GiD_Process Mescape Results Contourfill PRESS
    GiD_Process Mescape Results ContOptions ChangeTitleLeg $Text(lab035) escape
    GiD_Process Mescape
} 

proc AlyaPostFillMACHN { } {
    global Text
    catch {destroy .gid.faventfactor}    
    set rlist [GiD_Info post get cur_results_list Contour_Fill] 
    set Iresu [lsearch $rlist "MACNH"] 
    if { $Iresu == "-1" } {
        AlyaError $Text(lab100)
        return
    }
    GiD_Process Mescape utilities Comments Automatic No escape
    GiD_Process Mescape utilities render Normal
    GiD_Process Mescape Results ContOptions NumberOfColors 20
    GiD_Process Mescape Results ContOptions SetMinOptions ResetValue
    GiD_Process Mescape Results ContOptions SetMaxOptions ResetValue
    GiD_Process Mescape DisplayStyle BodyBound
    GiD_Process Mescape Results Contourfill MACHN
    GiD_Process Mescape Results ContOptions ChangeTitleLeg $Text(lab035) escape
    GiD_Process Mescape 
}

proc AlyaPostSurfTEMPE { } {
    global Text  textfactors

    set rlist [GiD_Info post get cur_results_list Result_Surface] 
    set Iresu [lsearch $rlist "TEMPE"] 
    if { $Iresu == "-1" } {
        AlyaError $Text(lab100)
        return
    }
    # Find factor
    set res_anal     [GiD_Info post get cur_analisis]
    set res_step     [GiD_Info post get cur_step]
    set res_view     "Result_surface"
    set res_res      "TEMPE"
    set res_comp     "TEMPE"
    if { $textfactors == "0" } {
        set textfactors  [GiD_Info post get cur_result_surface_factor \
                              $res_view $res_res $res_comp $res_anal $res_step ]
    } 
    AlyaPostFactos
    AlyaPostTEMPES
}

proc AlyaPostTEMPES { } {
   global textfactors Text

    set factors [expr double($textfactors)]
    GiD_Process Mescape utilities Comments Automatic No escape
    GiD_Process Mescape utilities render Smooth
    GiD_Process Mescape Results ContOptions NumberOfColors 20
    GiD_Process Mescape Results ContOptions SetMinOptions ResetValue
    GiD_Process Mescape Results ContOptions SetMaxOptions ResetValue
    GiD_Process Mescape DisplayStyle BodyBound
    GiD_Process Mescape Results ResultSurface TEMPE
    GiD_Process $factors
    GiD_Process Mescape Results ContOptions ChangeTitleLeg $Text(lab036) escape
    GiD_Process Mescape 
}

proc AlyaPostFacto { } {
    global   textfactor
    global   Text Color
    catch {destroy .gid.faventfactor}

    set      W .gid.faventfactor
    AlyaCreateGeneralWindow $W $Text(tit010) 400 140 .gid

    frame    $W.all             -bg $Color(border)
    frame    $W.all.top         -bg $Color(border)

    frame    $W.all.top.dada1   -borderwidth 3 -relief groove
    label    $W.all.top.dada1.l -text $Text(lab273) -font {-weight bold -size 8}
    entry    $W.all.top.dada1.e -width 10 -textvariable textfactor
    pack     $W.all.top.dada1.l -side left  -padx 2m -pady 2m
    pack     $W.all.top.dada1.e -side right -padx 2m -pady 2m
    pack     $W.all.top.dada1   -side top -padx 2m -pady 2m 

    pack     $W.all.top -side top -padx 2m

    frame    $W.all.buttons   -bg $Color(border)
    button   $W.all.buttons.1 -text $Text(bot003) -width 10 -relief raised -command "AlyaCheckFactor"
    pack     $W.all.buttons.1 -side left     -padx 1m
    button   $W.all.buttons.2 -text $Text(bot002) -width 10 -relief raised -command "destroy $W"
    pack     $W.all.buttons.2 -side left     -padx 1m
    pack     $W.all.buttons   -side bottom   -padx 2m -pady 2m 

    pack     $W.all
    bind     $W.all.top.dada1.e <Return> "tkButtonInvoke $W.all.buttons.1"
    bind     $W.all.top.dada1.e <Escape> "tkButtonInvoke $W.all.buttons.2"
    focus    $W.all.top.dada1.e  
}

proc AlyaPostFactos { } {
    global   textfactora
    global   Text Color
    catch {destroy .gid.faventfactors}

    set      W .gid.faventfactors
    AlyaCreateGeneralWindow $W $Text(tit010) 400 140 .gid

    frame    $W.all             -bg $Color(border)
    frame    $W.all.top         -bg $Color(border)

    frame    $W.all.top.dada1   -borderwidth 3 -relief groove
    label    $W.all.top.dada1.l -text $Text(lab274) -font {-weight bold -size 8}
    entry    $W.all.top.dada1.e -width 10 -textvariable textfactors
    pack     $W.all.top.dada1.l -side left  -padx 2m -pady 2m
    pack     $W.all.top.dada1.e -side right -padx 2m -pady 2m
    pack     $W.all.top.dada1   -side top -padx 2m -pady 2m 

    pack     $W.all.top -side top -padx 2m

    frame    $W.all.buttons   -bg $Color(border)
    button   $W.all.buttons.1 -text $Text(bot003) -width 10 -relief raised -command "AlyaCheckFactors"
    pack     $W.all.buttons.1 -side left     -padx 1m
    button   $W.all.buttons.2 -text $Text(bot002) -width 10 -relief raised -command "destroy $W"
    pack     $W.all.buttons.2 -side left     -padx 1m
    pack     $W.all.buttons   -side bottom   -padx 2m -pady 2m 

    pack     $W.all
    bind     $W.all.top.dada1.e <Return> "tkButtonInvoke $W.all.buttons.1"
    bind     $W.all.top.dada1.e <Escape> "tkButtonInvoke $W.all.buttons.2"
    focus    $W.all.top.dada1.e  
}

proc AlyaCheckFactor { } {
    global textfactor Text
    set patron {^((-)?[0-9]+(.[0-9]+)?)$}
    if { [regexp $patron $textfactor] != 1 } {
        AlyaError $Text(lab111)
        return
    }
     AlyaPostVeloV
}

proc AlyaCheckFactors { } {
    global textfactors Text
    set patron {^((-)?[0-9]+(.[0-9]+)?)$}
    if { [regexp $patron $textfactors] != 1 } {
        AlyaError $Text(lab111)
        return
    }
    AlyaPostTEMPES
}

proc AlyaCreateGeneralWindow { w title x y wbg } {
    global        Color 
    global        Alya

    AlyaGetPaths

    catch         {destroy $w }
    toplevel      $w
    wm geometry   $w   +$x+$y
    wm resizable  $w   0 0 
    wm transient  $w   $wbg
    wm title      $w   $title
   # wm iconbitmap $w   [file join $Alya(ptypepath) pics Alya.ico]
    focus         $w
    #background of the window    
    set  Color(background) [$w cget -background]
    #border of the window
    set  Color(border)     [AlyaCCColorActivo $Color(background)]
    #boton importante
    set  Color(boton1)     $Color(background)
    set  Color(boton1fg)   $Color(background)
    #boton general           
    set  Color(boton2)     $Color(background)

    bind $w   <Control-w> "ShowWizard"
    bind $w   <Alt-w>     "AlyaHideWindows"
    bind .gid <Alt-w>     "AlyaHideWindows"
}

# Assign Boundary 2D/3D
proc AlyaBeforeMeshGeneration { } {
    set Ndime [AlyaGetNdime]
    if { $Ndime == 2 } {
        GiD_UnAssignData condition DOMAIN_element_2D surfaces 1:end
        foreach surface_id [GiD_Geometry list surface 1:end] {
            set values [list $surface_id]
            GiD_AssignData condition DOMAIN_element_2D surfaces $values [list $surface_id]
        }              
        GiD_UnAssignData condition DOMAIN_boundary_2D lines 1:end
        foreach surface_id [GiD_Geometry list surface 1:end] {
            set InfoSurface [GiD_Info list_entities surfaces $surface_id]
            set error [regexp -line {Meshing=No} $InfoSurface]
            if { $error != 1 } {
                foreach line_id [GidUtils::GetSurfaceLines $surface_id] {
                    incr LineCriteria($line_id)
                }                
            }
        }
        foreach line_id [GiD_Geometry list line 1:end] {
            if { $LineCriteria($line_id)==1 } { 
                set values [list $line_id]
                GiD_AssignData condition DOMAIN_boundary_2D lines $values [list $line_id]
            }
        }
    } else {
        GiD_UnAssignData condition DOMAIN_boundary_3D surfaces 1:end                       
        foreach volume_id [GiD_Geometry list volume 1:end] {
            set InfoVolume [GiD_Info list_entities volumes $volume_id]
            set error [regexp -line {Meshing=No} $InfoVolume]
            if { $error != 1 } {
                foreach surface_id [GidUtils::GetVolumeSurfaces $volume_id] {
                    incr SurfCriteria($surface_id)
                }
            }
        }
        foreach isurf [GiD_Geometry list surface 1:end] {
            if { $SurfCriteria($isurf)==1 } {
                set values [list $isurf]
                GiD_AssignData condition DOMAIN_boundary_3D surfaces $values [list $isurf]
            }
        }
    }
}

proc AlyaAskParametersTurmu { } {
    global TurmuK TurmuO TurmuE TurmuMut TurmuReT TurmuMutMu 
    global Text Color Alya

    set   w .gid.faventturmu
    AlyaCreateGeneralWindow $w $Text(tit025) 400 140 .gid

    frame         $w.tout                     -bg $Color(border)     
    frame         $w.tout.wall                -bg $Color(border)      

    label         $w.tout.wall.title          -text $Text(lab378)\n$Text(lab379) -bg $Color(border) -font {-weight bold -size 8} -justify left
    pack          $w.tout.wall.title          -side top -anchor w  -pady 1m 

    frame         $w.tout.wall.up             -borderwidth 3 -relief groove     
    # Flow scales
    frame         $w.tout.wall.up.left 
    label         $w.tout.wall.up.left.veloca -text [concat $Text(lab380) " K=  "] 
    grid   config $w.tout.wall.up.left.veloca -row 0 -column 0 -sticky w  
    entry         $w.tout.wall.up.left.velocb -relief sunken -width 10 -textvariable TurmuK
    grid   config $w.tout.wall.up.left.velocb -row 0 -column 1 -sticky w  
    label         $w.tout.wall.up.left.lengta -text [concat $Text(lab381) " Eps="]
    grid   config $w.tout.wall.up.left.lengta -row 1 -column 0 -sticky w  
    entry         $w.tout.wall.up.left.lengtb -relief sunken -width 10 -textvariable TurmuE
    grid   config $w.tout.wall.up.left.lengtb -row 1 -column 1 -sticky w   
    #label         $w.tout.wall.up.left.tempea -text [concat $Text(lab382) " Ome="] 
    #grid   config $w.tout.wall.up.left.tempea -row 2 -column 0 -sticky w  
    #entry         $w.tout.wall.up.left.tempeb -relief sunken -width 10 -textvariable TurmuO
    #grid   config $w.tout.wall.up.left.tempeb -row 2 -column 1 -sticky w  
    pack          $w.tout.wall.up.left        -side left -anchor w -padx 2m
    # Calculate button
    frame         $w.tout.wall.up.center      

    image  create photo imarrow  -file [file join $Alya(ptypepath)/pics/favent-arrow.gif] 
    #canvas $w.tout.wall.up.center.canvas -width 33 -height 33 -bd 0 
    #$w.tout.wall.up.center.canvas create image 2 2 -image imarrow -anchor nw
    #pack   $w.tout.wall.up.center.canvas -side top -anchor w -padx 2m -pady 2m
 
    button        $w.tout.wall.up.center.text  -image imarrow \
                  -height 21 -width 21 -relief raised -command "AlyaUpdateTurmu $w"
    pack          $w.tout.wall.up.center.text  -side top

    pack          $w.tout.wall.up.center       -side left -anchor w -padx 2m
    # Dimensionless parameters
    frame         $w.tout.wall.up.righ        
    label         $w.tout.wall.up.righ.mu1    -text  [concat $Text(lab383) " mut=rho*Cmu*fmu*k^2/eps="] -font {-weight bold -size 8}
    grid   config $w.tout.wall.up.righ.mu1    -row 0 -column 0 -sticky w  
    label         $w.tout.wall.up.righ.mu2    -text  $TurmuMut
    grid   config $w.tout.wall.up.righ.mu2    -row 0 -column 1 -sticky w

    label         $w.tout.wall.up.righ.mm1    -text  [concat $Text(lab385) " mut/mu="] -font {-weight bold -size 8}
    grid   config $w.tout.wall.up.righ.mm1    -row 1 -column 0 -sticky w
    label         $w.tout.wall.up.righ.mm2    -text  $TurmuMutMu
    grid   config $w.tout.wall.up.righ.mm2    -row 1 -column 1 -sticky w

    label         $w.tout.wall.up.righ.re1    -text  [concat $Text(lab384) " ReT=k^2/(nu*eps)="] -font {-weight bold -size 8}
    grid   config $w.tout.wall.up.righ.re1    -row 2 -column 0 -sticky w
    label         $w.tout.wall.up.righ.re2    -text  $TurmuReT
    grid   config $w.tout.wall.up.righ.re2    -row 2 -column 1 -sticky w

    pack          $w.tout.wall.up.righ        -side left -anchor w  -padx 2m

    pack          $w.tout.wall.up             -side top -anchor w -ipady 2m 
    pack          $w.tout.wall                -padx 2m -pady 2m 
    pack          $w.tout 

    # Close button
    frame         $w.buttons        -bd 5 -bg [AlyaCCColorActivo [$w cget -background]]
    button        $w.buttons.close  -text $Text(bot002) -height 1 -width 10 -relief raised -command "destroy $w" 
    pack          $w.buttons.close  -pady 1m
    pack          $w.buttons        -side top -fill both -expand yes


    bind         $w <Return> "tkButtonInvoke $w.tout.wall.up.center.text"
    bind         $w <Escape> "destroy $w"

}
proc AlyaCalculateTurmu { } {
    global TurmuK TurmuO TurmuE TurmuMut TurmuReT TurmuMutMu
    global Text

    if { [catch {set results [expr $TurmuK]}] != 0 } { 
        AlyaError [= "K must be a real number"] 
        return
    }
    if { [catch {set results [expr $TurmuE]}] != 0 } { 
        AlyaError [= "Eps must be a real number"] 
        return
    }
    if { [catch {set results [expr $TurmuO]}] != 0 } { 
        AlyaError [= "Ome must be a real number"] 
        return
    } 

    set Prbdata             [GiD_Info gendata]
    set PrbdataDensity      [lindex $Prbdata [expr [lsearch $Prbdata "TUR_Density="]+1 ] ]
    set PrbdataViscosity    [lindex $Prbdata [expr [lsearch $Prbdata "TUR_Viscosity="]+1 ] ]
#    set PrbTurbModel        [lindex $Prbdata [expr [lsearch $Prbdata "TUR_model:"]+1 ] ]

    set TurmuReT   [expr $PrbdataDensity*$TurmuK*$TurmuK/($PrbdataViscosity*$TurmuE)]

    #if { $PrbTurbModel == "Launder_Sharma_k_epsilon" } {
        set pop [expr (1.0+$TurmuReT/50.0)*(1.0+$TurmuReT/50.0)]
        set fmu [expr exp(-3.4/$pop)]
        #set fmu 1.0
#        WarnWinText $fmu
#    } else {
#        set fmu 1.0
#        WarnWinText "hola 2"
#    }
    set TurmuMut   [expr 0.09*$fmu*$PrbdataDensity*$TurmuK*$TurmuK/$TurmuE]
    set TurmuMutMu [expr $TurmuMut/$PrbdataViscosity]

    set TurmuMut   [format "%5.3E" $TurmuMut] 
    set TurmuMutMu [format "%5.3E" $TurmuMutMu] 
    set TurmuReT   [format "%5.3E" $TurmuReT] 
}
proc AlyaUpdateTurmu { w } {
    global TurmuK TurmuO TurmuE TurmuMut TurmuReT TurmuMutMu
    AlyaCalculateTurmu
    $w.tout.wall.up.righ.mu2 conf -text $TurmuMut
    $w.tout.wall.up.righ.mm2 conf -text $TurmuMutMu
    $w.tout.wall.up.righ.re2 conf -text $TurmuReT

}
