proc AlyaNASTINPostprocess { {type "INSIDELEFT"} } {
    global GidPriv
    global ToolbarsPriv
    global Alya
    global AlyaCadBitmapsNames4 AlyaCadBitmapsCommands4 AlyaCadBitmapsHelp4 GIDDEFAULT

    set dir [ file join $Alya(ptypepath) styles ]

    set AlyaCadBitmapsNames4(0) "          \
                alya-convergence.gif   \
                ---                    \
                alya-VELOC-vector.gif  \
                alya-VELOC-fill.gif    \
                alya-STREA-line.gif   \
                alya-PRESS-fill.gif    "
    set AlyaCadBitmapsCommands4(0) [list   \
                { -np- AlyaGlobalConvergence}  \
                "" \
                {-np- AlyaPostVectVELOC}   \
                {-np- AlyaPostFillVELOC}   \
                {-np- AlyaPostLineSTREA}   \
                {-np- AlyaPostFillPRESS}   ] 
    set AlyaCadBitmapsHelp4(0) [list  \
                [_ "Plot TEMPER convergence"] \
                [_ "Plot global convergence"] \
                [_ ""] \
                [_ "VELOC vectors"]  \
                [_ "VELOC module fill"] [_ "STREA  contour lines"] \
                [_ "PRESS fill"] ]
    set prefix Post
    set geomname ${prefix}AlyaPostBar
    set Alya(toolbarwin,$geomname) [CreateOtherBitmaps AlyaPostBar "AlyaCad post bar" \
          AlyaCadBitmapsNames4 AlyaCadBitmapsCommands4 \
          AlyaCadBitmapsHelp4 $dir AlyaCadBitmaps4 \
          $type $prefix]
    AddNewToolbar "AlyaCad post bar" $geomname AlyaCadBitmaps4
}

proc AlyaTEMPERPostprocess { {type "INSIDELEFT"} } {
    global GidPriv
    global ToolbarsPriv
    global Alya
    global AlyaCadBitmapsNames5 AlyaCadBitmapsCommands5 AlyaCadBitmapsHelp5 GIDDEFAULT

    set dir [ file join $Alya(ptypepath) styles ]

    set AlyaCadBitmapsNames5(0) "              \
                alya-convergence.gif           \
                alya-TEMPE-cvg.gif             \
                ---                            \
                alya-TEMPE-fill.gif            \
                alya-TEMPE-surf.gif            "
    set AlyaCadBitmapsCommands5(0) [list       \
                { -np- AlyaGlobalConvergence}  \
                { -np- AlyaConvergence 2}      \
                ""                             \
                {-np- AlyaPostFillTEMPE}       \
                {-np- AlyaPostSurfTEMPE}       ] 
    set AlyaCadBitmapsHelp5(0) [list           \
                [_ "Plot global convergence"]  \
                [_ "Plot TEMPE convergence"]   \
                [_ ""]                         \
                [_ "TEMPER contour fill"]      \
                [_ "TEMPER result surface"]    ]
    set prefix Post
    set geomname ${prefix}AlyaPostBar
    set Alya(toolbarwin,$geomname) [CreateOtherBitmaps AlyaPostBar "AlyaCad post bar" \
          AlyaCadBitmapsNames5 AlyaCadBitmapsCommands5 \
          AlyaCadBitmapsHelp5 $dir AlyaCadBitmaps5 \
          $type $prefix]
    AddNewToolbar "AlyaCad post bar" $geomname AlyaCadBitmaps5
}
