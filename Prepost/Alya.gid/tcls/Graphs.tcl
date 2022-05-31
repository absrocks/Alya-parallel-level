#################################################################
#
#  This procedure sets variables for ShowGraphic
#
################################################################
proc   SetGraphic { } {
    global xpoints ypoints ngraph npoingraph 
    global messagegraphic messagexaxis messageyaxis nlabelsx nlabelsy \
	    ticksx ticksy formatx formaty plotgraphic showgrid color \
	    typeaxisX typeaxisY xmax xmin ymax ymin ncolu messagex messagey
    global xmaxX xminX yminY ymaxY nzery whichplot

    #
    # 1) Finds maximum and minimum values.
    #
    # 1.1) xmin
    if { $xminX == "Automatic" } {
	set xmin 1.0e10
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    if  { $npoingraph($k) == 0 } { set xmin 0.0 }
	    for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		if { "$xpoints($i,$k)" < $xmin } { set xmin $xpoints($i,$k) }
	    }
	}
    } elseif { $xminX == "Automatic2" } {
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    set npoin 0
	    for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
	    if  { $npoingraph($k) == 0 } { set xmin 0.0 }
		if { "$xpoints($i,$k)" >= $xmin } { 
		    set xpoints($npoin,$k) $xpoints($i,$k)
		    set ypoints($npoin,$k) $ypoints($i,$k)
		    incr npoin 1
		}
	    }
	    set npoingraph($k) $npoin
	}
    }
    # 1.2) xmax 
    if { $xmaxX == "Automatic" } {
	set xmax -1.0e10
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    if  { $npoingraph($k) == 0 } { set xmax 1.0 }
	    for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		if { "$xpoints($i,$k)" > $xmax } { set xmax $xpoints($i,$k) }
	    }
	}
    } elseif { $xmaxX == "Automatic2" } {
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    set npoin 0
	    if  { $npoingraph($k) == 0 } { set xmax 1.0 }
	    for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		if { "$xpoints($i,$k)" <= $xmax } { 
		    set xpoints($npoin,$k) $xpoints($i,$k)
		    set ypoints($npoin,$k) $ypoints($i,$k)
		    incr npoin 1
		}
	    }
	    set npoingraph($k) $npoin
	} 
    } 
    for {set k 0} {$k < $ngraph} {incr k 1} { set nzery($k) 0 }
    #if { $whichplot == 0 } {
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		if { $ypoints($i,$k) < 1.0e-15 } { 
		    incr nzery($k)
		}
	    }
	    if { $nzery($k) == $npoingraph($k) } { 
		set nzery($k) 1 
	    } else {
		set nzery($k) 0
	    }
	}
    #}
    # 1.3) ymin
    if { $yminY == "Automatic" } {
	set ymin 1.0e10
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    if { $nzery($k) == 0 } {
		if  { $npoingraph($k) == 0 } { set ymin 0.0 }
		for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		    if { "$ypoints($i,$k)" < "$ymin" } { set ymin $ypoints($i,$k) }
		}
	    }
	}
    } elseif { $yminY == "Automatic2" } {
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    set npoin 0
	    if  { $npoingraph($k) == 0 } { set ymin 0.0 }
	    for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		if { "$ypoints($i,$k)" >= $ymin } { 
		    set xpoints($npoin,$k) $xpoints($i,$k)
		    set ypoints($npoin,$k) $ypoints($i,$k)
		    incr npoin 1
		}
	    }
	    set npoingraph($k) $npoin
	}
    }

    # 1.4) ymax
    if { $ymaxY == "Automatic" } {
	set ymax -1.0e10
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    if  { $npoingraph($k) == 0 } { set ymax 0.0 }
	    for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		if { "$ypoints($i,$k)" > "$ymax" } { set ymax $ypoints($i,$k) }
	    }
	}
    } elseif { $ymaxY == "Automatic2" } {
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    set npoin 0
	    if  { $npoingraph($k) == 0 } { set ymax 0.0 }
	    for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		if { "$ypoints($i,$k)" <= $ymax } { 
		    set xpoints($npoin,$k) $xpoints($i,$k)
		    set ypoints($npoin,$k) $ypoints($i,$k)
		    incr npoin 1
		}
	    }
	    set npoingraph($k) $npoin
	} 
    } 

}
#######################################################################
#
#  This process plots graphics in a canvas widget
#
#  Following variables must be set prior to the call by another 
#  procedure (see SetGraphic):
#
#######################################################################
proc   ShowGraphic { W } {
    global xpoints ypoints ngraph npoingraph 
    global messagegraphic messagexaxis messageyaxis nlabelsx nlabelsy \
	    ticksx ticksy formatx formaty plotgraphic showgrid color \
	    typeaxisX typeaxisY xmax xmin ymax ymin ncolu messagex messagey 
    global whichplot nzery
    global graphcolor forecolor
    global factorx lxmin lxmax xi xf
    global factory lymin lymax yi yf

    #
    # Canvas size
    #
    set areawidth  530
    set areaheight 470
    set xorigin    70
    set yorigin    [expr $areaheight-30]
    set xi         $xorigin
    set xf         [expr $areawidth-10]
    set yi         $yorigin
    set yf         70
    set xl         10
    set yl         20
    #
    # Cleans canvas window (for refresh)
    #
    destroy $W
    canvas  $W -bg $graphcolor -width $areawidth -height $areaheight -bg $graphcolor -borderwidth 0 -highlightthickness 0 
    pack    $W -side right -expand true -ipady 2m -ipadx 2m
    #
    # Inicializations.
    #
    if {[array exists xposi]==1} {[array unset xposi]} 
    if {[array exists yposi]==1} {[array unset yposi]}
    set showgridx "No"
    set showgridy "No"
    if { $showgrid == "X Axis" } { set showgridx "Yes" }
    if { $showgrid == "Y Axis" } { set showgridy "Yes" }
    if { $showgrid == "Both"   } {  
	    set showgridx "Yes" 
	    set showgridy "Yes"
    }
	
    #set messagexaxis "Number of iterations"
    #set messageyaxis "L2 residual norm"
    
    #
    # In the case of logarithmic scale, checks for incompatibility of data and
    # changed the values of xmin, xmax, ymin and ymax.
    #
    set lxmin $xmin
    set lxmax $xmax
    set lymin $ymin
    set lymax $ymax

    if { $typeaxisX == "Logarithmic" } { 
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		if { "$xpoints($i,$k)" < "0.0" } { 
		    tk_messageBox -icon error -type ok -message "Zero or negative residual"
		    return
		}
	    }
	}
	if { $xmin > 0.0 } { set lxmin [expr log10($xmin)] } else { set lxmin 0.0 }
	if { $xmax > 0.0 } { set lxmax [expr log10($xmax)] } else { set lxmax 1.0 }
	set formatx  "%3.2f" 
	unset  messagexaxis 
	append messagexaxis "Log10 ( " $messagex " ) "
    }
    if { $typeaxisY == "Logarithmic" } { 
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    if { $nzery($k) == 0 } {
		for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		    if { $ypoints($i,$k) < 1.0e-15 } {
			set ypoints($i,$k) 1.0e-12
			#tk_messageBox -icon error -type ok -message "Zero or negative value"
			#return
		    }
		}
	    }
	}
	if { $ymin > 0.0 } { 
	    set lymin [expr log10($ymin)] 
	} else { 
	    set lymin 0.0 
	}
	if { $ymax > 0.0 } { 
	    set lymax [expr log10($ymax)] 
	} else { 
	    set lymax 1.0 
	}
	#set formaty  "%3.2f"
	unset  messageyaxis 
	append messageyaxis "Log10 ( " $messagey " ) "
    }

    set pgraph 0
    for {set k 0} {$k < $ngraph} {incr k 1} { 
	if { $nzery($k) == 0 } { incr pgraph }
    }

    #
    # Puts message.
    #
    set m 0
    set i 1
    if { $whichplot != 9999 } {
	#$W create text 10 15 -text "$messagegraphic(0)" -anchor w  -fill $forecolor
	if { $pgraph < 6 } {
	    for {set k 0} {$k < $ngraph} {incr k 1} {
		if { $nzery($k) == 0 } {
		    incr m
		    $W create rectangle [expr $xl+100*$m] $yl [expr $xl+10+100*$m] [expr $yl+10] \
			-fill $color($k) -tag graphlegend
		    $W create text      [expr $xl+15+100*$m]  [expr $yl+5] -text "$messagegraphic($i)" \
			-anchor w -fill $forecolor -tag graphlegend
		} 
		incr i 1
	    }
	} elseif { $pgraph < 8 } {
	    for {set k 0} {$k < $ngraph} {incr k 1} {
		if { $nzery($k) == 0 } {
		    incr m
		    $W create rectangle [expr $xl+60*$m] $yl [expr $xl+10+60*$m] [expr $yl+10] \
			-fill $color($k) -tag graphlegend
		    $W create text      [expr $xl+15+60*$m]  [expr $yl+5] -text "$messagegraphic($i)" \
			-anchor w -fill $forecolor -tag graphlegend
		}
		incr i 1
	    }
	    
	} else {
	    set  i 1
	    for {set k 0} {$k < 5} {incr k 1} {
		$W create rectangle [expr $xl+60*$k] [expr $yl-15] [expr $xl+10+60*$k] [expr $yl+10-15] \
                                    -fill $color($k) -tag graphlegend
		$W create text      [expr $xl+15+60*$k]  [expr $yl+5-15] -text "$messagegraphic($i)" \
                                    -anchor w -fill $forecolor -tag graphlegend
		incr i 1
	    }
	    for {set k 0} {$k < [expr $ngraph-5]} {incr k 1} {
		set k1 [expr $k+5]
		$W create rectangle [expr $xl+60*$k] [expr $yl] [expr $xl+10+60*$k] [expr $yl+10] \
                                    -fill $color($k1) -tag graphlegend
		$W create text      [expr $xl+15+60*$k]  [expr $yl+5] -text "$messagegraphic($i)" \
                                    -anchor w -fill $forecolor -tag graphlegend
		incr i 1
	    }
	    	}
	
    }
    #
    # Checks for constant plots. 
    #
    if { $lymin >= $lymax } {
	set lymin [expr $lymin-1]
	set lymax [expr $lymax+1]
	set nlabelsy 6
	set ticksy   0.5
    }
    if { $lxmin >= $lxmax } {
	set lxmin [expr $lxmin-1]
	set lxmax [expr $lxmax+1]
	set nlabelsx 6
	set ticksx   0.5
    }
    #
    # Scales values to fit them inside print area 
    #
    set factorx [expr (double($xf)-double($xi))/(double($lxmax)-double($lxmin))]
    set factory [expr (double($yf)-double($yi))/(double($lymax)-double($lymin))]

    if { $typeaxisX == "Linear" } { 
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		set xposi($i,$k) [ expr $factorx*($xpoints($i,$k)-$lxmin)+$xi ]
	    }
	}
    } else {
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		set xposi($i,$k) [ expr $factorx*(log10($xpoints($i,$k))-$lxmin)+$xi ]
	    }
	}
    }
    if { $typeaxisY == "Linear" } { 
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		set yposi($i,$k) [ expr $factory*($ypoints($i,$k)-$lymin)+$yi ]
	    }
	}
    } else {
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		set yposi($i,$k) [ expr $factory*(log10($ypoints($i,$k))-$lymin)+$yi ]
	    }
	}
    }
    #
    # Plots axis and messages.
    #
    # x axis.
    #
    if { "$lymin" >= 0.0 } {
	set xaxisposi $yi
    } elseif { "$lymax" <= 0.0 } {
	set xaxisposi $yf
    } else {
	set xaxisposi [ expr $yi-$factory*$lymin ]
    }
    $W create text [expr ($xi+$xf)/2] [expr $yi+25] -text $messagexaxis \
                   -fill $forecolor -tag messagexaxis

    #
    # y axis.
    #
    if { "$xmin" >= 0.0 } {
	set yaxisposi $xi
    } elseif { "$xmax" <= 0.0 } {
	set yaxisposi $xf
    } else {
	set yaxisposi [ expr $xi-$factorx*$lxmin ]
    }
    $W create line $xi $yi $xf $yi  -width 1 -fill $forecolor
    $W create line $xf $yi $xf $yf  -width 1 -fill $forecolor
    $W create line $xf $yf $xi $yf  -width 1 -fill $forecolor
    $W create line $xi $yf $xi $yi  -width 1 -fill $forecolor
    $W create text $yaxisposi [expr $yf-20] -text $messageyaxis -anchor w \
              -fill $forecolor -tag messageyaxis

    #
    # Prints line or symbols.
    #
    #$W config -scrollregion "$xi $yf $xf $yi"

    for {set k 0} {$k < $ngraph} {incr k 1} {
	if { $nzery($k)==0 } {
	    if {"$plotgraphic($k)" == "Line" } then {
		if {$npoingraph($k) > 0} {
		    set x2 $xposi(0,$k)
		    set y2 $yposi(0,$k)
		    for {set i 1} {$i < $npoingraph($k)} {incr i 1} {
			set x1 $x2
			set y1 $y2
			set x2 $xposi($i,$k)
			set y2 $yposi($i,$k)
			if { $x1 <= $xf && $x1 >= $xi && $x2 <= $xf && $x2 >= $xi && \
				 $y1 <= $yi && $y1 >= $yf && $y2 <= $yi && $y2 >= $yf } {
			    $W create line "$x1 $y1 $x2 $y2" -fill $color($k)
			}
		    }
		} 
	    } elseif {"$plotgraphic($k)" == "Points"} {
		for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		    set x1 $xposi($i,$k)
		    set y1 $yposi($i,$k)
		    if { $x1 <= $xf && $x1 >= $xi && $y1 <= $yi && $y1 >= $yf  } {
			$W create oval [expr $x1-2] [expr $y1-2] \
			    [expr $x1+2] [expr $y1+2] \
			    -fill $color($k) 
		    }		
		}
	    } else { }
	}
    }
    #
    # Number of labels
    #
    set npoin 0 
    for {set k 0} {$k < $ngraph} {incr k 1} {
	if { "$npoingraph($k)" > "$npoin" } { set npoin $npoingraph($k) } 
    } 
    if { $npoin <= 1 } { set npoin 2 }
    set nlabelsx $npoin
    if { $nlabelsx > 6 } { set nlabelsx 6 } 
    set ticksx [expr double($lxmax-$lxmin)/($nlabelsx-1) ] 
    
    set nlabelsy $npoin
    if { $nlabelsy > 6 } { set nlabelsy 6 } 
    set ticksy [expr double($lymax-$lymin)/($nlabelsy-1) ]
    #
    # Puts labels and dashed lines.
    # 
    for {set i 0} {$i < $nlabelsx} {incr i 1} {
	set xlabel [expr $lxmin+($i*$ticksx)]  
	set vlabel [format "$formatx"  $xlabel]
	set xlabel [expr $factorx*($xlabel-$lxmin)+$xi ]
	#$W create line "$xlabel $xaxisposi $xlabel [expr $xaxisposi-5]" -width 1
	$W create line "$xlabel $yi $xlabel [expr $yi-5]" -width 1
	if {$showgridx == "Yes" } {
	    $W create line " $xlabel $yi $xlabel $yf " -dash .
	}
	$W create text "$xlabel [expr $yi+12]" -text "$vlabel" -fill $forecolor
    }
    for {set i 0} {$i < $nlabelsy} {incr i 1} {
	set ylabel [expr $lymin+($i*$ticksy)]
	set vlabel [format "$formaty"  [expr ($ylabel)]]
	set ylabel [expr $factory*($ylabel-$lymin)+$yi ]
	$W create line "$yaxisposi $ylabel [expr $yaxisposi+5] $ylabel " -width 1
	if {$showgridy == "Yes" } {
	    $W create line " $xi $ylabel $xf $ylabel " -dash .
	}
	$W create text "[expr $xi-35] $ylabel" -text "$vlabel" -fill $forecolor
    }
    #pack   $W -side top  -expand true 
    #pack   .gid.postwin.plot         -side top    -padx 2  -pady 2 -fill both -ipadx 2 -ipady 2 

    #$W bind messagexaxis <Button-1>  {AlyaGraphsMark %x %y %W}
    #$W bind messagexaxis <B1-Motion> {AlyaGraphsDrag %x %y %W}
    #$W bind messageyaxis <Button-1>  {AlyaGraphsMark %x %y %W}
    #$W bind messageyaxis <B1-Motion> {AlyaGraphsDrag %x %y %W}
    #$W bind graphlegend  <Button-1>  {AlyaGraphsMark %x %y %W}
    #$W bind graphlegend  <B1-Motion> {AlyaGraphsDrag %x %y %W}

}

proc AlyaBind { W } {
    $W config -cursor crosshair
    bind $W  <ButtonPress>   "AlyaZoomPress   %x %y" 
    bind $W  <ButtonRelease> "AlyaZoomRelease %x %y" 
    bind $W  <B1-Motion>     "AlyaZoomMotion  %x %y" 
    bind $W  <Escape>        "AlyaZoomEscape $W"
    bind $W  <B1-Escape>     "AlyaZoomEscape $W"
}

proc AlyaZoomEscape { W } {
    bind $W  <ButtonPress>   "" 
    bind $W  <ButtonRelease> "" 
    bind $W  <B1-Motion>     "" 
    $W config -cursor {}
}

proc AlyaZoomMotion { x y } {
    global xzoom1 yzoom1
    global xi xf 
    global yi yf
    global forecolor

    if { $x >= $xf } { set x $xf }
    if { $x <= $xi } { set x $xi }
    if { $y >= $yi } { set y $yi }
    if { $y <= $yf } { set y $yf }
    catch {.gid.postwin.plot.graphic delete pBox}
    .gid.postwin.plot.graphic create rectangle $xzoom1 $yzoom1 $x $y -tag {pBox pBoxx} -stipple gray25 -fill gray25
}

proc AlyaZoomPress { x y } {
    global xzoom1 yzoom1
    global xi xf forecolor
    global yi yf

    .gid.postwin.plot.graphic config -cursor crosshair
    set xzoom1 $x
    set yzoom1 $y
    if { $x > $xf } { set xzoom1 $xf }
    if { $x < $xi } { set xzoom1 $xi }
    if { $y < $yf } { set yzoom1 $yf }
    if { $y > $yi } { set yzoom1 $yi }

}

proc AlyaZoomRelease { x y } {
    global xzoom1 yzoom1
    global xzoom2 yzoom2
    global factorx lxmin lxmax xi xf
    global factory lymin lymax yi yf
    global xminX xmaxX yminY ymaxY 
    global xmin  xmax  ymin  ymax 
    global AlyaAutomaticUpdate
    global typeaxisX typeaxisY

    .gid.postwin.plot.graphic config -cursor {}

    set xzoom2 $x
    set yzoom2 $y

    if { $xzoom2 > $xf } { set xzoom2 $xf }
    if { $xzoom2 < $xi } { set xzoom2 $xi }
    if { $yzoom2 < $yf } { set yzoom2 $yf }
    if { $yzoom2 > $yi } { set yzoom2 $yi }

    set xminX "Fixed"
    set xmaxX "Fixed"
    set yminY "Fixed"
    set ymaxY "Fixed"

    if { $typeaxisX == "Linear" } { 
	set xzoom1 [expr ($xzoom1-$xi)/$factorx+$lxmin]
	set xzoom2 [expr ($xzoom2-$xi)/$factorx+$lxmin]
    } else {
	set xzoom1 [expr pow(10,($xzoom1-$xi)/$factorx+$lxmin)]
	set xzoom2 [expr pow(10,($xzoom2-$xi)/$factorx+$lxmin)]	
    }
    if { $typeaxisY == "Linear" } { 
	set yzoom1 [expr ($yzoom1-$yi)/$factory+$lymin]
	set yzoom2 [expr ($yzoom2-$yi)/$factory+$lymin]
    } else {
	set yzoom1 [expr pow(10,($yzoom1-$yi)/$factory+$lymin)]
	set yzoom2 [expr pow(10,($yzoom2-$yi)/$factory+$lymin)]
    }

    if { $xzoom1 > $xzoom2 } {
	set xzoomt $xzoom1
	set xzoom1 $xzoom2
	set xzoom2 $xzoomt
    }
    if { $yzoom1 > $yzoom2 } {
	set yzoomt $yzoom1
	set yzoom1 $yzoom2
	set yzoom2 $yzoomt
    }

    set xmin  $xzoom1
    set xmax  $xzoom2
    set ymin  $yzoom1
    set ymax  $yzoom2

    #SetGraphic
    ShowGraphic .gid.postwin.plot.graphic
    if { [winfo exists .gid.postwin.plot.graphic]== 1 && $AlyaAutomaticUpdate == 1} { 
	after 8000 AlyaRefreshWindow .gid.postwin.plot.graphic
    }
    #AlyaRefreshWindow .gid.postwin.plot.graphic 
    #ShowGraphic .gid.postwin.plot.graphic 
}

proc AlyaSelectTorque { w } {
    global ListTorque
    global Ndime
    if { $Ndime == 2 } {
	set entities [GidUtils::PickEntities lines multiple]
    } else {
	set entities [GidUtils::PickEntities surfaces multiple]
    }
    set ListTorque ""
    foreach { aux } $entities {
	set aux [split $aux /:]
	set longaux [llength $aux]
	if { $longaux > 1 } {
	    for { set j [lindex $aux 0] } { $j <= [lindex $aux 1] } { incr j 1 } {
		set ListTorque [concat $ListTorque $j]
	    }
	} else {
	    set ListTorque [concat $ListTorque $aux]
	}
    }
    $w config -textvariable ListTorque
}

proc AlyaReadTorque { } {  
    global ListTorque FactorTorque
    global Alya
    global number
    global messagegraphic messagexaxis messageyaxis nlabelsx nlabelsy \
	   ticksx ticksy formatx formaty plotgraphic showgrid color \
	   typeaxisX typeaxisY xmax xmin ymax ymin ncolu messagex messagey
    global xpoints ypoints npoingraph ngraph
    global whichplot 
    global Text

    set Torque $ListTorque
    set TorqueNumber [llength $Torque]
    for {set k 0} {$k < $TorqueNumber} {incr k 1} {
	set Torque [lreplace $Torque $k $k [string trim [lindex $Torque $k] " "]]
	set TorqueLength [string length [lindex $Torque $k]]
	if { $TorqueLength == 1 } {
	    set CurrentTorque "000000"
	} elseif { $TorqueLength == 2 } {
	    set CurrentTorque "00000"
	} elseif { $TorqueLength == 3 } {
	    set CurrentTorque "0000"
	} elseif { $TorqueLength == 4 } {
	    set CurrentTorque "000"
	} elseif { $TorqueLength == 5 } {
	    set CurrentTorque "00"
	} elseif { $TorqueLength == 6 } {
	    set CurrentTorque "0"
	} elseif { $TorqueLength == 7 } {
	    set CurrentTorque ""
	}
	append CurrentTorque [lindex $Torque $k]
	set Torque [lreplace $Torque $k $k $CurrentTorque]
    }

    # Open convergence file
    set extension ".frc"
    AlyaGetPaths
    for {set k 0} {$k < $TorqueNumber} {incr k 1} {
	set filepath [file join $Alya(resultspath) $Alya(projectname)-[lindex $Torque $k]$extension]
	if {[file exists $filepath]!=1} {
	    for {set j 0} {$j < $ngraph} {incr j 1} { set npoingraph($j) 0.0 }
	    return
	} else {
	    set luout [open $filepath r 0600]
	} 
	set nline 0
	if { $k == 0 } {
	    while { ![eof $luout] } {
		set linea  [gets $luout]
		set ndime  [NumberDecode $linea]
		if { $ndime > "1" } { 
		    set nlin1  [expr $nline+1]
		    for {set j 0} {$j < $ngraph} {incr j 1} {
			set xpoints($nline,$j) $number(2)
			set ypoints($nline,$j) [expr $FactorTorque*$number($ncolu($j))]
		    }
		    incr nline 1
		}   
	    }
	} else {
	    while { ![eof $luout] } {
		set linea  [gets $luout]
		set ndime  [NumberDecode $linea]
		if { $ndime > "1" } { 
		    set nlin1  [expr $nline+1]
		    for {set j 0} {$j < $ngraph} {incr j 1} {
			set xpoints($nline,$j) $number(2)
			set ypoints($nline,$j) [expr $ypoints($nline,$j)+$FactorTorque*$number($ncolu($j))]
		    }
		    incr nline 1
		}   
	    }
	}
	for {set j 0} {$j < $ngraph} {incr j 1} {    
	    set npoingraph($j) $nline
	}
	close $luout
    }
}

#########################################################################
#
#  This procedure fills Alya canvas postprocess window
#
#########################################################################
proc AlyaRefreshWindow { w } { 
    global whichplot
    global AlyaAutomaticUpdate

    if { [winfo exists $w]!= 1} { return }

    # Gets data from a file
    set Iresu [AlyaReadConvergence]
    if { $Iresu==-1 } { return -1 }

    # Sets and show graphic data
    set Iresu [SetGraphic]
    if { $Iresu==-1 } { return -1 }
   
    # Show and refresh window
    ShowGraphic $w

    if { [winfo exists $w]== 1 && $AlyaAutomaticUpdate == 1} { 
	after 8000 AlyaRefreshWindow $w 
    }
}

#########################################################################
#
#  This procedure inicializes graphic settings. Some of the values can
#  be modified during execution using Alya graphic settings window.
#  Values than can be modified during execution are:
#
#  formatx formaty
#  showgrid
#  plotgraphic
#  linearaxisX linearaxisY
#  color
#
#  The rest of variables are parameters fixed for each particular problem.
#
#########################################################################
proc AlyaInitialSettings { } {
    global projectpath projectname faventtclpath practicaname
    global messagegraphic messagexaxis messageyaxis nlabelsx nlabelsy \
           ticksx ticksy formatx formaty plotgraphic showgrid color \
           typeaxisX typeaxisY xmax xmin ymax ymin ncolu messagex messagey
    global xpoints ypoints ngraph npoingraph graphtitle
    global eixX eixY choiceplot choicecolor graphcolor forecolor
    global formatxX formatyY showgridX showgridY 
    global xmaxX xminX yminY ymaxY
    global message1 message2
    global messagetaula
    global whichplot 
    global Text
    global GraphNumberZones
    global GraphNameZones


    # Convergence
    if { $whichplot == 0 } {
	set ngraph 0
	set messagegraphic(0)        ""
	set graphtitle($whichplot)   "Global convergence" 
	#
	#   NASTIN
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          4
	set  messagegraphic($ngraph) "NASTIN"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "red"
	#
	#   TEMPER
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          5
	set  messagegraphic($ngraph) "TEMPER"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "blue"
	#set  color($ngrap1)          "#FF00FF"
	#
	#   CODIRE
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          6
	set  messagegraphic($ngraph) "CODIRE"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "green"
	#set  color($ngrap1)          "#FF00FF"
	#
	#   TURBUL
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          7
	set  messagegraphic($ngraph) "TURBUL"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "black"
	#   
	#   EXMEDI
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          8
	set  messagegraphic($ngraph) "EXMEDI"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "orange"
	#   
	#   NASTAL
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          9
	set  messagegraphic($ngraph) "NASTAL"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "yellow"
	#   
	#   GOTITA
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          14
	set  messagegraphic($ngraph) "GOTITA"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "grey"
	# 
	#   WAVEQU
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          15
	set  messagegraphic($ngraph) "WAVEQU"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "darkgreen"
	# 
	#   LEVELS
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          16
	set  messagegraphic($ngraph) "LEVELS"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "orange"
  
    } elseif { $whichplot == 1 } {
	set ngraph 0
	set messagegraphic(0)        ""
	set graphtitle($whichplot)   "NASTIN convergence" 
	#
	#   Velocity
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          5
	set  messagegraphic($ngraph) "Velocity"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "red"
	#
	#   Pressure
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          6
	set  messagegraphic($ngraph) "Pressure"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "blue"

    } elseif { $whichplot == 2 } {
	set ngraph 0
	set messagegraphic(0)        ""
	set graphtitle($whichplot)   "TEMPER convergence" 
	#
	#   Temperature
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          5
	set  messagegraphic($ngraph) "Temperature"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "red"

    } elseif { $whichplot == 3 } {
	set ngraph 0
	set messagegraphic(0)        ""
	set graphtitle($whichplot)   "CODIRE convergence" 
	#
	#   CDR unknown
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          5
	set  messagegraphic($ngraph) "CDR Unknown"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "red"

    } elseif { $whichplot == 4 } {
	set ngraph 0
	set messagegraphic(0)        ""
	set graphtitle($whichplot)   "TURBUL convergence" 
	#
	#   Turbulence 1
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          5
	set  messagegraphic($ngraph) "Turbulence variable 1"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "red"
	#
	#   Turbulence 2
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          6
	set  messagegraphic($ngraph) "Turbulence variable 2"

    } elseif { $whichplot == 5 } {
    } elseif { $whichplot == 6 } {

    } elseif { $whichplot == 11 } {
	set ngraph 0
	set messagegraphic(0)        ""
	set graphtitle($whichplot)   "GOTITA convergence" 
	#
	#   Water volume fraction
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          5
	set  messagegraphic($ngraph) "Droplet velocity"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "red"
	#
	#   Droplet velocity
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          6
	set  messagegraphic($ngraph) "Water volume fraction"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "blue"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "blue"

    } elseif { $whichplot == 12 } {
	set ngraph 0
	set messagegraphic(0)        ""
	set graphtitle($whichplot)   "WAVEQU convergence" 
	#
	#   Wave amplitude
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          5
	set  messagegraphic($ngraph) "Wave amplitude"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "red"

    } elseif { $whichplot == 13 } {
	set ngraph 0
	set messagegraphic(0)        ""
	set graphtitle($whichplot)   "LEVELS convergence" 
	#
	#   Wave amplitude
	#
	incr ngraph                  1
	set  ngrap1                  [expr $ngraph-1]
	set  ncolu($ngrap1)          5
	set  messagegraphic($ngraph) "Level set"
	set  plotgraphic($ngrap1)    "Line"
	set  color($ngrap1)          "red"
    }
    #
    #   Axes for convergence graphs
    #
    set messagexaxis   "Number of iterations"
    set messageyaxis   "Residual norm"
    set messagex       " "
    set messagey       " "
    set formatx        "%5.0f"
    set formaty        "%3.2f"
    set typeaxisX      "Linear"
    set typeaxisY      "Logarithmic" 
    set xmaxX          "Automatic"
    set xminX          "Automatic"
    set ymaxY          "Automatic"
    set yminY          "Automatic"

    set nlabelsx       6
    set nlabelsy       6
    set graphcolor     "white"
    set forecolor      "black"
    set showgrid       "None"   
    set showgridX      "No" 
    set showgridY      "No" 
    set eixX           $typeaxisX     
    set eixY           $typeaxisY
    set formatxX       $formatx
    set formatyY       $formaty
    for {set igraph 0} {$igraph < $ngraph} {incr igraph} {
	set choicecolor($igraph) $color($igraph)
	set choiceplot($igraph)  $plotgraphic($igraph)
    }
}

proc AlyaReadConvergence { } {  
    global Alya
    global number
    global messagegraphic messagexaxis messageyaxis nlabelsx nlabelsy \
	   ticksx ticksy formatx formaty plotgraphic showgrid color \
	   typeaxisX typeaxisY xmax xmin ymax ymin ncolu messagex messagey
    global xpoints ypoints npoingraph ngraph
    global whichplot 
    global Text

    if { $whichplot == 0 } {
    	set extension ".cvg"
    } elseif { $whichplot == 1 } {
     	set extension ".nsi.cvg"
    } elseif { $whichplot == 2 } {
    	set extension ".tem.cvg"
    } elseif { $whichplot == 3 } {
    	set extension ".cdr.cvg"
    } elseif { $whichplot == 4 } {
    	set extension ".tur.cvg"
    } elseif { $whichplot == 11 } {
    	set extension ".got.cvg"
    } elseif { $whichplot == 12 } {
    	set extension ".wav.cvg"
    } elseif { $whichplot == 13 } {
    	set extension ".lev.cvg"
    }

    # Open convergence file  
    AlyaGetPaths
    set filepath [file join $Alya(resultspath) $Alya(projectname)$extension]
    if {[file exists $filepath]!=1} {
	for {set j 0} {$j < $ngraph} {incr j 1} { set npoingraph($j) 0.0 }
	return
    } else {
	set luout [open $filepath r 0600]
    } 
    # Read line by line
    set nline 0
    while { ![eof $luout] } {
       set linea  [gets $luout]
       set ndime  [NumberDecode $linea]
       if { $ndime > "1" } { 
	 set nlin1  [expr $nline+1]
	 for {set j 0} {$j < $ngraph} {incr j 1} {
	     set xpoints($nline,$j) $nlin1
	     set ypoints($nline,$j) $number($ncolu($j))
	     if { $ypoints($nline,$j) <= 0.0  }   { set ypoints($nline,$j) 1.0e-16 }
	     if { $ypoints($nline,$j) <= 1.e-16 } { set ypoints($nline,$j) 1.0e-16 }
	 }
	   incr nline 1
       }   
    }
    for {set j 0} {$j < $ngraph} {incr j 1} {    
       set npoingraph($j) $nline
    }
    # Close file
    close $luout
}

#########################################################################
#
#  This procedure opens favent postproces window.
#
########################################################################
proc AlyaPostprocesWindow { } { 
    global whichplot
    global Text Color 
    global graphcolor graphtitle
    global Alya

    set W .gid.postwin
    AlyaCreateGeneralWindow $W $Text(tit002) 100 100 .gid
    #
    # Graphic settings (each graphic within each experiment can have its own settings).
    #
    AlyaInitialSettings 
    AlyaSetGraphMenus $W
    # 
    # Creates canvas for graphic. Canvas window is always .postwin.plot.graphic
    # 

    frame  $W.plot -bg $graphcolor -bd 2 -relief sunken

    AlyaGetPaths
    image create photo grzoomin -file [file join $Alya(ptypepath) styles autocad-zoom-frame.gif] 
    image create photo grzoomfr -file [file join $Alya(ptypepath) styles autocad-zoom-extends.gif] 
    image create photo grclosew -file [file join $Alya(ptypepath) styles autocad-tool-quit.gif] 
    image create photo grredraw -file [file join $Alya(ptypepath) styles autocad-redraw.gif] 
    frame  $W.plot.buttons -background $Color(border) 
    button $W.plot.buttons.grzoomin -image grzoomin -relief groove -bd 0 -background $Color(border) \
           -command "AlyaBind      $W.plot.graphic"
    pack   $W.plot.buttons.grzoomin -side left
    button $W.plot.buttons.grzoomfr -image grzoomfr -relief groove -bd 0 -background $Color(border) \
           -command "AlyaAutomatic $W.plot.graphic"
    pack   $W.plot.buttons.grzoomfr -side left
    button $W.plot.buttons.grredraw -image grredraw -relief groove -bd 0 -background $Color(border) \
           -command "AlyaRefresh   $W.plot.graphic" 
    pack   $W.plot.buttons.grredraw -side left
    button $W.plot.buttons.grclosew -image grclosew -relief groove -bd 0 -background $Color(border) \
           -command "destroy $W" 
    pack   $W.plot.buttons.grclosew -side left 
    pack   $W.plot.buttons -side top -anchor w -expand true -fill both 

    label $W.plot.header  -text $graphtitle($whichplot) -justify center -bg $graphcolor -font {-weight bold -size 8}

    pack   $W.plot.header  -side top  -ipadx 2 -ipady 2 -padx 2 -pady 2 
    canvas $W.plot.graphic -width 500 -height 470 -bg $graphcolor -borderwidth 20 -highlightthickness 20 -relief flat
    pack   $W.plot.graphic -side top  -expand true 
    #button $W.plot.close   -text $Text(bot002)  -command "destroy $W" -underline 0 -width 10  
    #pack   $W.plot.close   -side bottom -padx 1m -pady 1m
    pack   $W.plot         -side top    -padx 2  -pady 2 -fill both -ipadx 2 -ipady 2 
    #
    # Fills graphic.
    #
    set Iresu [AlyaRefreshWindow $W.plot.graphic]
    if { $Iresu==-1 } { destroy $W}
}

#########################################################################
proc AlyaReadYplus { } {  
    global Alya
    global number
    global messagegraphic messagexaxis messageyaxis nlabelsx nlabelsy \
           ticksx ticksy formatx formaty plotgraphic showgrid color \
           typeaxisX typeaxisY xmax xmin ymax ymin ncolu messagex messagey  
    global xpoints ypoints npoingraph ngraph
    global whichplot 
    global Text
    # Opens convergence file
    AlyaGetPaths
    set extension ".cvg"
    set filepath [file join $Alya(resultspath) $Alya(projectname)$extension]
    if {[file exists $filepath]!=1} {
	for {set j 0} {$j < $ngraph} {incr j 1} { set npoingraph($j) 0.0 }
	return
    } else {
	set luout [open $filepath r 0600]
    } 
    # Read line by line
    set nline 0
    while { ![eof $luout] } {
       set linea  [gets $luout]
       set ndime  [NumberDecode $linea]
       if { $ndime > "1" } { 
	 set nlin1  [expr $nline+1]
         for {set j 14} {$j < 34} {incr j 1} {
             set l             [expr $j-14]
             set xpoints($l,2) $nlin1
             set ypoints($l,2) $number($j)
         }
  	 incr nline 1
       }   
    }
    set k             0
    set npoingraph(0) 80
    set xpoints(0,0)  0.0
    set ypoints(0,0)  0.0
    for {set j 14} {$j < 34} {incr j 1} {
        set  l             [expr $j-14]
	# point 2
        incr k             1
        set  k1            [expr $k-1]
        set  x             $xpoints($k1,0)
        set  xpoints($k,0) $x
        set  ypoints($k,0) $ypoints($l,2)
	# point 3
        incr k             1
        set  x             [expr $x+20]
        set  xpoints($k,0) $x
        set  ypoints($k,0) $ypoints($l,2)
	# point 4
        incr k             1
        set  xpoints($k,0) $x
        set  ypoints($k,0) 0.0
	# point 1, next rectangle
        incr k             1
        set  xpoints($k,0) $x
        set  ypoints($k,0) 0.0
    }
    set  npoingraph(1) 8
    set  xpoints(0,1) 40 
    set  ypoints(0,1) 0.0
    set  xpoints(1,1) 40 
    set  ypoints(1,1) $ypoints(2,2) 
    set  xpoints(2,1) 60 
    set  ypoints(2,1) $ypoints(2,2) 
    set  xpoints(3,1) 60 
    set  ypoints(3,1) 0.0
    set  xpoints(4,1) 60 
    set  ypoints(4,1) $ypoints(3,2) 
    set  xpoints(5,1) 60 
    set  ypoints(5,1) $ypoints(3,2) 
    set  xpoints(6,1) 80 
    set  ypoints(6,1) $ypoints(3,2) 
    set  xpoints(7,1) 80 
    set  ypoints(7,1) 0.0

    # Close file
    close $luout
}

#########################################################################
proc AlyaReadVeloc { } {  
    global Alya
    global number
    global messagegraphic messagexaxis messageyaxis nlabelsx nlabelsy \
           ticksx ticksy formatx formaty plotgraphic showgrid color \
           typeaxisX typeaxisY xmax xmin ymax ymin ncolu messagex messagey
    global xpoints ypoints npoingraph ngraph
    global whichplot 
    global Text
    #  Opens witness file
    set extension ".wit"
    AlyaGetPaths
    set filepath [file join $Alya(resultspath) $Alya(projectname)$extension]
    if {[file exists $filepath]!=1} {
	for {set j 0} {$j < $ngraph} {incr j 1} { set npoingraph($j) 0.0 }
	return
	#AlyaError $Text(lab107)
	#return -1
    } else {
     set luout [open $filepath r 0600]
    } 
    # Read line by line
    set nline 0
    while { ![eof $luout] } {
       set linea  [gets $luout]
       set ndime  [NumberDecode $linea]
       if { $ndime > "1" } { 
         for {set j 0} {$j < $ngraph} {incr j 1} {
             set xpoints($nline,$j) $number(1)
             set ypoints($nline,$j) $number($ncolu($j))
         }
  	 incr nline 1
       }   
    }
    for {set j 0} {$j < $ngraph} {incr j 1} {    
       set npoingraph($j) $nline
    }
    # Close file
    close $luout
}
#########################################################################
proc AlyaReadTempe { } {  
    global Alya
    global number
    global messagegraphic messagexaxis messageyaxis nlabelsx nlabelsy \
           ticksx ticksy formatx formaty plotgraphic showgrid color \
           typeaxisX typeaxisY xmax xmin ymax ymin ncolu messagex messagey
    global xpoints ypoints npoingraph ngraph
    global whichplot 
    global Text
    # Open witness file 
    set extension ".wit"
    AlyaGetPaths
    set filepath [file join $Alya(resultspath) $Alya(projectname)$extension]
    if {[file exists $filepath]!=1} {
	for {set j 0} {$j < $ngraph} {incr j 1} { set npoingraph($j) 0.0 }
	return
	#AlyaError $Text(lab107)
	#return -1
    } else {
     set luout [open $filepath r 0600]
    } 
    # Read line by line
    set nline 0
    while { ![eof $luout] } {
       set linea  [gets $luout]
       set ndime  [NumberDecode $linea]
       if { $ndime > "1" } { 
	 set nlin1  [expr $nline+1]
         for {set j 0} {$j < $ngraph} {incr j 1} {
             set xpoints($nline,$j) $number(1)
             set ypoints($nline,$j) $number($ncolu($j))
         }
  	 incr nline 1
       }   
    }
    for {set j 0} {$j < $ngraph} {incr j 1} {    
       set npoingraph($j) $nline
    }
    # Close file
    close $luout
}

proc AlyaSetGraphMenus { W } {
    global messagegraphic messagex messagey nlabelsx nlabelsy \
	ticksx ticksy formatx formaty plotgraphic showgrid color \
	typeaxisX typeaxisY xmax xmin ymax ymin 
    global xpoints ypoints ngraph npoingraph 
    global eixX eixY choiceplot choicecolor 
    global xmaxX xminX yminY ymaxY 
    global formatxX formatyY showgridX showgridY 
    global LanguageMaster   
    global graphcolor forecolor 
    global AlyaAutomaticUpdate
    global whichplot

    set e1 1
    set e2 1
    set e3 1
    set e4 1

    set LanguageMaster(FanlabGraphicSettings.title)    "Preferences "
    set LanguageMaster(FanlabGraphicSettings.eixos)    "Axis "
    set LanguageMaster(FanlabGraphicSettings.scale)    "Scale "
    set LanguageMaster(FanlabGraphicSettings.eixX)     "X-axis "
    set LanguageMaster(FanlabGraphicSettings.eixY)     "Y-axis " 
    set LanguageMaster(FanlabGraphicSettings.xyformat) "Format " 
    set LanguageMaster(FanlabGraphicSettings.grid)     "Show lines " 
    set LanguageMaster(FanlabGraphicSettings.line)     "linear "
    set LanguageMaster(FanlabGraphicSettings.loga)     "logarithmic " 
    set LanguageMaster(FanlabGraphicSettings.graf)     "Results "
    set LanguageMaster(FanlabGraphicSettings.styl)     "Style " 
    set LanguageMaster(FanlabGraphicSettings.lini)     "Line " 
    set LanguageMaster(FanlabGraphicSettings.colo)     "Color " 
    set LanguageMaster(FanlabGraphicSettings.punt)     "Points " 
    set LanguageMaster(FanlabGraphicSettings.void)     "Invisible " 
    set LanguageMaster(FanlabGraphicSettings.xmin)     "x min "
    set LanguageMaster(FanlabGraphicSettings.xmax)     "x max " 
    set LanguageMaster(FanlabGraphicSettings.ymin)     "y min "
    set LanguageMaster(FanlabGraphicSettings.ymax)     "y max "
    set LanguageMaster(FanlabGraphicSettings.auto)     "Automatic " 
    set LanguageMaster(FanlabGraphicSettings.fix)      "Fix to : "
    set LanguageMaster(FanlabGraphicSettings.but1)     "Accept " 
    set LanguageMaster(FanlabGraphicSettings.but2)     "Cancel " 
    set LanguageMaster(AcceptGraphicSettings.err1)     "Logarithmic scale is impossible " 

#
# Menus
# 
    menu $W.menu -tearoff 0 
     
# File
    set m $W.menu.files
    menu $m -tearoff 0
    $W.menu add cascade -label "Results"  -menu $m
    
    $m add check -label "Automatic update"    -command "AlyaGraphAccept" -variable AlyaAutomaticUpdate
    $m add separator 
    $m add radio -label "Global convergence"  -command "destroy $W;AlyaConvergence  0" -variable whichplot -value  0
    $m add radio -label "NASTIN convergence"  -command "destroy $W;AlyaConvergence  1" -variable whichplot -value  1
    $m add radio -label "TEMPER convergence"  -command "destroy $W;AlyaConvergence  2" -variable whichplot -value  2
    $m add radio -label "CODIRE convergence"  -command "destroy $W;AlyaConvergence  3" -variable whichplot -value  3
    $m add radio -label "TURBUL convergence"  -command "destroy $W;AlyaConvergence  4" -variable whichplot -value  4
    $m add radio -label "GOTITA convergence"  -command "destroy $W;AlyaConvergence 11" -variable whichplot -value 11
    $m add radio -label "WAVEQU convergence"  -command "destroy $W;AlyaConvergence 12" -variable whichplot -value 12
    $m add radio -label "LEVELS convergence"  -command "destroy $W;AlyaConvergence 13" -variable whichplot -value 13
    $m add separator 
    $m add command -label "Quit"              -command "destroy $W"

#

    set m $W.menu.eixos
    menu $m -tearoff 0
    $W.menu add cascade -label "Axis"  -menu $m

    $m add cascade -label " Scale "     -menu $W.menu.eixos.scale
    $m add cascade -label " Format "    -menu $W.menu.eixos.format
    $m add cascade -label " Labels "    -menu $W.menu.eixos.labels
    $m add cascade -label " Show grid " -menu $W.menu.eixos.grid
    $m add separator 
    $m add command -label " Limits "    -command {AlyaEditLimits} 

#   Scale
    set m $W.menu.eixos.scale
    menu $m  -tearoff 0
    $m add cascade -label " X-axis "    -menu $W.menu.eixos.scale.x
    $m add cascade -label " Y-axis "    -menu $W.menu.eixos.scale.y       
    set m $W.menu.eixos.scale.x
    menu $m -tearoff 0
    $m add radio -label $LanguageMaster(FanlabGraphicSettings.line) -variable eixX -value "Linear"  -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label $LanguageMaster(FanlabGraphicSettings.loga) -variable eixX -value "Logarithmic" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    set m $W.menu.eixos.scale.y 
    menu $m -tearoff 0
    $m add radio -label $LanguageMaster(FanlabGraphicSettings.line) -variable eixY -value "Linear"  -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label $LanguageMaster(FanlabGraphicSettings.loga) -variable eixY -value "Logarithmic" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 

#   Format
    set m $W.menu.eixos.format
    menu $m  -tearoff 0
    $m add cascade -label " X-axis "    -menu $W.menu.eixos.format.fx
    $m add cascade -label " Y-axis "    -menu $W.menu.eixos.format.fy       
    set m $W.menu.eixos.format.fx 
    menu $m -tearoff 0
    $m add radio -label "%3.0f"  -variable formatxX -value "%3.0f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%3.1f"  -variable formatxX -value "%3.1f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%3.2f"  -variable formatxX -value "%3.2f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%4.0f"  -variable formatxX -value "%4.0f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%4.1f"  -variable formatxX -value "%4.1f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%4.2f"  -variable formatxX -value "%4.2f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%5.0f"  -variable formatxX -value "%5.0f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%5.3f"  -variable formatxX -value "%5.3f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%4.1e"  -variable formatxX -value "%4.1e" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%4.2e"  -variable formatxX -value "%4.2e" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%5.2e"  -variable formatxX -value "%5.2e" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%5.3e"  -variable formatxX -value "%5.3e" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    set m $W.menu.eixos.format.fy
    menu $m -tearoff 0
    $m add radio -label "%3.0f"  -variable formatyY -value "%3.0f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%3.1f"  -variable formatyY -value "%3.1f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%3.2f"  -variable formatyY -value "%3.2f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%4.0f"  -variable formatyY -value "%4.0f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%4.1f"  -variable formatyY -value "%4.1f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%4.2f"  -variable formatyY -value "%4.2f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  
    $m add radio -label "%5.0f"  -variable formatyY -value "%5.0f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  
    $m add radio -label "%5.3f"  -variable formatyY -value "%5.3f" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  
    $m add radio -label "%4.1e"  -variable formatyY -value "%4.1e" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%4.2e"  -variable formatyY -value "%4.2e" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%5.2e"  -variable formatyY -value "%5.2e" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    $m add radio -label "%5.3e"  -variable formatyY -value "%5.3e" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"   

#   Labels
    set m $W.menu.eixos.labels
    menu $m  -tearoff 0
    $m add cascade -label " X-axis "    -menu $W.menu.eixos.labels.lx
    $m add cascade -label " Y-axis "    -menu $W.menu.eixos.labels.ly       
    set m $W.menu.eixos.labels.lx 
    menu $m -tearoff 0 
    $m add radio -label "2"  -variable nlabelsx -value "2" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  
    $m add radio -label "4"  -variable nlabelsx -value "4" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  
    $m add radio -label "6"  -variable nlabelsx -value "6" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  
    $m add radio -label "8"  -variable nlabelsx -value "8" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  
    $m add radio -label "10" -variable nlabelsx -value "10" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  
    set m $W.menu.eixos.labels.ly 
    menu $m -tearoff 0 
    $m add radio -label "2"  -variable nlabelsy -value "2" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  
    $m add radio -label "4"  -variable nlabelsy -value "4" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  
    $m add radio -label "6"  -variable nlabelsy -value "6" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  
    $m add radio -label "8"  -variable nlabelsy -value "8" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  
    $m add radio -label "10" -variable nlabelsy -value "10" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  

#   Show Grid
    set m $W.menu.eixos.grid
    menu $m  -tearoff 0
    $m add cascade -label " X-axis "    -menu $W.menu.eixos.grid.gx
    $m add cascade -label " Y-axis "    -menu $W.menu.eixos.grid.gy       
    set m $W.menu.eixos.grid.gx 
    menu $m -tearoff 0 
    $m add radio -label "Yes" -variable showgridX -value "Yes" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  
    $m add radio -label "No " -variable showgridX -value "No"  -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  
    set m $W.menu.eixos.grid.gy 
    menu $m -tearoff 0 
    $m add radio -label "Yes" -variable showgridY -value "Yes" -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  
    $m add radio -label "No " -variable showgridY -value "No"  -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4"  

 
    set m $W.menu.grafic
    menu $m -tearoff 0
    $W.menu add cascade -label "Options" -menu $m
    set k 1
    for {set i 0} {$i < $ngraph} {incr i 1} {
	$m add cascade -label " $messagegraphic($k) " -command {} -menu $W.menu.grafic.var$i
	#if { $i < [expr $ngraph-1] } { $m add separator }
	incr k 1
    }
    set k 1
    for {set i 0} {$i < $ngraph} {incr i 1} {
	set m $W.menu.grafic.var$i
	menu $m  -tearoff 0
	#$m add command -label "$messagegraphic($k)" -command {} 
	$m add cascade -label " Style " -menu $W.menu.grafic.var$i.l$i
	$m add cascade -label " Color " -menu $W.menu.grafic.var$i.c$i
	#if { $i < [expr $ngraph-1] } { $m add separator }
	incr k 1
    }
    for {set i 0} {$i < $ngraph} {incr i 1} {
	set m $W.menu.grafic.var$i.l$i
	menu $m -tearoff 0
	$m add radio -label $LanguageMaster(FanlabGraphicSettings.lini) -variable choiceplot($i) \
                                    -value "Line"  -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
	$m add radio -label $LanguageMaster(FanlabGraphicSettings.punt) -variable choiceplot($i) \
                                    -value "Points"  -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
	$m add radio -label $LanguageMaster(FanlabGraphicSettings.void) -variable choiceplot($i) \
                                    -value "No"  -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
	set m $W.menu.grafic.var$i.c$i
	menu $m -tearoff 0
	$m add radio -label "" -background black  -variable choicecolor($i) -value "black" \
                               -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
	$m add radio -label "" -background white  -variable choicecolor($i) -value "white" \
                               -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
	$m add radio -label "" -background red    -variable choicecolor($i) -value "red" \
                               -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
	$m add radio -label "" -background blue   -variable choicecolor($i) -value "blue" \
                               -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
	$m add radio -label "" -background green  -variable choicecolor($i) -value "green"  \
                               -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4 " 
	$m add radio -label "" -background yellow -variable choicecolor($i) -value "yellow" \
                               -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
	$m add radio -label "" -background orange -variable choicecolor($i) -value "orange" \
                               -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
	$m add radio -label "" -background cyan1   -variable choicecolor($i) -value "cyan1" \
                                    -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
	$m add radio -label "" -background purple1 -variable choicecolor($i) -value "purple1" \
                                    -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
	$m add radio -label "" -background RosyBrown1 -variable choicecolor($i) -value "RosyBrown1" \
                                    -command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    }

# menu.eixos.style

    set m $W.menu.style
    menu $m -tearoff 0
    $W.menu add cascade -label "Style"  -menu $m
    $m add cascade -label " Background Color "     -menu $W.menu.style.backcolor
    $m add cascade -label " Foreground Color "     -menu $W.menu.style.forecolor

    set m $W.menu.style.backcolor
    menu $m  -tearoff 0   
    $m add radio -label "" -foreground black  -background black -variable graphcolor -value "black" \
	    -command "AlyaGraphsChangeColor" 
    $m add radio -label "" -foreground white  -background white -variable graphcolor -value "white" \
	    -command "AlyaGraphsChangeColor" 
    $m add radio -label "" -foreground red    -background red   -variable graphcolor -value "red" \
	    -command "AlyaGraphsChangeColor" 
    $m add radio -label "" -foreground blue   -background blue  -variable graphcolor -value "blue" \
	    -command "AlyaGraphsChangeColor" 
    $m add radio -label "" -foreground green  -background green -variable graphcolor -value "green"  \
	    -command "AlyaGraphsChangeColor" 
    $m add radio -label "" -foreground yellow -background yellow -variable graphcolor -value "yellow" \
	    -command "AlyaGraphsChangeColor" 
    $m add radio -label "" -foreground orange -background orange -variable graphcolor -value "orange" \
	    -command "AlyaGraphsChangeColor" 
    $m add radio -label "" -foreground grey   -background grey   -variable graphcolor -value "grey" \
	    -command "AlyaGraphsChangeColor" 
    $m add radio -label "" -foreground cyan1  -background cyan1  -variable graphcolor -value "cyan1" \
	    -command "AlyaGraphsChangeColor" 
    $m add radio -label "" -foreground purple1   -background purple1 -variable graphcolor -value "purple1" \
	    -command "AlyaGraphsChangeColor" 
    $m add radio -label "" -foreground RosyBrown1   -background RosyBrown1  -variable graphcolor -value "RosyBrown1" \
	    -command "AlyaGraphsChangeColor" 
    
    set m $W.menu.style.forecolor
    menu $m  -tearoff 0   
    $m add radio -label "" -background black  -variable forecolor -value "black" \
	    -command "AlyaFontChangeColor" -font {-weight bold -size 8}
    $m add radio -label "" -background white  -variable forecolor -value "white" \
	    -command "AlyaFontChangeColor" -font {-weight bold -size 8}
    $m add radio -label "" -background red    -variable forecolor -value "red" \
	    -command "AlyaFontChangeColor" 	-font {-weight bold -size 8}	     
    $m add radio -label "" -background blue   -variable forecolor -value "blue" \
	    -command "AlyaFontChangeColor" 	-font {-weight bold -size 8}	     
    $m add radio -label "" -background green  -variable forecolor -value "green"  \
	    -command "AlyaFontChangeColor" 	-font {-weight bold -size 8}	     
    $m add radio -label "" -background yellow -variable forecolor -value "yellow" \
	    -command "AlyaFontChangeColor" 	-font {-weight bold -size 8}	     
    $m add radio -label "" -background orange -variable forecolor -value "orange" \
	    -command "AlyaFontChangeColor" 	-font {-weight bold -size 8}	     
    $m add radio -label "" -background grey   -variable forecolor -value "grey" \
	    -command "AlyaFontChangeColor" 
    $m add radio -label "" -background cyan1   -variable forecolor -value "cyan1" \
	    -command "AlyaFontChangeColor" 
    $m add radio -label "" -background purple1   -variable forecolor -value "purple1" \
	    -command "AlyaFontChangeColor" 
    $m add radio -label "" -background RosyBrown1   -variable forecolor -value "RosyBrown1" \
	    -command "AlyaFontChangeColor" 
    
    
    $W configure -menu $W.menu

#    tk_messageBox -icon error  -type ok -message $xmin
#
# Limits.
#
    #frame  $W.limts
# Eix x
    #frame  $W.limts.x -relief groove -bd 2 
# xmin
    #frame       $W.limts.x.min 
    #set m       $W.limts.x.min 
    #label       $m.txt -text $LanguageMaster(FanlabGraphicSettings.xmin) -font {-weight bold -size 8}
    #pack        $m.txt -side top 
    #frame       $m.radio
    #radiobutton $m.radio.rb1 -text $LanguageMaster(FanlabGraphicSettings.auto) -variable xminX -value "Automatic" 
    #pack        $m.radio.rb1 -side top -anchor w
    #frame       $m.radio.rb2
    #radiobutton $m.radio.rb2.b -text $LanguageMaster(FanlabGraphicSettings.fix) -variable xminX -value "Fixed" 
    #pack        $m.radio.rb2.b -side left -anchor n
    #set e1 [entry $m.radio.rb2.e -width 8 -textvariable valuee1]
    #pack        $m.radio.rb2.e -side right -anchor n
    #pack        $m.radio.rb2 -side bottom -anchor w
    #pack        $m.radio     -side top -anchor w
    #pack        $m -side left -anchor n 
# xmax
    #frame       $W.limts.x.max 
    #set m       $W.limts.x.max 
    #label       $m.txt -text $LanguageMaster(FanlabGraphicSettings.xmax) -font {-weight bold -size 8}
    #pack        $m.txt -side top 
    #frame       $m.radio
    #radiobutton $m.radio.rb1 -text $LanguageMaster(FanlabGraphicSettings.auto) -variable xmaxX \
    #            -value "Automatic"
    #pack        $m.radio.rb1 -side top -anchor w
    #frame       $m.radio.rb2
    #radiobutton $m.radio.rb2.b -text $LanguageMaster(FanlabGraphicSettings.fix) -variable xmaxX \
    #            -value "Fixed" 
    #pack        $m.radio.rb2.b -side left -anchor n
    #set e2 [entry $m.radio.rb2.e -width 8 -textvariable valuee2]
    #pack        $m.radio.rb2.e -side right -anchor n
    #pack        $m.radio.rb2 -side bottom -anchor w
    #pack        $m.radio     -side top -anchor w
    #pack        $m -side right -anchor n -ipadx 2m -ipady 2m
    #pack        $W.limts.x -fill both -padx 2m -pady 2m -side left
# Eix y
    #frame       $W.limts.y -relief groove -bd 2 
# ymin
    #frame       $W.limts.y.min 
    #set m       $W.limts.y.min 
    #label       $m.txt -text $LanguageMaster(FanlabGraphicSettings.ymin) -font {-weight bold -size 8}
    #pack        $m.txt -side top 
    #frame       $m.radio
    #radiobutton $m.radio.rb1 -text $LanguageMaster(FanlabGraphicSettings.auto) -variable yminY -value "Automatic" 
    #pack        $m.radio.rb1 -side top -anchor w
    #frame       $m.radio.rb2
    #radiobutton $m.radio.rb2.b -text $LanguageMaster(FanlabGraphicSettings.fix) -variable yminY -value "Fixed"
    #pack        $m.radio.rb2.b -side left -anchor n
    #set e3 [entry $m.radio.rb2.e -width 8 -textvariable valuee3]
    #pack        $m.radio.rb2.e -side right -anchor n
    #pack        $m.radio.rb2 -side bottom -anchor w
    #pack        $m.radio     -side top -anchor w
    #pack        $m -side left -anchor n 
# ymax
    #frame       $W.limts.y.max 
    #set m       $W.limts.y.max 
    #label       $m.txt -text $LanguageMaster(FanlabGraphicSettings.ymax) -font {-weight bold -size 8}
    #pack        $m.txt -side top 
    #frame       $m.radio
    #radiobutton $m.radio.rb1 -text $LanguageMaster(FanlabGraphicSettings.auto) -variable ymaxY -value "Automatic"
    #pack        $m.radio.rb1 -side top -anchor w
    #frame       $m.radio.rb2
    #radiobutton $m.radio.rb2.b -text $LanguageMaster(FanlabGraphicSettings.fix) -variable ymaxY -value "Fixed" 
    #pack        $m.radio.rb2.b -side left -anchor n
    #set e4 [entry $m.radio.rb2.e -width 8 -textvariable valuee4]
    #pack        $m.radio.rb2.e -side right -anchor n
    #pack        $m.radio.rb2 -side bottom -anchor w
    #pack        $m.radio     -side top -anchor w
    #pack        $m -side right -anchor n -ipadx 2m -ipady 2m
    #pack        $W.limts.y -fill both -padx 2m -pady 2m -side left
    #
    #pack        $W.limts -side top
#   #
# Bottom frame with buttons.
# 
    #frame  $W.buttons 
    #button $W.buttons.1 -text $LanguageMaster(FanlabGraphicSettings.but1) -relief raised \
	#-command "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    #pack   $W.buttons.1 -side left 
    #button $W.buttons.2 -text $LanguageMaster(FanlabGraphicSettings.but2) -relief raised -font {-weight bold -size 8} \
	#-command "destroy $W"
    #pack   $W.buttons.2 -side left 
    #pack   $W.buttons   -expand true -ipadx 2 -ipady 2 -padx 2 -pady 2

    #set m $W.menu
    #bind $m   <FocusIn> {}
    #bind $m   <ButtonPress> "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    #for {set i 0} {$i < $ngraph} {incr i 1} {
	#set m $W.menu.grafic.l$i    
	#bind $m <Button>  "AcceptGraphicSettings $W $e1 $e2 $e3 $e4" 
    #}
}

proc AcceptGraphicSettings { W e1 e2 e3 e4 } {  
    global messagegraphic messagex messagey nlabelsx nlabelsy \
	ticksx ticksy formatx formaty plotgraphic showgrid color \
	typeaxisX typeaxisY xmax xmin ymax ymin 
    global xpoints ypoints ngraph npoingraph 
    global eixX eixY choiceplot choicecolor 
    global xmaxX xminX yminY ymaxY 
    global formatxX formatyY showgridX showgridY 
    global LanguageMaster   
    global valuee1 valuee2 valuee3 valuee4


#
# Checks is logarithminc scale is possible.
#
    if { $eixX != "Linear" } { 
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		if { "$xpoints($i,$k)" <= "0.0" } { 
		    AlyaError "Logarithmic scale is impossible"
		    return
		}
	    }
	}
    }
    if { $eixY != "Linear" } { 
	for {set k 0} {$k < $ngraph} {incr k 1} {
	    for {set i 0} {$i < $npoingraph($k)} {incr i 1} {
		if { "$ypoints($i,$k)" <= "0.0" } { 
		    AlyaError "Logarithmic scale is impossible"
		    return
		}
	    }
	}
    }
    
    set typeaxisX      $eixX
    set typeaxisY      $eixY
    for {set i 0} {$i < $ngraph} {incr i 1} {
	set plotgraphic($i) $choiceplot($i)  
	set color($i)       $choicecolor($i) 
    } 
#
# Axis format.
# 
    set formatx $formatxX
    set formaty $formatyY
#
# Grid.
# 
    if { $showgridX == "Yes" } { 
	if { $showgridY == "Yes" } { 
	    set showgrid "Both" 
	} else { 
	    set showgrid "X Axis" 
	} 
    } else { 
	if { $showgridY == "Yes" } {  
	    set showgrid "Y Axis" 
	} else { 
	    set showgrid "None" 
	} 
    }  
#
# Actualizes changes.
#
    SetGraphic
    ShowGraphic .gid.postwin.plot.graphic
}

proc AlyaEditLimits { } {
    global ticksx ticksy formatx formaty plotgraphic showgrid color \
	    typeaxisX typeaxisY xmax xmin ymax ymin 
    global xpoints ypoints ngraph npoingraph 
    global eixX eixY choiceplot choicecolor 
    global xmaxX xminX yminY ymaxY 
    global formatxX formatyY showgridX showgridY 
    global Color
    global valuee1 valuee2 valuee3 valuee4
    global e1 e2 e3 e4

    set valuee1 $xmin 
    set valuee2 $xmax 
    set valuee3 $ymin 
    set valuee4 $ymax 
    
    set LanguageMaster(FanlabGraphicSettings.title)    "Preferences "
    set LanguageMaster(FanlabGraphicSettings.eixos)    "Axis "
    set LanguageMaster(FanlabGraphicSettings.scale)    "Scale "
    set LanguageMaster(FanlabGraphicSettings.eixX)     "X-axis "
    set LanguageMaster(FanlabGraphicSettings.eixY)     "Y-axis " 
    set LanguageMaster(FanlabGraphicSettings.xyformat) "Format " 
    set LanguageMaster(FanlabGraphicSettings.grid)     "Show lines " 
    set LanguageMaster(FanlabGraphicSettings.line)     "linear "
    set LanguageMaster(FanlabGraphicSettings.loga)     "logarithmic " 
    set LanguageMaster(FanlabGraphicSettings.graf)     "Graphs "
    set LanguageMaster(FanlabGraphicSettings.styl)     "Style " 
    set LanguageMaster(FanlabGraphicSettings.lini)     "Line " 
    set LanguageMaster(FanlabGraphicSettings.colo)     "Color " 
    set LanguageMaster(FanlabGraphicSettings.punt)     "Points " 
    set LanguageMaster(FanlabGraphicSettings.void)     "Invisible " 
    set LanguageMaster(FanlabGraphicSettings.xmin)     "x min "
    set LanguageMaster(FanlabGraphicSettings.xmax)     "x max " 
    set LanguageMaster(FanlabGraphicSettings.ymin)     "y min "
    set LanguageMaster(FanlabGraphicSettings.ymax)     "y max "
    set LanguageMaster(FanlabGraphicSettings.auto)     "Automatic " 
    set LanguageMaster(FanlabGraphicSettings.fix)      "Fix to : "
    set LanguageMaster(FanlabGraphicSettings.but1)     "Accept " 
    set LanguageMaster(FanlabGraphicSettings.but2)     "Cancel " 
    set LanguageMaster(AcceptGraphicSettings.err1)     "Logarithmic scale is impossible " 

    # 
    # Creates canvas for graphic. Canvas window is always .postwin.plot.graphic
    # 
    set W .gid.editlimits
    AlyaCreateGeneralWindow $W "Alya - Edit graph limits" 100 100 .gid.postwin
#
# Limits.
#
    frame       $W.all   -bg $Color(border) 

    frame       $W.all.title     -bg $Color(border) 
    label       $W.all.title.txt -text "Set the minimum and maximum axis values" -font {-weight bold -size 8} -bg $Color(border) 
    pack        $W.all.title.txt -side top -anchor w
    pack        $W.all.title     -side top -padx 2 -pady 2

    frame       $W.all.limts -bg $Color(border) 

# Eix x
    frame       $W.all.limts.x -relief groove -bd 2 
# xmin
    frame       $W.all.limts.x.min 
    set m       $W.all.limts.x.min 
    label       $m.txt -text $LanguageMaster(FanlabGraphicSettings.xmin) -font {-weight bold -size 8}
    pack        $m.txt -side top
    frame       $m.radio
    radiobutton $m.radio.rb1 -text $LanguageMaster(FanlabGraphicSettings.auto) -variable xminX -value "Automatic" 
    pack        $m.radio.rb1 -side top -anchor w
    frame       $m.radio.rb2
    radiobutton $m.radio.rb2.b -text $LanguageMaster(FanlabGraphicSettings.fix) -variable xminX -value "Fixed" 
    pack        $m.radio.rb2.b -side left -anchor n
    set e1 [entry $m.radio.rb2.e -width 8 -textvariable valuee1]
    pack        $m.radio.rb2.e -side right -anchor n
    pack        $m.radio.rb2 -side bottom -anchor w
    pack        $m.radio     -side top -anchor w
    pack        $m -side left -anchor n 
# xmax
    frame       $W.all.limts.x.max 
    set m       $W.all.limts.x.max 
    label       $m.txt -text $LanguageMaster(FanlabGraphicSettings.xmax) -font {-weight bold -size 8}
    pack        $m.txt -side top 
    frame       $m.radio
    radiobutton $m.radio.rb1 -text $LanguageMaster(FanlabGraphicSettings.auto) -variable xmaxX \
                -value "Automatic"
    pack        $m.radio.rb1 -side top -anchor w
    frame       $m.radio.rb2
    radiobutton $m.radio.rb2.b -text $LanguageMaster(FanlabGraphicSettings.fix) -variable xmaxX \
                -value "Fixed" 
    pack        $m.radio.rb2.b -side left -anchor n
    set e2 [entry $m.radio.rb2.e -width 8 -textvariable valuee2]
    pack        $m.radio.rb2.e -side right -anchor n
    pack        $m.radio.rb2 -side bottom -anchor w
    pack        $m.radio     -side top -anchor w
    pack        $m -side right -anchor n -ipadx 2m -ipady 2m
    pack        $W.all.limts.x -fill both -padx 2m -pady 2m 
# Eix y
    frame       $W.all.limts.y -relief groove -bd 2 
# ymin
    frame       $W.all.limts.y.min 
    set m       $W.all.limts.y.min 
    label       $m.txt -text $LanguageMaster(FanlabGraphicSettings.ymin) -font {-weight bold -size 8}
    pack        $m.txt -side top 
    frame       $m.radio
    radiobutton $m.radio.rb1 -text $LanguageMaster(FanlabGraphicSettings.auto) -variable yminY -value "Automatic" 
    pack        $m.radio.rb1 -side top -anchor w
    frame       $m.radio.rb2
    radiobutton $m.radio.rb2.b -text $LanguageMaster(FanlabGraphicSettings.fix) -variable yminY -value "Fixed"
    pack        $m.radio.rb2.b -side left -anchor n
    set e3 [entry $m.radio.rb2.e -width 8 -textvariable valuee3]
    pack        $m.radio.rb2.e -side right -anchor n
    pack        $m.radio.rb2 -side bottom -anchor w
    pack        $m.radio     -side top -anchor w
    pack        $m -side left -anchor n 
# ymax
    frame       $W.all.limts.y.max 
    set m       $W.all.limts.y.max 
    label       $m.txt -text $LanguageMaster(FanlabGraphicSettings.ymax) -font {-weight bold -size 8}
    pack        $m.txt -side top 
    frame       $m.radio
    radiobutton $m.radio.rb1 -text $LanguageMaster(FanlabGraphicSettings.auto) -variable ymaxY -value "Automatic"
    pack        $m.radio.rb1 -side top -anchor w
    frame       $m.radio.rb2
    radiobutton $m.radio.rb2.b -text $LanguageMaster(FanlabGraphicSettings.fix) -variable ymaxY -value "Fixed" 
    pack        $m.radio.rb2.b -side left -anchor n
    set e4 [entry $m.radio.rb2.e -width 8 -textvariable valuee4]
    pack        $m.radio.rb2.e -side right -anchor n
    pack        $m.radio.rb2 -side bottom -anchor w
    pack        $m.radio     -side top -anchor w
    pack        $m -side right -anchor n -ipadx 2m -ipady 2m
    pack        $W.all.limts.y -fill both -padx 2m -pady 2m 
    
    pack        $W.all.limts -side top

    frame       $W.all.buttons -bg $Color(border) 
    button      $W.all.buttons.1 -text "Accept" -relief raised \
	                     -command "AlyaSetAcceptMinMax $W $e1 $e2 $e3 $e4" -width 10
    pack        $W.all.buttons.1 -side left    -padx 1m
    button      $W.all.buttons.2 -text "Close" -relief raised -command "destroy $W" -width 10
    pack        $W.all.buttons.2 -side left    -padx 1m
    pack        $W.all.buttons   -side bottom  -padx 2m -pady 2m 

    pack        $W.all

    bind        $W <Return> "tkButtonInvoke $W.all.buttons.1"
    bind        $W <Escape> "tkButtonInvoke $W.all.buttons.2"
    
}

proc AlyaSetAcceptMinMax { W e1 e2 e3 e4 } {
    global ticksx ticksy formatx formaty plotgraphic showgrid color \
	    typeaxisX typeaxisY xmax xmin ymax ymin 
    global xpoints ypoints ngraph npoingraph 
    global eixX eixY choiceplot choicecolor 
    global xmaxX xminX yminY ymaxY 
    global formatxX formatyY showgridX showgridY 
    global Color
    global valuee1 valuee2 valuee3 valuee4
    global AlyaAutomaticUpdate
#
#  Maximum and minimum values.
#
    set patron {^((-)?[0-9]+(.[0-9]+)?)$}
    if { $xminX == "Fixed" } { 
	if { [ catch { set xmin [expr double([$e1 get])] } ] } {
	    AlyaError "Limits must be real numbers"
	    return
	}
	if { [regexp $patron $xmin] != 1 } {
	    AlyaError "Limits must be real numbers"
	    return
	}
    }
    if { $xmaxX == "Fixed" } { 
	if { [ catch { set xmax [expr double([$e2 get])] } ] } {
	    AlyaError "Limits must be real numbers"
	    return
	}
	if { [ catch {expr $xmax} ] } {
	    AlyaError "Limits must be real numbers"
	    return
	}
    }
    if { $yminY == "Fixed" } { 
	if { [ catch { set ymin [expr double([$e3 get])] } ] } {
	    AlyaError "Limits must be real numbers"
	    return
	}
	if { [regexp $patron $ymin] != 1 } {
	    AlyaError "Limits must be real numbers"
	    return
	}
    }
    if { $ymaxY == "Fixed" } { 
	if { [ catch { set ymax [expr double([$e4 get])] } ] } {
	    AlyaError "Limits must be real numbers"
	    return
	}
	if { [regexp $patron $ymax] != 1 } {
	    AlyaError "Limits must be real numbers"
	    return
	}
    }
    SetGraphic
    ShowGraphic .gid.postwin.plot.graphic
    if { [winfo exists .gid.postwin.plot.graphic]== 1 && $AlyaAutomaticUpdate == 1} { 
	after 8000 AlyaRefreshWindow .gid.postwin.plot.graphic
    }
    #AlyaRefreshWindow .gid.postwin.plot.graphic 
    
}

proc AlyaGraphsChangeColor { } {
    global  graphcolor   
    .gid.postwin.plot.graphic conf  -bg $graphcolor
    .gid.postwin.plot.header  conf  -bg $graphcolor
    .gid.postwin.plot         conf  -bg $graphcolor
}

proc AlyaFontChangeColor { } {
    global  forecolor   
    .gid.postwin.plot.header  conf  -fg $forecolor
    AlyaRefreshWindow .gid.postwin.plot.graphic 
}

proc AlyaCanvasZoom {canvas scalx scaly} {
    $canvas scale all 0 0 $scalx $scaly
    if { $scalx > $scaly } {
	set scalm $scaly
    } else {
	set scalm $scalx	
    }
    foreach item [$canvas find all] {
	if {[$canvas type $item] == "text"} {
	    set font [font actual [$canvas itemcget $item -font]]
	    set index [lsearch -exact $font -size]
	    incr index
	    set size [lindex $font $index]
	    set size [expr {round($size * $scalm)}]
	    set font [lreplace $font $index $index $size]
	    $canvas itemconfigure $item -font $font
	}
    }
    
    $canvas configure -scrollregion [$canvas bbox all]
}

proc AlyaAutomatic { W } {
    global xmaxX xminX yminY ymaxY
    set xminX  "Automatic"
    set xmaxX  "Automatic"
    set yminY  "Automatic"
    set ymaxY  "Automatic"
    SetGraphic
    ShowGraphic $W  
}

proc AlyaRefresh { W } {
    global xmaxX xminX yminY ymaxY
    #set xminX  "Automatic"
    #set xmaxX  "Automatic"
    #set yminY  "Automatic"
    #set ymaxY  "Automatic"
    #SetGraphic
    AlyaRefreshWindow $W    
}

proc AlyaGraphAccept { } {
    #set W .gid.graphpostscript
    #AlyaCreateGeneralWindow $W "Alya - Postscript" 100 100 .gid
    #Setup $W
    #return
    #Setup .gid.postwin.plot.graphic "C:/nada.ps"
    #return
    global AlyaAutomaticUpdate
    ShowGraphic .gid.postwin.plot.graphic
    if { [winfo exists .gid.postwin.plot.graphic]== 1 && $AlyaAutomaticUpdate == 1} { 
	after 8000 AlyaRefreshWindow .gid.postwin.plot.graphic
    }
}

proc AlyaGraphConfigure { } {
    global old_width
    global old_height

    set wwidth  [winfo width  .gid.postwin]
    set wheight [winfo height .gid.postwin]
    #WarnWinText $wwidth,$wheight
    #return
    set scalx [expr $wwidth/$old_width]
    set scaly [expr $wheight/$old_height]
    WarnWinText $scalx,$scaly
    AlyaCanvasZoom .gid.postwin.plot.graphic $scalx $scalx
    set old_width  $wwidth
    set old_height $wheight
    #set old_width  [winfo width  .gid.postwin]
    #set old_height [winfo height .gid.postwin]

}

proc AlyaGraphsMark { x y w } {
    global state
    # Find the object
    set state($w,obj) [$w find closest $x $y]
    set state($w,x) $x
    set state($w,y) $y
}

proc AlyaGraphsDrag { x y w } {
    global state
    set dx [expr $x - $state($w,x)]
    set dy [expr $y - $state($w,y)]
    $w move $state($w,obj) $dx $dy
    set state($w,x) $x
    set state($w,y) $y
}

proc Setup { w file } {
    global fontMap

    # Create several strings in different fonts and sizes
    foreach family {times courier helvetica "MS Sans Serif"} {
	set weight bold
	switch -- $family {
	    times { set fill blue; set psfont Times}
	    courier { set fill green; set psfont Courier }
	    helvetica { set fill red; set psfont Helvetica }
	    "MS Sans Serif" { set fill red; set psfont Times }
	}
	foreach size {8 10 14 24} {
	    # Guard against missing fonts
	    set fontMap(-*-$family-$weight-*-*-*-$size-*) [list $psfont $size]
	}
    }
    set size 24
    set family times
    WarnWinText $fontMap(-*-$family-$weight-*-*-*-$size-*)
    set fontMap(fixed) [list Courier 24]
    set colorMap(blue)  {0.1 0.1 0.9 setrgbcolor}
    set colorMap(green) {0.0 0.9 0.1 setrgbcolor}
    $w postscript -fontmap fontMap(-*-$family-$weight-*-*-*-$size-*) -colormap colorMap \
	    -file $file \
	    -pagex 0.i -pagey 11.i -pageanchor nw
}

proc Postscript { c file } {
    global fontMap
    # Tweak the output color
    set colorMap(blue) {0.1 0.1 0.9 setrgbcolor}
    set colorMap(green) {0.0 0.9 0.1 setrgbcolor}
    # Position the text at the upper-left corner of
    # an 8.5 by 11 inch sheet of paper
    $c postscript -fontmap fontMap -colormap colorMap \
	    -file $file \
	    -pagex 0.i -pagey 11.i -pageanchor nw
}

proc AlyaReadForce { } {  
    global ListForce FactorForce
    global Alya
    global number
    global messagegraphic messagexaxis messageyaxis nlabelsx nlabelsy \
	   ticksx ticksy formatx formaty plotgraphic showgrid color \
	   typeaxisX typeaxisY xmax xmin ymax ymin ncolu messagex messagey
    global xpoints ypoints npoingraph ngraph
    global whichplot 
    global Text

    #set Force [split $ListForce " "]
    set Force $ListForce
    set ForceNumber [llength $Force]
    for {set k 0} {$k < $ForceNumber} {incr k 1} {
	set Force [lreplace $Force $k $k [string trim [lindex $Force $k] " "]]
	set ForceLength [string length [lindex $Force $k]]
	if { $ForceLength == 1 } {
	    set CurrentForce "000000"
	} elseif { $ForceLength == 2 } {
	    set CurrentForce "00000"
	} elseif { $ForceLength == 3 } {
	    set CurrentForce "0000"
	} elseif { $ForceLength == 4 } {
	    set CurrentForce "000"
	} elseif { $ForceLength == 5 } {
	    set CurrentForce "00"
	} elseif { $ForceLength == 6 } {
	    set CurrentForce "0"
	} elseif { $ForceLength == 7 } {
	    set CurrentForce ""
	}
	append CurrentForce [lindex $Force $k]
	set Force [lreplace $Force $k $k $CurrentForce]
    }

    # Open convergence file
    set extension ".frc"
    AlyaGetPaths
    for {set k 0} {$k < $ForceNumber} {incr k 1} {
	set filepath [file join $Alya(resultspath) $Alya(projectname)-[lindex $Force $k]$extension]
	if {[file exists $filepath]!=1} {
	    for {set j 0} {$j < $ngraph} {incr j 1} { set npoingraph($j) 0.0 }
	    return
	} else {
	    set luout [open $filepath r 0600]
	} 
	set nline 0
	if { $k == 0 } {
	    while { ![eof $luout] } {
		set linea  [gets $luout]
		set ndime  [NumberDecode $linea]
		if { $ndime > "1" } { 
		    set nlin1  [expr $nline+1]
		    for {set j 0} {$j < $ngraph} {incr j 1} {
			set xpoints($nline,$j) $number(2)
			set ypoints($nline,$j) [expr $FactorForce*$number($ncolu($j))]
		    }
		    incr nline 1
		}   
	    }
	} else {
	    while { ![eof $luout] } {
		set linea  [gets $luout]
		set ndime  [NumberDecode $linea]
		if { $ndime > "1" } { 
		    set nlin1  [expr $nline+1]
		    for {set j 0} {$j < $ngraph} {incr j 1} {
			set xpoints($nline,$j) $number(2)
			set ypoints($nline,$j) [expr $ypoints($nline,$j)+$FactorForce*$number($ncolu($j))]
		    }
		    incr nline 1
		}   
	    }
	}
	for {set j 0} {$j < $ngraph} {incr j 1} {    
	    set npoingraph($j) $nline
	}
	close $luout
    }
}


proc AlyaReadACH { } {  
    global Alya
    global number
    global messagegraphic messagexaxis messageyaxis nlabelsx nlabelsy \
           ticksx ticksy formatx formaty plotgraphic showgrid color \
           typeaxisX typeaxisY xmax xmin ymax ymin ncolu messagex messagey
    global xpoints ypoints npoingraph ngraph
    global whichplot 
    global Text
    global GraphNumberZones
    global GraphNameZones

    #  Opens witness file
    set extension ".ach"
    AlyaGetPaths

    set filepath [file join $Alya(resultspath) $Alya(projectname)$extension]
    if {[file exists $filepath]!=1} {
	for {set j 0} {$j < $ngraph} {incr j 1} { set npoingraph($j) 0.0 }
	return
    } else {
	set luout [open $filepath r 0600]
    } 
    # Read zone names
    set nline 0
    set linea  [gets $luout]
    set linea  [gets $luout]
    set linea  [gets $luout]
    set linea  [gets $luout]
    while { ![eof $luout] } {
       set linea  [gets $luout]
       set ndime  [NumberDecode $linea]
       if { $ndime > "1" } { 
         for {set j 0} {$j < $ngraph} {incr j 1} {
             set xpoints($nline,$j) $number(1)
             set ypoints($nline,$j) $number($ncolu($j))
         }
  	 incr nline 1
       }   
    }
    for {set j 0} {$j < $ngraph} {incr j 1} {    
       set npoingraph($j) $nline
    }
    # Close file
    close $luout
}

proc AlyaReadZonesACH { } {  
    global Alya
    global number
    global messagegraphic messagexaxis messageyaxis nlabelsx nlabelsy \
           ticksx ticksy formatx formaty plotgraphic showgrid color \
           typeaxisX typeaxisY xmax xmin ymax ymin ncolu messagex messagey
    global xpoints ypoints npoingraph ngraph
    global whichplot 
    global Text
    global GraphNumberZones
    global GraphNameZones

    #  Opens witness file
    set extension ".ach"
    AlyaGetPaths
    set GraphNumberZones 0
    set GraphNumberZones 0

    set filepath [file join $Alya(resultspath) $Alya(projectname)$extension]
    if {[file exists $filepath]!=1} {
	for {set j 0} {$j < $ngraph} {incr j 1} { set npoingraph($j) 0.0 }
	return
    } else {
	set luout [open $filepath r 0600]
    } 
    # Read zone names
    set nline 0
    set linea  [gets $luout]
    set linea  [gets $luout]
    set linea  [gets $luout]
    set linea  [string trim $linea "$"]
    set GraphNumberZones [NumberDecode $linea]
    for {set Izone 1} {$Izone <= $GraphNumberZones} {incr Izone} {
	set GraphNameZones($Izone) $number($Izone)
    }
    close $luout
    return
}
