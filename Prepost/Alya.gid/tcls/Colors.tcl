proc AlyaColor { } {
    global Color
    #border of the window
    set    Color(border)        "#EFEBE7"   
    #background of the window    
    set    Color(background)    "#D6D3CE"      
    #background of the entries
    set    Color(entry)         white
    #color del text     
    set    Color(text)          black   
    #color etapa corriente
    set    Color(highlight)     "#E8DB90"  
    #set    Color(highlight)     "#E7F7B5"  
    #boton importante
    set    Color(boton1)        "#D6D3CE"   
    set    Color(boton1fg)      "black"   
    #boton general
    set    Color(boton2)        "#D6D3CE"   

}

##########################################################################################
### CCGetRGB CCColorActivo CCColorSombra -> para ContourColores
##########################################################################################

proc AlyaCCGetRGB { w color} {
    set ret $color
    set n [ scan $color #%2x%2x%2x r g b]
    if { $n != 3} {
	set rgb [ winfo rgb $w $color]
	set r [ expr int( 0.5 + [ lindex $rgb 0]/256.0)]
	set g [ expr int( 0.5 + [ lindex $rgb 1]/256.0)]
	set b [ expr int( 0.5 + [ lindex $rgb 2]/256.0)]
	set ret [ format #%2x%2x%2x $r $g $b]
    }
    return $ret
}

proc AlyaCCColorActivo { color_usuario { factor 17} } {
    set ret ""
    set color_nuevo [ AlyaCCGetRGB . $color_usuario]
    set n [ scan $color_nuevo #%2x%2x%2x r g b]
    if { $n == 3} {
	set r [ expr $r + $factor]
	if { $r > 255} { set r 255}
	set g [ expr $g + $factor]
	if { $g > 255} { set g 255}
	set b [ expr $b + $factor]
	if { $b > 255} { set b 255}
	set ret [ format #%2x%2x%2x $r $g $b]
    }
    return $ret
}

proc AlyaCCColorSombra { color_usuario { factor 17} } {
    set ret ""
    set color_nuevo [ AlyaCCGetRGB . $color_usuario]
    set n [ scan $color_nuevo #%2x%2x%2x r g b]
    if { $n == 3} {
	set r [ expr $r - $factor]
	if { $r < 0} { set r 0}
	set g [ expr $g - $factor]
	if { $g < 0} { set g 0}
	set b [ expr $b - $factor]
	if { $b < 0} { set b 0}
	set ret [ format #%2x%2x%2x $r $g $b]
    }
    return $ret
}
