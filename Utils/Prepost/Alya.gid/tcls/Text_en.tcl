proc AlyaText { } {
    # bot ... boton labels
    # lab ... label texts
    # tit ... window titles
    global Text

    ######################
    # Alya information #
    ######################
    set Text(lab200) "About Alya"
    set Text(lab201) "Alya: A Package for the Numerical Simulation of Ventilation Problems"
    set Text(lab202) "Program written by:"
    set Text(lab203) "Please send any comment or suggestion to houzeaux\@cimne.upc.es"
    set Text(lab204) "Alya project definition"
    set Text(lab205) "Name"
    set Text(lab206) "Path"
    set Text(lab207) "Mesh characteristics"
    set Text(lab208) "Element type"
    set Text(lab209) "Number of nodes"
    set Text(lab210) "Number of elements"
    set Text(lab212) "Dimension"
    set Text(lab213) "Geometry not meshed"

    ##################
    # General botons #
    ##################
    set Text(bot001) "Ok"
    set Text(bot002) "Close"
    set Text(bot003) "Ok"
    set Text(bot004) "Yes"
    set Text(bot005) "No"
    set Text(bot006) "Finish"
    set Text(bot007) "Delete"
    set Text(bot008) "Select"
    set Text(bot009) "Draw"
    set Text(bot010) "Unselect"
    set Text(bot013) "Accept"
    set Text(bot014) "Select All"
    set Text(bot015) "Apply"
    set Text(bot016) "Pick"
    set Text(bot017) "Check"

    ####################
    # Pre-process menu #
    ####################
    set Text(lab001) "Help"
    set Text(lab002) "Wizard"
    set Text(lab003) "Log file"
    set Text(lab004) "Graphs"
    set Text(lab005) "Convergence"
    set Text(lab006) "Wall law"
    set Text(lab007) "Witness velocity"
    set Text(lab008) "Witness pressure"
    set Text(lab009) "Witness temperature"
    set Text(lab010) "Utilities"
    set Text(lab011) "Create a section"
    set Text(lab012) "Import from Fatherm"
    set Text(lab013) "Create a Landscape"
    set Text(lab014) "Compute parameters"
    set Text(lab015) "Flatten geometry"
    set Text(lab016) "About this project"
    set Text(lab017) "About Alya"
    set Text(lab018) "Change language"
    set Text(lab019) "Change Skin"
    set Text(lab0001) "Goto Fatherm"
    set Text(lab0002) "Force"
    set Text(lab0003) "Torque"

    #####################
    # Post-process menu #
    #####################
    set Text(lab020) "Help"
    set Text(lab021) "Calculate results"
    set Text(lab022) "PMV & PPD"
    set Text(lab023) "View results"
    set Text(lab024) "Select time step"
    set Text(lab025) "Velocity module"
    set Text(lab026) "Velocity vectors"
    set Text(lab027) "Pressure"
    set Text(lab028) "Temperature"
    set Text(lab029) "PMV"
    set Text(lab030) "PPD"
    set Text(lab031) "Wind chill equivalent temperature"
    set Text(lab032) "Wind chill index"
    set Text(lab033) "Wind chill equivalent temperature"
    set Text(lab034) "VELOC"
    set Text(lab035) "PRESS"
    set Text(lab036) "TEMPE"

    ##########
    # Errors #
    ##########
    set Text(lab100) "Result not available for selected time step"
    set Text(lab101) "Maybe you have to calculate it before..."
    set Text(lab102) "There is no Logfile for this project" 
    set Text(lab103) "There is no convergence for this project" 
    set Text(lab104) "There is no open project" 
    set Text(lab105) "There is no result file for this project"
    set Text(lab106) "There is no wall law file for this project"
    set Text(lab107) "There is no witness file for this project"
    set Text(lab108) "Wrong plane: the two points coincide"
    set Text(lab109) "An Error in the equation of the section plane has been found"
    set Text(lab110) "An error has been detected while intersecting surface number"
    set Text(lab111) "You must give a number"
    set Text(lab112) "Could not find Fatherm correspondance file. Results cannot be imported"
    set Text(lab113) "There is no Fatherm result file for this project" 
    set Text(lab114) "Invalid 2D section"
    set Text(lab115) "A possible error has been detected when importing the Fatherm results.\nWe continue anyway..." 
    set Text(lab116) "Velocity scale must be positive!"
    set Text(lab117) "Length scale must be positive!"
    set Text(lab118) "Temperature scale must be positive!"
    set Text(lab119) "Error when drawing! Restart GiD without saving..."
    set Text(lab120) "There is no geometry to create a section"
    set Text(lab121) "There is no geometry to create a landscape"
    set Text(lab122) "You must give a positive number"
    set Text(lab123) "Special characters and whitespaces are not allowed"
    set Text(lab124) "There is no associated Fatherm problem"
    set Text(lab125) "There are no zone to be meshed.\nYou must select the airflow zones before meshing..."

    ##################
    # Windows titles #
    ##################
    set Text(tit001) "Alya - Log file"
    set Text(tit002) "Alya - Graph"
    # Wizard for the CFD steps
    set Text(tit003) "Alya - Wizard"  
    # Select time step
    set Text(tit005) "Alya - Select Time"
    # Calculate PMV and PPD
    set Text(tit006) "Alya - PMV and PPD"
    # Calculate Section
    set Text(tit004) "Alya - Section"  
    set Text(tit007) "Alya - Section equation"
    set Text(tit008) "Alya - Project section"
    set Text(tit009) "Alya - Rotate section"
    # Postprocess
    set Text(tit010) "Alya - Velocity factor"
    # Landscape
    set Text(tit011) "Alya - Landscape"
    # Import results
    set Text(tit012) "Alya - Import results from Fatherm"
    # Parameters
    set Text(tit013) "Alya - Dimensionless parameters"
    # Zones
    set Text(tit014) "Alya - Select zones"
    # Insert Objects
    set Text(tit015) "Alya - Insert objects"
    # About Alya
    set Text(tit016) "Alya - About"
    # About this project
    set Text(tit017) "Alya - About this project"
    # Language
    set Text(tit018) "Alya - Change language"
    # Calculate Chill
    set Text(tit019) "Alya - Equivalent temperature"
    # About Alya
    set Text(tit020) "Alya - Error"
    # Layers
    set Text(tit1022) "Alya - Layers"
    # Wizard for the postprocess
    set Text(tit021) "Alya - Postprocess Wizard"  
    # Problem data
    set Text(tit022) "Alya - Problem data"  
    # Problem data
    set Text(tit023) "Alya - Generate mesh"  
    # Force and Torque
    set Text(tit024)  "Alya - Plot force"  
    set Text(tit1024) "Alya - Plot torque"  
    set Text(tit028)  "Alya - Create more surfaces" 
    # Turbulent Parameters
    set Text(tit025) "Alya - Turbulent parameters (based on Launder-Sharma)"

    ##########
    # Graphs #
    ##########
    set Text(lab300) "" 
    set Text(lab301) "Velocity" 
    set Text(lab302) "Pressure" 
    set Text(lab303) "Temperature" 
    set Text(lab304) "Number of iterations" 
    set Text(lab305) "Residual norm" 
    set Text(lab306) "y+ range"
    set Text(lab307) "Percentage"
    set Text(lab308) "Witness velocity"
    set Text(lab309) "Witness"
    set Text(lab310) "Time/Iteration"
    set Text(lab311) "Witness temperature"
    set Text(lab312) "Convergence history"
    set Text(lab313) "Wall law"
    set Text(lab314) "Witness history"
    set Text(lab315) "Evolution of force components"
    set Text(lab316) "Give the list of lines on which to sum the forces, by putting a space \nbetween them or selecting them with the Select button"
    set Text(lab317) "Give the list of surfaces on which to sum the forces, by putting a space \nbetween them or selecting them with the Select button"
    set Text(lab318) "List"
    set Text(lab319) "Scale factor"
    set Text(lab1315) "Evolution of torque components"
    set Text(lab1316) "Give the list of lines on which to sum the torques, by putting a space \nbetween them or selecting them with the Select button"
    set Text(lab1317) "Give the list of surfaces on which to sum the forces, by putting a space \nbetween them or selecting them with the Select button"
    set Text(lab1318) "List"
    set Text(lab1319) "Scale factor"
    set Text(lab1320) "Force"
    set Text(lab1321) "Torque"

    ##############
    # CFD Wizard #
    ##############
    set Text(lab150) "Data"
    set Text(lab151) "Edit the physical and numerical data"
    set Text(lab152) "Select Zones"
    set Text(lab153) "Select the connected zones to solve the airflow" 
    set Text(lab154) "Conditions"
    set Text(lab155) "Impose conditions on the domain boundary"
    set Text(lab156) "Meshing"
    set Text(lab157) "Generate the finite element mesh"
    set Text(lab158) "Save project" 
    set Text(lab159) "Save current project" 
    set Text(lab160) "Analysis"
    set Text(lab161) "Perform the simulation" 
    set Text(lab162) "Postprocess"
    set Text(lab163) "Visualize the results of the simulation" 
    set Text(lab164) "Create surfaces"
    set Text(lab165) "Create the surfaces automatically"
    set Text(lab166) "Geometry"
    set Text(lab167) "Create the computational domain" 
    set Text(lab168) "Draw geometry"
    set Text(lab169) "Import a dxf file"
    set Text(lab170) "Create a section"
    set Text(lab171) "Assign conditions"
    set Text(lab172) "Import from Fatherm"
    set Text(lab173) "Draw some graphs"
    set Text(lab174) "Generate mesh"

    ##################
    # Wizard Section #
    ##################
    set Text(lab248) "Rotate the geometry through a 90 degrees CW angle?"
    set Text(lab249) "Confirm you want to save this section? \nSelect No to cancel the section"
    set Text(lab250) "This wizard helps you to generate a section of a 3D geometry"
    set Text(lab251) "Define section plane"
    set Text(lab252) "Define the vertical or horizontal section"
    set Text(lab253) "Create section"
    set Text(lab254) "Create and visualize the section"
    set Text(lab255) "Save section"
    set Text(lab256) "Confirm and save the section in 2D"
    set Text(lab257) "Rotate section" 
    set Text(lab258) "Rotate the section in the XY plane. This is useful only for visualization purpose" 
    set Text(lab259) "Pick or choose two points on the XY plane. \nThe 2D section is defined by the corresponding \nline and the orthogonal direction to the screen."
    set Text(lab1259) "Enter the height of the horizontal section."
    set Text(lab246) "1st point"
    set Text(lab247) "2nd point"
    set Text(lab245) "Pick"
    set Text(lab244) "Point"
    set Text(lab1260) "Height relative to the ground floor (m)"

    ######################
    # Select a time step #
    ######################
    set Text(lab260) "Steps:" 
    set Text(lab261) "View:" 
    set Text(lab262) "factor:" 
    set Text(lab263) "Component:" 
    
    #############################
    # PMV and Chill Calculation #
    #############################
    set Text(lab265) "Enter the constants for the PMV/PPD calculation"
    set Text(lab266) "Clothing factor"
    set Text(lab267) "Metabolism" 
    set Text(lab268) "Wet metabolism" 
    set Text(lab269) "Radiative temperature"
    set Text(lab270) "Relative Humidity"
    set Text(lab271) "Ambient temperature" 
    set Text(lab272) "Enter the constant for the equivalent temperature (wind chill) calculation"

    ###############
    # Postprocess #
    ###############
    set Text(lab273) "Vector scale factor"
    set Text(lab274) "Result surface scale factor"

    #############
    # Landscape #
    #############
    set Text(lab275) "Edit the size of the landscape"

    ##################
    # Import results #
    ##################
    set Text(lab276) "Import boundary conditions from Fatherm"
    set Text(lab277) "Import"
    set Text(lab287) "Choose a result by pressing the corresponding arrow:"

    ############################
    # Dimensionless parameters #
    ############################
    set Text(lab278) "Compute some dimensionless parameters."
    set Text(lab279) "If U is set to zero, then it is calculated using the Boussinesq velocity scale" 
    set Text(lab280) "Velocity scale"
    set Text(lab281) "Length scale"
    set Text(lab282) "Temperature scale"
    set Text(lab283) "Reynolds number"
    set Text(lab284) "Raleigh number"
    set Text(lab285) "Grashof number"
    set Text(lab286) "Richardson number"

    ################
    # Select zones #
    ################
    set Text(lab290) "surfaces"
    set Text(lab291) "volumes"
    set Text(lab292) "Select the connected" 
    set Text(lab293) "where you want to solve the airflow"
    set Text(lab294) "Press 'Draw' to paint the selected zones"
    set Text(lab295) "Press 'Finish' \nwhen finished"
    set Text(lab296) "Press 'Finish' \nwhen finished"
    set Text(lab297) "Press 'Finish' \nwhen finished"

    ###############
    # Put objects #
    ###############
    set Text(lab320) "Enter the final position of the object"

    #############
    # Languages #
    #############
    set Text(lab330) "Language"
    set Text(lab331) "Select a language in the following list"

    #############
    # Start GiD #
    #############
    set Text(lab332) "This is Alya Version"
    set Text(lab333) ""
    set Text(lab334) "Version"
    set Text(lab335) "in English" 
    set Text(lab336) "Press Ctrl-w to show/hide Alya Windows." 

    #######################
    # Import from Fatherm #
    #######################
    set Text(lab340) "Boundary Conditions imported from Fatherm:"
    set Text(lab341) "Walls:   "
    set Text(lab342) "Inflows: "
    set Text(lab343) "No boundary condition imported from Fatherm"
    set Text(lab344) "Outflows: "
    set Text(lab345) "Opened:   "
    ############################
    # Dimensionless parameters #
    ############################
    set Text(lab378) "Compute some turbulent parameters."
    set Text(lab379) " " 
    set Text(lab380) "Turbulent kinetic energy"
    set Text(lab381) "Turbulent dissipation"
    set Text(lab382) "Turbulent dis. specific rate"
    set Text(lab383) "Turbulent viscosity"
    set Text(lab384) "Turbulent Reynolds number"
    set Text(lab385) "Viscosity ratio"

    ##########
    # Layers #
    ##########
    set Text(lab1096) ""
    set Text(lab1097) "Colour"
    set Text(lab1098) "On"
    set Text(lab1099) "Off"
    set Text(lab1100) "Changes the color used to draw the entities in render mode."
    set Text(lab1101) "Entities belonging to a layer ON will be drawn."
    set Text(lab1102) "Entities belonging to a layer OFF will not be drawn."
    set Text(lab1103) "Entities belonging to a freeze layer will be drawn but cannot be selected or modified." 
    set Text(lab1104) "Entities belonging to a unfreeze layer will be selected and modified."
    set Text(lab1105) "Select layers with filters (for example: 'b*' select all the layers with name beginning by 'b')"
    set Text(lab1106) "Sel"
    set Text(lab1107) "Layer To use "
    set Text(lab1108) "Set the current layer (the new entities are created in this layer)."
    set Text(lab1109) "This option permmits to select one entity and retrieve its layer."
    set Text(lab1110) "New"
    set Text(lab1111) "Create a new layer."
    set Text(lab1112) "Delete"
    set Text(lab1113) "Delete a layer."
    set Text(lab1114) "Sort the list of layers by name."
    set Text(lab1115) "To back"
    set Text(lab1116) "Send entities to back/front of a layer (entities in back are not drawn and frozen)"
    set Text(lab1117) "Also lower entities"
    set Text(lab1118) "Opposite"
    set Text(lab1119) "Send To"
    set Text(lab1120) "Send the selected entities to a layer. "
    set Text(lab1121) "If 'Also lower entities' is checked, the subentities are also send to the layer"
    set Text(lab1122) "before applying operation delete it is necessary to select layers"
    set Text(lab1123) "Layer '"
    set Text(lab1124) "' can't be erased.\n(see the message bar for details)"
    set Text(lab1125) "OK to all"
    set Text(lab1126) "before applying operation colour it is necessary to select layers"
    set Text(lab1127) "Select color"
    set Text(lab1128) "before applying operation "
    set Text(lab1129) " it is necessary to select layers"
    set Text(lab1130) "before applying 'Entities' it is necessary to select 1 layer"
    set Text(lab1131) "Press 'Finish' when selection\nis finished"
    set Text(lab1132) "Points Lines Surfaces Volumes Dimensions All"
    set Text(lab1133) "Nodes Elements Dimensions All" 
    set Text(lab1134) "Selected entity does not belong to one layer"
    set Text(lab1135) "Point Line Surface Volume"
    set Text(lab1136) "Node Element"
    set Text(lab1137) "Points Lines Surfaces Volumes All"
    set Text(lab1138) "Nodes Elements All"
    set Text(lab1139) "Bring to front"
    set Text(lab1140) "Bring ALL to front"
    set Text(lab1141) "A layer name must be written before selecting 'New'"
    set Text(lab1142) "Layer'"
    set Text(lab1143) "' already exists"
    set Text(lab1144) "Create a new storey"
    set Text(lab1145) "Storey number:"
    set Text(lab1146) "Storey number must be a numerical value"
    set Text(lab1147) "Create a new layer"
    set Text(lab1148) "Create a new component in "
    set Text(lab1149) "Insert a new component"
    set Text(lab1150) "Choose a component to insert in"
    set Text(lab1151) "Create a new component"
    set Text(lab1152) "Component name:" 
    set Text(lab1153) "On/Off:"
    set Text(lab1154) "Freeze/Unfreeze:"
    set Text(lab1155) "Change color:"
    set Text(lab1156) "New/Delete:"
    set Text(lab1157) "Send to:"
    set Text(lab1158) "To back:"
    set Text(lab1159) "New"
    set Text(lab1160) "Del"
    set Text(lab1161) "Create a new layer"
    set Text(lab1162) "Name layer:"
    set Text(lab1163) "Create a new layer in "
    set Text(lab1164) "Insert a new layer"
    set Text(lab1165) "Ren"
    set Text(lab1166) "Rename a layer."
    set Text(lab1167) "Change a name layer"

    ##########################
    # CFD Postprocess Wizard #
    ##########################
    set Text(lab2150) "Graphs"
    set Text(lab2151) "Plot the selected graph"
    set Text(lab2152) "Show result"
    set Text(lab2153) "Show the selected variable"
    set Text(lab2154) "Velocity vectors"
    set Text(lab2155) "Plot the velocity vectors" 
    set Text(lab2156) "Preprocess"
    set Text(lab2157) "Go back to the preprocess" 

    ################
    # Problem data #
    ################
    set Text(lab2180) "Do you want to use default data values?\nPress No if you want to edit your own data"

    ###################
    # Mesh generation #
    ###################
    set Text(lab2181) "Do you want Alya to calculate the mesh size automatically?"
    set Text(lab2182) "Choose the density of the mesh"
    set Text(lab2183) "High"
    set Text(lab2184) "Low"

    #####################
    # Friendly messages #
    #####################
    set Text(lab2200) "All zones have been selected"
    set Text(lab2201) "You have selected STEDI-Vent default problem data"

    ### Create new surfaces
    set Text(lab352)  "If Alya hasn't automaticaly created all the surfaces, press 'Pick' and select an interior point to the surface you want to create"

} 

