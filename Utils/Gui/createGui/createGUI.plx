#!/usr/bin/perl

#usage modules
use FindBin;                 # locate this script
use lib "$FindBin::Bin/libraries";  # use the parent directory
use File::Copy;
use File::Path;
use XML::Simple;
use Data::Dumper;
use POSIX;

#global variables
$basePath = "../AlyaGuiSource";

#----------------------------------------------------------------------
#-------------------------------functions------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#--Function createUIOuPos
#--Description:
#----Genereates the ui file part corresponding to an "output and postprocess group"
#--Parameters:
#----group: The group metadata
#--Returns: the ui file part
sub createUIOuPos {
	local($group) = @_;

	#open the ui template file
	open FILE, "<resources/outputGroup.ui";
	my $groupUI = do { local $/; <FILE> };
	close FILE;
	
	#postProcess vars
	#get the postProcess vars and replaces it in the template
	$vars = $group->{postProcess}->{var};
	@varsList = $vars;
	if (ref $vars eq 'ARRAY') {
		@varsList = @{$vars};
	}
	my $varsData;
	foreach $var (@varsList) {
		open FILE, "<resources/textItem.ui";
		my $textItem = do { local $/; <FILE> };
		close FILE;
		$textItem =~ s/<<text>>/$var/;
		$varsData = $varsData . $textItem
	}
	$groupUI =~ s/<<postProcessVars>>/$varsData/;

	#elementSet vars
	#get the elementSet vars and replaces it in the template
	$vars = $group->{elementSet}->{var};
	@varsList = $vars;
	if (ref $vars eq 'ARRAY') {
		@varsList = @{$vars};
	}
	$varsData = "";
	foreach $var (@varsList) {
		open FILE, "<resources/textItem.ui";
		my $textItem = do { local $/; <FILE> };
		close FILE;
		$textItem =~ s/<<text>>/$var/;
		$varsData = $varsData . $textItem
	}
	$groupUI =~ s/<<elementSetVars>>/$varsData/;

	#boundarySet vars
	#get the boundarySet vars and replaces it in the template
	$vars = $group->{boundarySet}->{var};
	@varsList = $vars;
	if (ref $vars eq 'ARRAY') {
		@varsList = @{$vars};
	}
	$varsData = "";
	foreach $var (@varsList) {
		open FILE, "<resources/textItem.ui";
		my $textItem = do { local $/; <FILE> };
		close FILE;
		$textItem =~ s/<<text>>/$var/;
		$varsData = $varsData . $textItem
	}
	$groupUI =~ s/<<boundarySetVars>>/$varsData/;

	#nodeSet vars
	#get the nodeSet vars and replaces it in the template
	$vars = $group->{nodeSet}->{var};
	@varsList = $vars;
	if (ref $vars eq 'ARRAY') {
		@varsList = @{$vars};
	}
	$varsData = "";
	foreach $var (@varsList) {
		open FILE, "<resources/textItem.ui";
		my $textItem = do { local $/; <FILE> };
		close FILE;
		$textItem =~ s/<<text>>/$var/;
		$varsData = $varsData . $textItem
	}
	$groupUI =~ s/<<nodeSetVars>>/$varsData/;

	#witnessPoints vars
	#get the witnessPoints vars and replaces it in the template
	$vars = $group->{witnessPoints}->{var};
	@varsList = $vars;
	if (ref $vars eq 'ARRAY') {
		@varsList = @{$vars};
	}
	$varsData = "";
	foreach $var (@varsList) {
		open FILE, "<resources/textItem.ui";
		my $textItem = do { local $/; <FILE> };
		close FILE;
		$textItem =~ s/<<text>>/$var/;
		$varsData = $varsData . $textItem
	}
	$groupUI =~ s/<<witnessPointsVars>>/$varsData/;

	return $groupUI;	
}

#----------------------------------------------------------------------
#--Function createUIBouCo
#--Description:
#----Genereates the ui file part corresponding to an "boundary conditions group"
#--Parameters:
#----group: The group metadata
#--Returns: the ui file part
sub createUIBouCo {
	local($group) = @_;

	#open the ui template file
	open FILE, "<resources/boucoGroup.ui";
	my $groupUI = do { local $/; <FILE> };
	close FILE;
	
	#postProcess vars
	#get the postProcess vars and replaces it in the template
=pod
	$vars = $group->{postProcess}->{var};
	@varsList = $vars;
	if (ref $vars eq 'ARRAY') {
		@varsList = @{$vars};
	}
	my $varsData;
	foreach $var (@varsList) {
		open FILE, "<resources/textItem.ui";
		my $textItem = do { local $/; <FILE> };
		close FILE;
		$textItem =~ s/<<text>>/$var/;
		$varsData = $varsData . $textItem
	}
	$groupUI =~ s/<<postProcessVars>>/$varsData/;
=cut
	return $groupUI;	
}

sub getInputLines {
	local($subGroup) = @_;
	my $xml = new XML::Simple;
	my $data = $xml->XMLin("alyaGui.xml");
	my $subGroupType = $subGroup->{subGroupType};
	my @inputLinesList;

	if ($subGroupType && ($subGroupType ne "service") && ($subGroupType ne "module")) {
		#obtain the inputLines from a subGroupType.
		#find the subgroup
		$subGroupTypes = $data->{subGroupTypes}->{subGroupType};
		@subGroupTypesList = $subGroupTypes;
		if (ref $subGroupTypes eq 'ARRAY') {
			@subGroupTypesList = @{$subGroupTypes};
		}
		foreach $subGroupTypeClass (@subGroupTypesList) {
			if ($subGroupTypeClass->{subGroupTypeName} eq $subGroupType) {
				#----inputLine
				$inputLines = $subGroupTypeClass->{inputLine};
				@inputLinesList = $inputLines;
				if (ref $inputLines eq 'ARRAY') {
					@inputLinesList = @{$inputLines};
				}			
			}
		}
	}
	else {
		#----inputLine
		$inputLines = $subGroup->{inputLine};
		@inputLinesList = $inputLines;
		if (ref $inputLines eq 'ARRAY') {
			@inputLinesList = @{$inputLines};
		}	
	}
	return @inputLinesList;
}

#----------------------------------------------------------------------
#--Function createUIGroup
#--Description:
#----Genereates the ui file corresponding to an "standard group"
#--Parameters:
#----groupName: The group name
#----group: The group metadata
#----gidx: The group index(the group's list is ordered, this is the index)
#--Returns: the ui file part
sub createUIGroup {	
	local($groupName, $group, $gidx) = @_;
	print "group name:" . $groupName . "\n";

	#open the group template ui file
	open FILE, "<resources/group.ui";
	my $groupUI = do { local $/; <FILE> };
	close FILE;

	#replaces the main tokens on the template by the values from the xml
	$groupUI =~ s/#tabName/$groupName/;
	$groupUI =~ s/#tabLabel/$groupName/;
	my $subs = $groupName . "Layout";
	$groupUI =~ s/#layoutName/$subs/;
	$subs = $groupName . "Scroll";
	$groupUI =~ s/#scrollAreaName/$subs/;
	$subs = $groupName . "ScrollContent";
	$groupUI =~ s/#scrollAreaContentName/$subs/;
	$subs = $groupName . "VertLayout";
	$groupUI =~ s/#verticalLayoutName/$subs/;
	
	#generates the subgroup ui data file
	my $subGroupData = "";
	#----iterate over each subgroup
	$subGroups = $group->{subGroup};
	@subGroupsList = $subGroups;
	if (ref $subGroups eq 'ARRAY') {
		@subGroupsList = @{$subGroups};
	}	
	#-------------------------------------
	my $sgidx = 0;
	foreach $subGroup (@subGroupsList) {
		$subGroupName = $subGroup->{subGroupName};
		$subGroupValue = $subGroup->{subGroupValue};
		$subGroupCheckeable = $subGroup->{subGroupCheckeable};

		#creates an unique subgroup identifier
		$subGroupIdx = "sg" . $subGroupName . $gidx . $sgidx;
		#opens the subroup template
		open FILE, "<resources/subGroup.ui";
		my $subGroupUI = do { local $/; <FILE> };
		close FILE;
		
		#replaces the subgroup template tokens by its corresponding values
		$subGroupUI =~ s/#itemName/$subGroupIdx/;
		$auxName = $subGroupName;
		if ($subGroupValue) {
			$auxName = $auxName . ": " . $subGroupValue;
		}
		$subGroupUI =~ s/#itemTitle/$auxName/;
		$subs = $subGroupIdx . "Layout";
		$subGroupUI =~ s/#itemLayoutName/$subs/;
		#Determines if its checkeable or not
		$checkable = "false";
		if ($subGroupCheckeable eq "true") {
			$checkable = "true";
		}
		$subGroupUI =~ s/<<checkable>>/$checkable/;
		@inputLinesList = getInputLines($subGroup);
		
		#creates each input line ui data
		my $inputLineData = "";
		my $ilidx = 0;
		foreach $inputLine (@inputLinesList) {
			if ($inputLine) {
				$inputLineName = $inputLine->{inputLineName};
				$inputLineIdx = $inputLineName . $gidx . $sgidx . $ilidx;
				#opens the ui inputLine template file
				open FILE, "<resources/inputLine.ui";
				my $inputLineUI = do { local $/; <FILE> };
				close FILE;
				#replaces the template tokens by the corresponding values
				$subs = $inputLineIdx . "Layout";
				$inputLineUI =~ s/#lineLayoutName/$subs/;
				$subs = $inputLineIdx . "CheckBox";
				$inputLineUI =~ s/#inputLineName/$subs/;
				$inputLineUI =~ s/#inputLineText/$inputLineName/;
				$subs = $inputLineIdx . "Frame";
				$inputLineUI =~ s/#frameName/$subs/;
				$subs = $inputLineIdx . "FrameLayout";
				$inputLineUI =~ s/#frameLayoutName/$subs/;
				#----inputElement
				my $inputElementData = "";
				$inputElements = $inputLine->{inputElement};
				@inputElementsList = $inputElements;
				if (ref $inputElements eq 'ARRAY') {
					@inputElementsList = @{$inputElements};
				}	
				#creates the ui inputElement data
				my $inputElementsData = "";
				my $ieidx = 0;
				my $editIndex = 1;
				my $withList = 0;
				foreach $inputElement (@inputElementsList) {
					$type = $inputElement->{inputElementType};
					$inputName = $inputElement->{inputLineEditName};
					$inputElementName = $inputElement->{inputElementName};
          my $inputElementValueType = $inputElement->{inputElementValueType};
					$inputElementIdx = "ie" . $inputName . $gidx . $sgidx . $ilidx . $ieidx;
					#creates the ui for the edit input type: a line edit with name
					if ($type eq "edit") {
						open FILE, "<resources/edit.ui";
						my $editUI = do { local $/; <FILE> };
						close FILE;
						$subs = $inputElementIdx . "Label";
						$editUI =~ s/#labelName/$subs/;
						$editUI =~ s/#labelText/$inputName/;
            $editUI =~ s/#tooltip/$inputElementValueType/;
						$subs = $inputElementIdx . "Edit";
						$editUI =~ s/#editName/$subs/;
						$inputElementsData = $inputElementsData . $editUI;
					}
					#creates the ui for the edite input type: a line edit without name, only inser commas
					elsif ($type eq "edit2") {
						open FILE, "<resources/edit.ui";
						my $editUI = do { local $/; <FILE> };
						close FILE;
						$subs = $inputElementIdx . "Label";
						$editUI =~ s/#labelName/$subs/;
						if ($editIndex != 1) {
							$editUI =~ s/#labelText/,/;
						}
						else {
							$editUI =~ s/#labelText//;
						}
            $editUI =~ s/#tooltip/$inputElementValueType/;					
						$subs = $inputElementIdx . "Edit";
						$editUI =~ s/#editName/$subs/;
						$inputElementsData = $inputElementsData . $editUI;
						$editIndex = $editIndex + 1;
					}
					#creates the ui for the combo input type
					elsif ($type eq "combo") {							
						open FILE, "<resources/combo.ui";
						my $comboUI = do { local $/; <FILE> };
						close FILE;
						$subs = $inputElementIdx . "Combo";
						$comboUI =~ s/#comboName/$subs/;
						#check label						
						my $label;
						if ($inputElementName) {
							open FILE, "<resources/label.ui";
							my $labelUI = do { local $/; <FILE> };
							close FILE;
							$subs = $inputElementIdx . "Label";
							$labelUI =~ s/#labelName/$subs/;
							$labelUI =~ s/#labelText/$inputElementName/;
							$label = $labelUI;
						}
						$comboUI =~ s/<<label>>/$label/;

						my $comboOptions = "";					
						$items = $inputElement->{item};
						@itemsList = $items;
						if (ref $items eq 'ARRAY') {
							@itemsList = @{$items};
						}	
						#creates the combo values
						foreach $item (@itemsList) {
							my $itemName = $item->{itemName};
							$option = "<item>\n" . 
								  "<property name=\"text\">\n" .
								  "<string>$itemName</string>\n" .
								  "</property>\n" .
								  "</item>\n";
							$comboOptions = $comboOptions . $option;
						}
						$comboUI =~ s/#item/$comboOptions/;
						$inputElementsData = $inputElementsData . $comboUI;
					}
					elsif ($type eq "list") {
						$withList = 1;
						open FILE, "<resources/list.ui";
						my $listUI = do { local $/; <FILE> };
						close FILE;
						my $subs = $inputElementIdx . "List";
						$listUI =~ s/#listName/$subs/g;					
						my $listOptions = "";					
						$items = $inputElement->{item};
						@itemsList = $items;
						if (ref $items eq 'ARRAY') {
							@itemsList = @{$items};
						}	
						#creates the combo values
						foreach $item (@itemsList) {
							my $itemName = $item->{itemName};
							$option = "<column>\n".
	                                   "<property name=\"text\">\n".
		                               "<string>$itemName</string>\n".
	                                   "</property>\n".
                                       "</column>";
							$listOptions = $listOptions . $option;
						}
						$listUI =~ s/<<columns>>/$listOptions/;
						$inputElementsData = $inputElementsData . $listUI;
					}
					$ieidx = $ieidx + 1;
				}
				#creates the helb button ui, if the inputline has help text
				if ($inputLine->{inputLineHelp}) {
					#add the help button
					open FILE, "<resources/helpButton.ui";
					my $helpUI = do { local $/; <FILE> };
					close FILE;
					$subs = $inputLineIdx . "HelpButton";
					$helpUI =~ s/#helpButtonName/$subs/;
					$inputLineUI =~ s/#helpButton/$helpUI/;
				}
				else {
					$inputLineUI =~ s/#helpButton//;
				}
				#special characteristics for lists:
				my $inputLineHeight = "50";
				my $alignleft = " alignment=\"Qt::AlignLeft\"";
				if ($withList) {
					$inputLineHeight = "500";
					$alignleft = "";					
				}
				$inputLineUI =~ s/#alignleft/$alignleft/;
				$inputLineUI =~ s/#height/$inputLineHeight/;
				
				#add to inputLine the inputElements
				$inputLineUI =~ s/#inputElement/$inputElementsData/;
				$inputLineData = $inputLineData . $inputLineUI;
				$ilidx = $ilidx + 1;
			}						
		}
		#Add to subgroup the inputLine
		$subGroupUI =~ s/#items/$inputLineData/;
		$subGroupData = $subGroupData . $subGroupUI; 
		$sgidx = $sgidx + 1;
	}
	#Add to group the subGroupData info
	$groupUI =~ s/#subgroup/$subGroupData/;
	return $groupUI;
}


#----------------------------------------------------------------------
#--Function createUIFile
#--Description:
#----Genereates the ui file corresponding to a module
#--Parameters:
#----moduleName: The module name
#----module: The module metadata
#--Returns: the ui file part corresponding to a module
sub createUIFile {

local($moduleName, $module) = @_;
	#-------------------------------------------------------------------------------
        #---------------------------create the UI file----------------------------------
	#-------------------------------------------------------------------------------
	#Open the module ui template file
	open FILE, "<resources/module.ui";
	my $moduleUI = do { local $/; <FILE> };
	close FILE;

	#replaces the template tokens by its corresponding values
	my $classModuleName = ucfirst($moduleName);
	$moduleUI =~ s/#class/$classModuleName/;
	$moduleUI =~ s/#name/$classModuleName/;
	#---group
	$groups = $module->{group};
	@groupsList = $groups;
	if (ref $groups eq 'ARRAY') {
		@groupsList = @{$groups};
	}	
	#creates the ui group file part
	my $groupData = "";
	my $gidx = 0;
	foreach $group (@groupsList) {
		$groupType = $group->{groupType};
		$groupName = $group->{groupName};
		#group type oupos, generates an output and postprocess group ui
		if ($groupType eq "oupos") {
			$groupData = $groupData . createUIOuPos($group);
		}
		elsif ($groupType eq "bouco") {
			$groupData = $groupData . createUIBouCo($group);
		}
		#other group type, generates standard group ui
		else {
			$groupData = $groupData . createUIGroup($groupName, $group, $gidx);
		}
		$gidx = $gidx + 1;
	}
	#Add group data to module
	$moduleUI =~ s/#group/$groupData/;

	return $moduleUI;

}

#----------------------------------------------------------------------
#--Function createSavePart
#--Description:
#----Generates the save part code for a default group, this part is in charge of generating the alya files from the ui values
#--Parameters:
#----groupName: The group name
#----gidx: The group index(the group's list is ordered, this is the index)
#--Returns: the save file part
sub createSavePart {
	local($groupName, $group, $gidx) = @_;

	#savesData contains the c++ code to generate the alya files from the ui group values
  my $savesData = "file << \"$groupName\" << std::endl;\n"; 
	
	my $subGroupData = "";
	#----subgroup
	$subGroups = $group->{subGroup};
	@subGroupsList = $subGroups;
	if (ref $subGroups eq 'ARRAY') {
		@subGroupsList = @{$subGroups};
	}	
	#iterates over each subGroup and generates his save part
	$sgidx = 0;
	foreach $subGroup (@subGroupsList) {
		$subGroupName = $subGroup->{subGroupName};
		$subGroupType = $subGroup->{subGroupType};
		$subGroupValue = $subGroup->{subGroupValue};
		$subGroupCheckeable = $subGroup->{subGroupCheckeable};
		$subGroupIdx = "sg" . $subGroupName . $gidx . $sgidx;
		
		#if the subgroup is checkeable then creates the "if" statement to check if it is checked
		if ($subGroupCheckeable eq "true") {
			$savesData = $savesData . "if (ui->$subGroupIdx->isChecked()) {\n";
		}
		if ($subGroupName) {
			#is subgroup type is a module or service then it generates the "if statement" to check if it is enabled
			#also generates a ": On" on the alya file
			if ($subGroupType eq "service" || $subGroupType eq "module") {				
				$savesData = $savesData . "if (ui->$subGroupIdx->isEnabled()) {\n";
				$savesData = $savesData . "    file << \"$subGroupName: On\" << std::endl;\n";
			}
			#for other subgroup types it generates a default save code
			else {
				if ($subGroupValue) {
					$savesData = $savesData . "file << \"$subGroupName: $subGroupValue\" << std::endl;\n";
				}
				else {
					$savesData = $savesData . "file << \"$subGroupName\" << std::endl;\n";
				}				
			}			
		}
		#check if we have to print the subGroupType value
		if ($subGroupType && ($subGroupType ne "service") && ($subGroupType ne "module")) {
			$savesData = $savesData . "file << \"$subGroupType\" << std::endl;\n";
		}			
		#----inputLine
		@inputLinesList = getInputLines($subGroup);	
		#creates a save statement for each subGroup inputLine
		my $inputLineData = "";
		my $ilidx = 0;
		foreach $inputLine (@inputLinesList) {
			if ($inputLine) {
				$inputLineName = $inputLine->{inputLineName};
				$inputLineIdx = $inputLineName . $gidx . $sgidx . $ilidx;
			
				open FILE, "<resources/checkBoxSave.cpp";
				my $template = do { local $/; <FILE> };
				close FILE;
				$template =~ s/<<checkName>>/$inputLineName/g;
				$template =~ s/<<checkIdx>>/$inputLineIdx/g;
				if (!$inputLine->{inputElement}) {
					$template =~ s/<<checkValue>>/ON/g;	
				}
				else {
					$template =~ s/<<checkValue>>//g;
				}

				#----inputElement
				my $inputElementData = "";
				$inputElements = $inputLine->{inputElement};
				@inputElementsList = $inputElements;
				if (ref $inputElements eq 'ARRAY') {
					@inputElementsList = @{$inputElements};
				}
				#creates a save statement for each inputLine inputElement
				$inputIndex = 0;
				$ieidx = 0;
				foreach $inputElement (@inputElementsList) {
					$type = $inputElement->{inputElementType};
					$inputName = $inputElement->{inputLineEditName};
					$inputElementName = $inputElement->{inputElementName};
					$inputElementIdx = "ie" . $inputName . $gidx . $sgidx . $ilidx . $ieidx;
					#creates save statement for the edit type
					if ($type eq "edit") {
						open FILE, "<resources/edit.cpp";
						my $edit = do { local $/; <FILE> };
						close FILE;
						$edit =~ s/<<editIdx>>/$inputElementIdx/g;
						$edit =~ s/<<editName>>/$inputName/g;
						if ($inputIndex == 0) {
							$edit =~ s/<<comma>>//g;
						}
						else {
							$edit =~ s/<<comma>>/,/g;
						}					            
						$inputElementData = $inputElementData . $edit . "\n";
					}
					#creates save statement for the edit2 type(without name)
					elsif ($type eq "edit2") {						
						open FILE, "<resources/edit2.cpp";
						my $edit = do { local $/; <FILE> };
						close FILE;
						$edit =~ s/<<editIdx>>/$inputElementIdx/g;
						if ($inputIndex == 0) {
							$edit =~ s/<<comma>>//g;
						}
						else {
							$edit =~ s/<<comma>>/,/g;
						}					            
						$inputElementData = $inputElementData . $edit . "\n";
					}
					#creates save statement for the combo type
					elsif ($type eq "combo") {
						open FILE, "<resources/combo.cpp";
						my $combo = do { local $/; <FILE> };
						close FILE;
						$combo =~ s/<<comboIdx>>/$inputElementIdx/g;
						my $comboName;
						if ($inputElementName) {
							$comboName = $inputElementName . ": ";
						}
						if ($inputIndex == 0) {
							$combo =~ s/<<comma>>/ $comboName/g;
						}
						else {
							$combo =~ s/<<comma>>/, $comboName/g;
						}
						$inputElementData = $inputElementData . $combo . "\n";
					}
					#creates save statement for the list type
					elsif ($type eq "list") {
						open FILE, "<resources/listSave.cpp";
						my $list = do { local $/; <FILE> };
						close FILE;
						my $listName = $inputElementIdx . "List";
						$list =~ s/#listName/$listName/g;
						$list =~ s/#elementName/$inputLineName/g;
						
						$items = $inputElement->{item};
						@itemsList = $items;
						if (ref $items eq 'ARRAY') {
							@itemsList = @{$items};
						}	
						my $rowNum = 0;
						my $listRows;
						foreach $item (@itemsList) {
							my $itemName = $item->{itemName};
							open FILE, "<resources/listRowSave.cpp";
							my $listRow = do { local $/; <FILE> };
							close FILE;
							$listRow =~ s/#listName/$listName/g;
							$listRow =~ s/#rownum/$rowNum/g;
							$listRow =~ s/#rowname/$itemName/g;
							if ($rowNum == 0) {
								$listRow =~ s/#comma//g;
						    }
						    else {
								$listRow =~ s/#comma/,/g;
							}						
							$rowNum = $rowNum + 1;
							$listRows = $listRows . $listRow;
						}
						$list =~ s/#rowSaves/$listRows/;
						
						$inputElementData = $inputElementData . $list . "\n";
					}
					$inputIndex = $inputIndex + 1;
					$ieidx = $ieidx + 1;
				}
				$template =~ s/<<inputLine>>/$inputElementData/g;
				$savesData = $savesData . $template;
				$ilidx = $ilidx + 1;
			}
		}
		#check if we have to print the subGroupType value
		if ($subGroupType && ($subGroupType ne "service") && ($subGroupType ne "module")) {
			$savesData = $savesData . "file << \"END_$subGroupType\" << std::endl;\n";
		}
		if ($subGroupName) {
			#end subgroup save statement
			$savesData = $savesData . "    file << \"END_$subGroupName\" << std::endl;\n";
			if ($subGroupType eq "service" || $subGroupType eq "module") {
				$savesData = $savesData . "}";
			}
		}
		if ($subGroupCheckeable eq "true") {
			$savesData = $savesData . "}\n";
		}
		$sgidx = $sgidx + 1;
	}
	#end group save statement
	$savesData = $savesData . "file << \"END_$groupName\" << std::endl;\n";
	return $savesData;
}

#----------------------------------------------------------------------
#--Function createOuPosSave
#--Description:
#----Generates the save part code for a oupos group, output and postprocess special group,
#----this part is in charge of generating the alya files from the ui values
#--Parameters:
#----group: The group metadata
#--Returns: the save file part
sub createOuPosSave {
	local($group) = @_;

	open FILE, "<resources/ouposSave.cpp";
	my $ouposSave = do { local $/; <FILE> };
	close FILE;

	return $ouposSave;
	
}

#----------------------------------------------------------------------
#--Function createBouCoSave
#--Description:
#----Generates the save part code for a bouco group, boundary conditions special group,
#----this part is in charge of generating the alya files from the ui values
#--Parameters:
#----group: The group metadata
#--Returns: the save file part
sub createBouCoSave {
	local($group) = @_;

	open FILE, "<resources/boucoSave.cpp";
	my $ouposSave = do { local $/; <FILE> };
	close FILE;

	return $ouposSave;
	
}

#----------------CONDITIONAL COMBOS------------------------
sub conditionalCombos {
	local($inputLine, $gidx, $sgidx, $ilidx, $moduleClass, $part) = @_;

	my $inputElementData = "";
	$inputElements = $inputLine->{inputElement};
	@inputElementsList = $inputElements;
	if (ref $inputElements eq 'ARRAY') {
		@inputElementsList = @{$inputElements};
	}	
	#creates the ui inputElement data
	$inputElementsData = "";
	$ieidx = 0;
	$editIndex = 1;

	#montage estructuras
        my $depCombos = {};
	my $optionalElements = {};
	foreach $inputElement (@inputElementsList) {
		$type = $inputElement->{inputElementType};
		$inputName = $inputElement->{inputLineEditName};
		$inputElementName = $inputElement->{inputElementName};
		$inputElementIdx = "ie" . $inputName . $gidx . $sgidx . $ilidx . $ieidx;

		#-------Rellenado lista de combos con dependencias
		if ($type eq "combo") {
			$items = $inputElement->{item};
			@itemsList = $items;
			if (ref $items eq 'ARRAY') {
				@itemsList = @{$items};
			}	
			#creates the combo values
			my $depItems = {};
			foreach $item (@itemsList) {
				my $itemName = $item->{itemName};
				if (exists $item->{itemDependence}) {
					my $itemDep = $item->{itemDependence};
					$depItems->{$itemName} = $itemDep;
				}
			}			
			my $numKeys = 0;
			$numKeys += keys %$depItems;
			if ($numKeys > 0) {
				my $comboName = $inputElementIdx. "Combo"; 
				$depCombos->{$comboName}=$depItems;
			}
		}
		#-------

		#-------Rellenado lista elementos opcionales
		if (exists $inputElement->{inputElementGroup}) {
			my $eleGroup =  $inputElement->{inputElementGroup};
			
			if ($type eq "edit") {
				my $labelName = $inputElementIdx . "Label";
				my $editName = $inputElementIdx . "Edit";
				push( @{$optionalElements->{$eleGroup}}, $labelName);
				push( @{$optionalElements->{$eleGroup}}, $editName);
			}
			#creates the ui for the edite input type: a line edit without name, only inser commas
			elsif ($type eq "edit2") {
				my $labelName = $inputElementIdx . "Label";
				my $editName = $inputElementIdx . "Edit";
				push( @{$optionalElements->{$eleGroup}}, $labelName);
				push( @{$optionalElements->{$eleGroup}}, $editName);
			}
			#creates the ui for the combo input type
			elsif ($type eq "combo") {
				my $editName = $inputElementIdx . "Combo";
				push( @{$optionalElements->{$eleGroup}}, $editName);				
				#check label
				if ($inputElementName) {
					my $labelName = $inputElementIdx . "Label";
					push( @{$optionalElements->{$eleGroup}}, $labelName);
				}				
			}
			#creates the ui for the list input type
			elsif ($type eq "list") {
				my $editName = $inputElementIdx . "List";
				push( @{$optionalElements->{$eleGroup}}, $editName);
				#add the list plus and minus buttons
				my $listNamePlus = $inputElementIdx . "ListPlus";
				push( @{$optionalElements->{$eleGroup}}, $listNamePlus);
				my $listNameMinus = $inputElementIdx . "ListMinus";
				push( @{$optionalElements->{$eleGroup}}, $listNameMinus);								
			}
		}
		#-------------------------------------------
		$ieidx = $ieidx + 1;
	}

	#RECORRIDO ESTRUCTURAS
	print "Structure:" . Dumper($depCombos);
	print "OPTELEMENTS:" . Dumper($optionalElements);
	#Contiene la implementacion de las slots
	my $slotFunctions;
	#Contiene la definición de las slots
	my $slotFunctionsDef;
	#Contiene la llamada al slot cuando se crea el modulo
	my $slotCall;
	for $comboKey ( keys %$depCombos ) {
		$moduleClassUp = ucfirst($moduleClass);
		my $slotFunction = "void $moduleClassUp" . "::on_" . $comboKey . "_currentIndexChanged(const QString &arg1) {\n";
		$slotFunctionsDef = $slotFunctionsDef . "void on_" . $comboKey . "_currentIndexChanged(const QString &arg1);\n";
		$slotCall = $slotCall . "on_".$comboKey."_currentIndexChanged(ui->".$comboKey."->currentText());";
		for $itemKey ( keys %{$depCombos->{$comboKey}} ) {
			my $value = $depCombos->{$comboKey}->{$itemKey};
			for $i ( 0 .. $#{ $optionalElements->{$value} } ) {
				my $elementName = $optionalElements->{$value}[$i];
				$slotFunction = $slotFunction . "   ui->$elementName->hide();\n";
			}
			$slotFunction = $slotFunction . "   if (arg1 == \"$itemKey\") {\n";
			for $i ( 0 .. $#{ $optionalElements->{$value} } ) {
				my $elementName = $optionalElements->{$value}[$i];
				$slotFunction = $slotFunction . "      ui->$elementName->show();\n";
			}
			$slotFunction = $slotFunction . "   }\n";
		}
		$slotFunction = $slotFunction . "}\n";
		$slotFunctions = $slotFunctions . $slotFunction;
	}
	$finalResult;
#	print "SLOT FUNCTIONS:" . $slotFunctions . "\n";
#	print "SLOT FUNCTIONS DEF:" . $slotFunctionsDef . "\n";
#	print "SLOT INIT CALLS:" . $slotCall . "\n";
	if ($part == 1) {
		return ($slotFunctions, $slotCall);
	}
	else {
		return $slotFunctionsDef;
	}
}
#---------------FIN CONDITIONAL COMBOS---------------------


#---------------INPUT ELEMENTS SLOTS------------------------
sub inputElementSlots {
	local($inputLine, $gidx, $sgidx, $ilidx, $moduleClass, $part) = @_;

	my $inputElementData = "";
	$inputElements = $inputLine->{inputElement};
	@inputElementsList = $inputElements;
	if (ref $inputElements eq 'ARRAY') {
		@inputElementsList = @{$inputElements};
	}	
	#creates the ui inputElement data
	$inputElementsData = "";
	$ieidx = 0;
	$editIndex = 1;

	#montage estructuras
	my $slots;
	my $slotsDef;
	foreach $inputElement (@inputElementsList) {
		$type = $inputElement->{inputElementType};
		$inputName = $inputElement->{inputLineEditName};
		$inputElementName = $inputElement->{inputElementName};
		$inputElementIdx = "ie" . $inputName . $gidx . $sgidx . $ilidx . $ieidx;

		#-------Rellenado lista de combos con dependencias
		if ($type eq "list") {
			#listSlots implementation
			open FILE, "<resources/listSlots.cpp";
			my $listSlots = do { local $/; <FILE> };
			close FILE;
			my $subs = $inputElementIdx . "List";
			my $moduleClassUp = ucfirst($moduleClass);
			$listSlots =~ s/#listName/$subs/g;
			$listSlots =~ s/<<moduleClass>>/$moduleClassUp/g;
        #-----init each column separately
        my $items = $inputElement->{item};
        my $itemIdx = 0;
        @itemsList = $items;
        if (ref $items eq 'ARRAY') {
          @itemsList = @{$items};
        }
        my $newItemCols;        
        foreach $item (@itemsList) {          
          my $valueType = $item->{itemValueType};
          print "item value type:" . $valueType . "\n";
          open FILE, "<resources/listNewItem.cpp";
          my $newItem = do { local $/; <FILE> };
          #VALIDATOR CODE
          my $validatorCode;
          if ($valueType eq "REAL") {
            $validatorCode = "realValidator";
          }
          elsif ($valueType eq "INT") {
            $validatorCode = "intValidator";
          }
          elsif ($valueType eq "REALSEQ") {
            $validatorCode = "realseqValidator";
          }
          else {
            $validatorCode = "0";
          }
          print "item validator code:" . $validatorCode . "\n";
          $newItem =~ s/#validator/$validatorCode/g;
          $newItem =~ s/#colNum/$itemIdx/g;
          $newItem =~ s/#tooltip/$valueType/g;                    
          $itemIdx = $itemIdx + 1;
          $newItemCols = $newItemCols . $newItem;
        }
      $listSlots =~ s/<<listColumns>>/$newItemCols/g;
      #-----end init each column separately
			#listSlots declaration
			my $listSlotsDef;
			$listSlotsDef = "void on_". $subs ."Plus_clicked();\n";
			$listSlotsDef = $listSlotsDef . "void on_". $subs ."Minus_clicked();\n";
			$slots = $slots . $listSlots;
			$slotsDef = $slotsDef . $listSlotsDef;
		}
		#-------------------------------------------
		$ieidx = $ieidx + 1;
	}

	if ($part == 1) {
		return ($slots);
	}
	else {
		return $slotsDef;
	}
}
#---------------FIN LIST BUTTONS SLOTS---------------------

#---------------INPUT ELEMENTS VALIDATIONS------------------------
sub inputElementValidations {
	local($inputLine, $gidx, $sgidx, $ilidx) = @_;

	my $inputElementData = "";
	$inputElements = $inputLine->{inputElement};
	@inputElementsList = $inputElements;
	if (ref $inputElements eq 'ARRAY') {
		@inputElementsList = @{$inputElements};
	}	
	#creates the ui inputElement data
	$ieidx = 0;
	$editIndex = 1;

	#montage estructuras
	my $validations;
	foreach $inputElement (@inputElementsList) {
		$type = $inputElement->{inputElementType};
		$inputName = $inputElement->{inputLineEditName};
		$inputElementName = $inputElement->{inputElementName};
    my $valueType = $inputElement->{inputElementValueType};
		$inputElementIdx = "ie" . $inputName . $gidx . $sgidx . $ilidx . $ieidx;
    my $editName;
    if ($type eq "edit") {
				$editName = $inputElementIdx . "Edit";
    }
    #creates the ui for the edite input type: a line edit without name, only inser commas
		elsif ($type eq "edit2") {
			$editName = $inputElementIdx . "Edit";
		}
		#creates the ui for the combo input type
		elsif ($type eq "combo") {
			$editName = $inputElementIdx . "Combo";
		}
		#creates the ui for the list input type
		elsif ($type eq "list") {
			$editName = $inputElementIdx . "List";
		}
    
    my $validatorCode;
    if ($valueType eq "REAL") {
      $validatorCode = "ui->$editName->setValidator(realValidator);\n";
    }
    elsif ($valueType eq "INT") {
      $validatorCode = "ui->$editName->setValidator(intValidator);\n";
    }
    elsif ($valueType eq "REALSEQ") {
      $validatorCode = "ui->$editName->setValidator(realseqValidator);\n";
    }
    $validations = $validations . $validatorCode;		
		$ieidx = $ieidx + 1;
	}
  return $validations;
}
#---------------FIN INPUT ELEMENTS VALIDATIONS---------------------

#----------------------------------------------------------------------
#--Function createCPPFile
#--Description:
#----Generates the c++ code functionality for a module
#--Parameters:
#----moduleName: The module name
#----module: The module metadata
#--Returns: the module c++ file part
sub createCPPFile {
	local($moduleName, $module) = @_;
	#open the module c++ template
	open FILE, "<resources/module.cpp";
	my $moduleCPP = do { local $/; <FILE> };
	close FILE;
	#replace the template tokens by the values
	$moduleCPP =~ s/<<moduleName>>/$moduleName/g;
	my $moduleNameUC = ucfirst($moduleName);
	$moduleCPP =~ s/<<className>>/$moduleNameUC/g;
	my $fileExtension = $module->{fileExtension};
	$moduleCPP =~ s/<<fileExtension>>/$fileExtension/g;
	#----Create the save statements
	my $savesData = "";
        #---------------------------------------------------------------------------------
	$groups = $module->{group};
	@groupsList = $groups;
	if (ref $groups eq 'ARRAY') {
		@groupsList = @{$groups};
	}
	my $gidx = 0;
	my $withOupos = 0;
	my $withBouco = 0;
	foreach $group (@groupsList) {
		$groupName = $group->{groupName};
		$groupType = $group->{groupType};
		#outpos, special group type "outpout and postprocess"
		if ($groupType eq "oupos") {
			$savesData = $savesData . createOuPosSave($group);
			$withOupos = 1;
		}
		elsif ($groupType eq "bouco") {
			$savesData = $savesData . createBouCoSave($group);
			$withBouco = 1;
		}
		#default group type
		else {
			$savesData = $savesData . createSavePart($groupName, $group, $gidx);
		}		
		$gidx = $gidx + 1;
	}
	#---------------------------------------------------------------------------------

	$moduleCPP =~ s/<<uisaves>>/$savesData/g;
	#-----------------------------PRIVATE SLOTS-----------------------------------------
	#----Create the checkBox action functions
	#This are event handlers for each checkBox of the inputLine, each handler enables or disables the corresponding inputLine
	my $slotsData = "";
	my $condSlotsDef = "";
  my $validationsCode = "";
	#---------------------------------------------------------------------------------
	$groups = $module->{group};
	@groupsList = $groups;
	if (ref $groups eq 'ARRAY') {
		@groupsList = @{$groups};
	}
	my $gidx = 0;
	foreach $group (@groupsList) {
		$subGroups = $group->{subGroup};
		@subGroupsList = $subGroups;
		if (ref $subGroups eq 'ARRAY') {
			@subGroupsList = @{$subGroups};
		}	
		#iterate over each subGroup
		$sgidx = 0;
		foreach $subGroup (@subGroupsList) {
			@inputLinesList = getInputLines($subGroup);	
			#iterate over each inputLine
			my $inputLineData = "";
			my $ilidx = 0;
			foreach $inputLine (@inputLinesList) {
				if ($inputLine) {			
					$inputLineName = $inputLine->{inputLineName};
					$inputLineIdx = $inputLineName . $gidx . $sgidx . $ilidx;
					#creates the handler for the inputLine checkbox
					if ($inputLine->{inputElement}) {
						open FILE, "<resources/checkClicked.cpp";
						my $slot = do { local $/; <FILE> };
						close FILE;
						$slot =~ s/<<className>>/$moduleNameUC/g;
						$slot =~ s/<<lineInputName>>/$inputLineIdx/g;
						$subs = $inputLineIdx . "Frame";
						$slot =~ s/<<frameName>>/$subs/g;
						$slotsData = $slotsData . $slot;
					}
					#creates the handler for the inputLine help button
					if ($inputLine->{inputLineHelp}) {
						#help action button
						open FILE, "<resources/helpButton.cpp";
						my $slot = do { local $/; <FILE> };
						close FILE;
						$slot =~ s/<<className>>/$moduleNameUC/g;
						$slot =~ s/<<lineInputName>>/$inputLineIdx/g;
						$helpText = $inputLine->{inputLineHelp};
						$helpText =~ s/\"/\\\"/g;
						#If text has multiple lines.
						my @lines = split /\n/, $helpText;
						my $helpTextAux;
						foreach my $line (@lines) {
							#Quitamos los espacios iniciales de la linea
							$line =~ s/^\s+//;
							#encapsulamos cada linea entre comillas
							$helpTextAux = $helpTextAux ."\"" . $line . "\"\n";
						}
						$slot =~ s/<<helpText>>/$helpTextAux/g;
						$slotsData = $slotsData . $slot;
					}
					#--------------------------------------------------------------
					#------------------------COMBOS CONDICIONALES------------------
					#--------------------------------------------------------------
					my ($one, $two) = conditionalCombos($inputLine, $gidx, $sgidx, $ilidx, $moduleName, 1);
					$slotsData = $slotsData . $one;
					$condSlotsDef = $condSlotsDef . $two;
					#--------------------------------------------------------------
					#------------------------ELEMENT SLOTS--------------
					#--------------------------------------------------------------			        
					my ($one) = inputElementSlots($inputLine, $gidx, $sgidx, $ilidx, $moduleName, 1);
					$slotsData = $slotsData . $one;
          #--------------------------------------------------------------
					#------------------------ELEMENT VALIDATIONS--------------
					#--------------------------------------------------------------			        
					my $one = inputElementValidations($inputLine, $gidx, $sgidx, $ilidx);
					$validationsCode = $validationsCode . $one;
					#--------------------------------------------------------------
					#------------------------FIN ELEMENT VALIDATIONS--------------
					#--------------------------------------------------------------	
					$ilidx = $ilidx + 1;
				}				
			}
			$sgidx = $sgidx + 1;
		}
		$gidx = $gidx + 1;
	}

	#-----------------event slots--------------------
	#add the oupos events
        if ($withOupos) {
		open FILE, "<resources/ouposSlots.cpp";
		my $ouposSlots = do { local $/; <FILE> };
		close FILE;
		$ouposSlots =~ s/<<moduleClass>>/$moduleNameUC/g;
		$slotsData = $slotsData . $ouposSlots;
	}
	#add the bouco events
        if ($withBouco) {
		open FILE, "<resources/boucoSlots.cpp";
		my $ouposSlots = do { local $/; <FILE> };
		close FILE;
		$ouposSlots =~ s/<<moduleClass>>/$moduleNameUC/g;
		$slotsData = $slotsData . $ouposSlots;
	}
	$moduleCPP =~ s/<<privateSlots>>/$slotsData/g;

	#------------loadSave gui values part
	my $loadSave = "";
	my $saveGui = "";
	my $loadGui = "";
	if ($withBouco) {
		open FILE, "<resources/boucoLoadSaveGui.cpp";
		my $loadSaveGui = do { local $/; <FILE> };
		close FILE;
		$loadSaveGui =~ s/<<moduleClass>>/$moduleNameUC/g;
		$loadSave = $loadSave . $loadSaveGui;
		$saveGui = "saveBoundCondGui(filePath);";
		$loadGui = "loadBoundCondGui(filePath);";
	}
	$moduleCPP =~ s/<<loadSaveGui>>/$loadSave/g;
	$moduleCPP =~ s/<<saveGui>>/$saveGui/g;
	$loadGui = $loadGui . $condSlotsDef;
	$moduleCPP =~ s/<<loadGui>>/$loadGui/g;

	#------------loadCodes part
	my $loadCodes = "";
	if ($withBouco) {
		open FILE, "<resources/loadCodes.cpp";
		my $loadCodesFile = do { local $/; <FILE> };
		close FILE;
		$loadCodes = $loadCodes . $loadCodesFile;
	}
	$moduleCPP =~ s/<<loadCodes>>/$loadCodes/g;

	#----------calls in the module constructor
	$moduleCPP =~ s/<<constructorCalls>>/$condSlotsDef/g;
  
  #---------code to validate inputs
  $moduleCPP =~ s/<<setValidators>>/$validationsCode/g;
  
	

	#
	#----Add the modules and services groups that have to be disabled or enabled:
	#for each subGroup of type "module" or "service" creates the code to disable this module
	my $modToDisable = "";
        #---------------------------------------------------------------------------------
	$groups = $module->{group};
	@groupsList = $groups;
	if (ref $groups eq 'ARRAY') {
		@groupsList = @{$groups};
	}
	my $gidx = 0;
	foreach $group (@groupsList) {
		$subGroups = $group->{subGroup};
		@subGroupsList = $subGroups;
		if (ref $subGroups eq 'ARRAY') {
			@subGroupsList = @{$subGroups};
		}	
		$sgidx = 0;
		foreach $subGroup (@subGroupsList) {
			$subGroupType = $subGroup->{subGroupType};
			$subGroupName = $subGroup->{subGroupName};
			$subGroupIdx = "sg" . $subGroupName . $gidx . $sgidx;
			if ($subGroupType eq "module" || $subGroupType eq "service") {
				open FILE, "<resources/disableSubGroups.cpp";
				my $disableSubGroups = do { local $/; <FILE> };
				close FILE;
				$disableSubGroups =~ s/<<subGroupIdx>>/$subGroupIdx/g;
				$modToDisable = $modToDisable . $disableSubGroups;
			}
			$sgidx = $sgidx + 1;
		}
		$gidx = $gidx + 1;
	}
	#---------------------------------------------------------------------------------
	$moduleCPP =~ s/<<modulesToDisable>>/$modToDisable/g;
	
	return $moduleCPP;
}

#----------------------------------------------------------------------
#--Function createHFile
#--Description:
#----Generates the header code functionality for a module
#--Parameters:
#----moduleName: The module name
#----module: The module metadata
#--Returns: the module file header part
sub createHFile {

	local($moduleName, $module) = @_;
	#open the header file template
        open FILE, "<resources/module.h";
	my $moduleH = do { local $/; <FILE> };
	close FILE;
	#replaces the template tokens by values
	my $moduleNameUC = uc($moduleName);
	$moduleH =~ s/<<defname>>/$moduleNameUC/g;
	my $moduleNameUCF = ucfirst($moduleName);
	$moduleH =~ s/<<className>>/$moduleNameUCF/g;
	$groups = $module->{group};
	@groupsList = $groups;
	if (ref $groups eq 'ARRAY') {
		@groupsList = @{$groups};
	}
	#Variabe to create the code to declare the event handlers of the module
	my $slotsData = "";
	my $gidx = 0;
	my $withOupos = 0;
	my $withBouco = 0;
	foreach $group (@groupsList) {
		$subGroups = $group->{subGroup};
		$groupType = $group->{groupType};
		if ($groupType eq "oupos") {
			$withOupos = 1;
		}
		if ($groupType eq "bouco") {
			$withBouco = 1;
		}
		
		@subGroupsList = $subGroups;
		if (ref $subGroups eq 'ARRAY') {
			@subGroupsList = @{$subGroups};
		}	
		#iterates over each subGroup
		my $sgidx = 0;
		foreach $subGroup (@subGroupsList) {
			@inputLinesList = getInputLines($subGroup);	
			#iterates over each inputLine
			my $inputLineData = "";
			my $ilidx = 0;
			foreach $inputLine (@inputLinesList) {
				#creates the declaration of the checkbox and help event handlers
				if ($inputLine) {
					$inputLineName = $inputLine->{inputLineName};
					$inputLineIdx = $inputLineName . $gidx . $sgidx . $ilidx;
					if ($inputLine->{inputElement}) {
						$slot = "    void on_" . $inputLineIdx . "CheckBox_toggled(bool checked);\n";
						$slotsData = $slotsData . $slot;
					}
					if ($inputLine->{inputLineHelp}) {
						$slot = "    void on_" . $inputLineIdx . "HelpButton_clicked();\n";
						$slotsData = $slotsData . $slot;
					}
					#--------------------------------------------------------------
					#------------------------COMBOS CONDICIONALES------------------
					#--------------------------------------------------------------
					my $slotsCombos = conditionalCombos($inputLine, $gidx, $sgidx, $ilidx, $moduleName, 2);
					$slotsData = $slotsData . $slotsCombos;
					#--------------------------------------------------------------
					#------------------------FIN COMBOS CONDICIONALES--------------
					#--------------------------------------------------------------
					my $slotsCombos = inputElementSlots($inputLine, $gidx, $sgidx, $ilidx, $moduleName, 2);
					$slotsData = $slotsData . $slotsCombos;
					#----------------------------------------------------------------
					$ilidx = $ilidx + 1;
				}
			}
			$sgidx = $sgidx + 1;
		}
		$gidx = $gidx + 1;
	}
	
	#private slots--------------------
	#add the oupos events
        if ($withOupos) {
		open FILE, "<resources/ouposSlots.h";
		my $ouposSlots = do { local $/; <FILE> };
		close FILE;
		$slotsData = $slotsData . $ouposSlots;
	}
	#add the bouco events
        if ($withBouco) {
		open FILE, "<resources/boucoSlots.h";
		my $boucoSlots = do { local $/; <FILE> };
		close FILE;
		$slotsData = $slotsData . $boucoSlots;
	}

	$moduleH =~ s/<<privateSlots>>/$slotsData/;


	#private methods--------------------
	$privateMethods = "";
        if ($withBouco) {
		$privateMethods = $privateMethods . "    void insertTimeFunction(int code);\n";
		$privateMethods = $privateMethods . "    void saveBoundCondGui(QString filePath);\n";
		$privateMethods = $privateMethods . "    void loadBoundCondGui(QString filePath);\n";
		$privateMethods = $privateMethods . "    void addbcCNRow(QString codes);\n";
		$privateMethods = $privateMethods . "    void addbcCBRow(QString codes);\n";
	}
	$moduleH =~ s/<<privateMethods>>/$privateMethods/;
	
	return $moduleH;	
}

#----------------------------------------------------------------------
#--Function createAlyawidget
#--Description:
#----Generates the header code, the c++ code and the ui code of the alyaWidget module.
#----This is the main module of the application, it contains the other gui modules.
#--Parameters:
#----modulesList: The modules list
#--Returns: the module file header part
sub createAlyawidget {
	local(@modulesList) = @_;

	#---------------------alyawidget cpp file
	#open the alyawidget cpp template file
	open FILE, "<resources/main/alyawidget.cpp";
	my $alyawidget = do { local $/; <FILE> };
	close FILE;
	$modulesActivation = "";
	$modulesKernelActivation;
	$newGUI;
	$disableModulesServicesParams = "";
	$enableModule = "";
	$enableService = "";
	foreach $module (@modulesList) {	
		$moduleName = $module->{moduleName};
		$isModuleKernel = $module->{isModuleKernel};
		my $moduleNameUC = ucfirst($moduleName);
		
		#modules activation
		#if its a kernelModule follows a special activation and sincronitzacion functions of his subGroups 			
		if ("yes" eq $isModuleKernel) {
			open FILE, "<resources/main/activeKernelModule.cpp";
			my $activeModule = do { local $/; <FILE> };
			close FILE;
			$activeModule =~ s/<<moduleName>>/$moduleName/g;
			$modulesKernelActivation = $modulesKernelActivation . $activeModule;
			#code to disable his modules and services
			$disableModulesServicesParams = $disableModulesServicesParams . "     ui->$moduleName->disableModulesServicesParams();\n";
			#code to enable a module
			$enableModule = $enableModule . "        ui->$moduleName->enableModule(text);\n";
			#code to enable a service
			$enableService = $enableService . "        ui->$moduleName->enableService(text);\n";
		}
		#normal module activation code
		else {
			open FILE, "<resources/main/activeModule.cpp";
			my $activeModule = do { local $/; <FILE> };
			close FILE;
			$activeModule =~ s/<<moduleName>>/$moduleName/g;
			$modulesActivation = $modulesActivation . $activeModule;
		}		

		#newGUI: Code to create a new empty gui dinamically
		open FILE, "<resources/main/newModule.cpp";
		my $newModule = do { local $/; <FILE> };
		close FILE;
		$newModule =~ s/<<moduleName>>/$moduleName/g;
		$newModule =~ s/<<moduleClass>>/$moduleNameUC/g;
		$newGUI = $newGUI . $newModule;		


	}	
	$alyawidget =~ s/<<activeModule>>/$modulesActivation/;
	$alyawidget =~ s/<<activeKernelModule>>/$modulesKernelActivation/;
	$alyawidget =~ s/<<newModule>>/$newGUI/;
	$alyawidget =~ s/<<disableModulesServicesParams>>/$disableModulesServicesParams/;
	$alyawidget =~ s/<<enableModule>>/$enableModule/;
	$alyawidget =~ s/<<enableService>>/$enableService/;

	#save the module to a file
	open FILE, ">$basePath/alyawidget.cpp" or die $!;
	print FILE $alyawidget;
	close FILE;

	#----------------------------alyawidget ui file
	#opens the alyawidget ui template
	open FILE, "<resources/main/alyawidget.ui";
	my $alyawidgetUI = do { local $/; <FILE> };
	close FILE;
	$alyaModules;
	$alyaTabs;
	$customWidgets;
	foreach $module (@modulesList) {	
		$moduleName = $module->{moduleName};
		$isModuleKernel = $module->{isModuleKernel};	
		my $moduleNameUC = ucfirst($moduleName);

		#alyaModules
		#if is not a kernelModule add the module name to the modules list
		if ("yes" ne $isModuleKernel) {
			open FILE, "<resources/main/item.ui";
			my $item = do { local $/; <FILE> };
			close FILE;
			$item =~ s/<<text>>/$moduleName/g;
			$alyaModules = $alyaModules . $item;
		}

		#alyaTabs
		#create one tab for each module
		open FILE, "<resources/main/alyaTab.ui";
		my $alyaTab = do { local $/; <FILE> };
		close FILE;
		$alyaTab =~ s/<<moduleName>>/$moduleName/g;
		$alyaTab =~ s/<<moduleClass>>/$moduleNameUC/g;
		$alyaTabs = $alyaTabs . $alyaTab;

		#customWidgets
		#for each module creates a customWidget: a widget of the module type
		open FILE, "<resources/main/customWidget.ui";
		my $customWidget = do { local $/; <FILE> };
		close FILE;
		$customWidget =~ s/<<moduleName>>/$moduleName/g;
		$customWidget =~ s/<<moduleClass>>/$moduleNameUC/g;
		$customWidgets = $customWidgets . $customWidget;
	}
	$alyawidgetUI =~ s/<<alyaModules>>/$alyaModules/;
	$alyawidgetUI =~ s/<<alyaTabs>>/$alyaTabs/;
	$alyawidgetUI =~ s/<<customWidgets>>/$customWidgets/;

	open FILE, "<resources/main/item.ui";
	my $item = do { local $/; <FILE> };
	close FILE;
	$item =~ s/<<text>>/parall/g;
	$alyawidgetUI =~ s/<<alyaServices>>/$item/;
	#save the module to a file
	open FILE, ">$basePath/alyawidget.ui" or die $!;
	print FILE $alyawidgetUI;
	close FILE;
	
	#copy the h file
	copy("resources/main/alyawidget.h","$basePath/alyawidget.h") or die "Copy failed: $!";	
}

#----------------------------------------------------------------------
#--Function createCMakeFile
#--Description:
#----creates the cMakeLists.txt file used to compile the application with the make command
#--Parameters:
#----modulesList: The modules list
sub createCMakeFile {
	local(@modulesList) = @_;
	#opens the CMakeLists.txt template
	open FILE, "<resources/main/CMakeLists.txt";
	my $cmakelists = do { local $/; <FILE> };
	close FILE;
	$hModuleFiles;
	$uiModuleFiles;
	$sources;
	$modulesCount = @modulesList;
	$modulesIdx = 1;
	#for each module creates a reference to its ui, h, and cpp file.
	foreach $module (@modulesList) {	
		$moduleName = $module->{moduleName};
		
		$hModuleFiles = $hModuleFiles . $moduleName . "/" . $moduleName . ".h";
		$uiModuleFiles = $uiModuleFiles . $moduleName . "/" . $moduleName . ".ui";
		$sources = $sources . $moduleName . "/" . $moduleName . ".cpp\n" . $moduleName . "/" . $moduleName . ".h";

		if ($modulesIdx < $modulesCount) {
			$hModuleFiles = $hModuleFiles . "\n";
			$uiModuleFiles = $uiModuleFiles . "\n";
			$sources = $sources . "\n";
		}
		$modulesIdx = $modulesIdx + 1;
	}	
	$cmakelists =~ s/<<hModuleFiles>>/$hModuleFiles/;
	$cmakelists =~ s/<<uiModuleFiles>>/$uiModuleFiles/;
	$cmakelists =~ s/<<sources>>/$sources/;

	#save the module to a file
	open FILE, ">$basePath/CMakeLists.txt" or die $!;
	print FILE $cmakelists;
	close FILE;
}
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#--------------------------------MAIN-----------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#xml data containing the gui description
$xml = new XML::Simple;
$data = $xml->XMLin("alyaGuiAuto.xml");
rmtree($basePath );
mkdir($basePath );

$modules = $data->{module};
@modulesList = $modules;
if (ref $modules eq 'ARRAY') {
	@modulesList = @{$modules};
}

#Iterate over each module
foreach $module (@modulesList) {	
	$moduleName = $module->{moduleName};
	mkdir($basePath . "/" . $moduleName );
	
	#-------------------------------------------------------------------------------
        #---------------------------create the UI file----------------------------------
	#-------------------------------------------------------------------------------
	my $uiFile = createUIFile ($moduleName, $module);
	#save the module ui file
	open FILE, ">$basePath/$moduleName/$moduleName.ui" or die $!;
	print FILE $uiFile;
	close FILE;

        #-------------------------------------------------------------------------------
        #---------------------------Create the h file-----------------------------------
	#-------------------------------------------------------------------------------
	my $hFile = createHFile ($moduleName, $module);
	#save the module to a file
	open FILE, ">$basePath/$moduleName/$moduleName.h" or die $!;
	print FILE $hFile;
	close FILE;

        #-------------------------------------------------------------------------------
        #---------------------------Create the cpp file-----------------------------------
	#-------------------------------------------------------------------------------
	my $CPPFile = createCPPFile ($moduleName, $module);
	#save the module to a file
	open FILE, ">$basePath/$moduleName/$moduleName.cpp" or die $!;
	print FILE $CPPFile;
	close FILE;	
}

#-------------------------------------------------------------------------------
#---------------------------Create the main gui files---------------------------
#-------------------------------------------------------------------------------

#create the alyawidget files
createAlyawidget(@modulesList);
#create the cmake file
createCMakeFile(@modulesList);
#copy main files
copy("resources/main/alya.qrc","$basePath/alya.qrc") or die "Copy failed: $!";
copy("resources/main/help.img","$basePath/help.img") or die "Copy failed: $!";
copy("resources/main/module.h","$basePath/module.h") or die "Copy failed: $!";
copy("resources/main/myMainWindow.cxx","$basePath/myMainWindow.cxx") or die "Copy failed: $!";
copy("resources/main/myMainWindow.h","$basePath/myMainWindow.h") or die "Copy failed: $!";
copy("resources/main/myMainWindow.ui","$basePath/myMainWindow.ui") or die "Copy failed: $!";
copy("resources/main/ParaViewFilters.xml","$basePath/ParaViewFilters.xml") or die "Copy failed: $!";
copy("resources/main/ParaViewReaders.xml","$basePath/ParaViewReaders.xml") or die "Copy failed: $!";
copy("resources/main/ParaViewSources.xml","$basePath/ParaViewSources.xml") or die "Copy failed: $!";
copy("resources/main/ParaViewWriters.xml","$basePath/ParaViewWriters.xml") or die "Copy failed: $!";
copy("resources/plus.png","$basePath/plus.png") or die "Copy failed: $!";
copy("resources/minus.png","$basePath/minus.png") or die "Copy failed: $!";
exit 0
