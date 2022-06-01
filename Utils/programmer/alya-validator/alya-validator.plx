#!/usr/bin/perl
# configuration
use FindBin;                # locate this script
use lib "$FindBin::Bin/libraries";  # use the parent directory
use File::Copy;
use File::Path;
use File::Compare;
use File::Copy::Recursive qw(dirmove);
use File::Find;
use XML::Simple;
use Data::Dumper;
use POSIX;


#---------------functions------------------------

sub getVars {
        local($line) = @_;

	my @locvars;
	# do this if line is not starting with a ! meaning a fortran comment
	if ($line!~/^\!/) {
		$line =~ s/.*:://;
		@vars = split(",", $line);
		for $i ( 0 .. $#vars ) {
			$var = $vars[$i];
			#remove spaces
			$var =~ s/^\s+//;
			$var =~ s/\s+$//;
			#remove post data
			if($var =~ m/^(\w+)\(*/) {
				push @locvars, $1;
			}
		}
	}
	
	return @locvars;	
}


sub haveDocumentation {
	my($lines, $lineIdx) = @_;
  	
	my $warnings;
	my $finish = 0;
	my $finded = 0;
	my $idx = $lineIdx;
	my $commentPart;
	#!> @addtogroup groupName
	while (!$finish && $idx >= 0) {
		$line = $lines->[$idx];
		$commentPart = $commentPart . $line . "\n";
		#si es lo que buscamos salimos
		if ($line =~ m/^\t*\s*\!\>\s*\@addtogroup\s+\w+/) {
			$finded = 1;
			$finish = 1;
		}
		#si no es una linea de comentario o una linea en blanco salimos
		elsif (!($line =~ m/^\t*\s*\!/) && !($line =~ m/^\t*\s*$/)) {
			$finish = 1;
		}

		if (!$finish) {
			$idx = $idx - 1;
		}		
	}
	
	if (!$finded) {
		$warnings = "The \@addtogroup comment is not present\n";
	}
	#en el trozo de código que hemos obtenido checkeamos el resto de tags
	if (!($commentPart =~ m/^\t*\s*\!\>\s*\@\{/m)) {
		$warnings = $warnings . "The \@{ comment is not present\n";
	}
	if (!($commentPart =~ m/^\t*\s*\!\>\s*\@file\s+\w+/m)) {
		$warnings = $warnings . "The \@file comment is not present\n";
	}
	if (!($commentPart =~ m/^\t*\s*\!\>\s*\@author\s+\w+/m)) {
		$warnings = $warnings . "The \@author comment is not present\n";
	}
	if (!($commentPart =~ m/^\t*\s*\!\>\s*\@date\s+\w+/m)) {
		$warnings = $warnings . "The \@date comment is not present\n";
	}
	if (!($commentPart =~ m/^\t*\s*\!\>\s*\@brief\s+\w+/m)) {
		$warnings = $warnings . "The \@brief comment is not present\n";
	}
	if (!($commentPart =~ m/^\t*\s*\!\>\s*\@details/m)) {
		$warnings = $warnings . "The \@details comment is not present\n";
	}
	if (!($commentPart =~ m/^\t*\s*\!\>\s*\@}/m)) {
		$warnings = $warnings . "The \@} comment is not present\n";
	}

	return $warnings;
}

sub getSubroutineParams {
	my($lines, $lineIdx) = @_;
  	
	my $warnings;
	my $finish = 0;
	my $finded = 0;
	my $idx = $lineIdx;
	my $max = $lineIdx + 20;
	my @parameters;
	# get the parameters list
        # iterate over each parameter subroutine line
	while (!$finish && $idx < $max) {
		$line = $lines->[$idx];
		while($line =~ m/(\w+)\s*,/g) {
			push @parameters, $1;
		}
		if ($line =~ m/(\w+)\s*\)/) {
			push @parameters, $1;
		}
		if ($line =~ m/(\w+)\s*\&/) {
			push @parameters, $1;
		}
		if ($line =~ m/\)/) {
			$finish = 1;
		}
		if (!$finish) {
			$idx = $idx + 1;
		}		
	}

	return @parameters;
}

#dada una declaración de este tipo, en fortran, extrae la lista de variables:
=pod
  character(150)           :: &
       fil_conve_nsi,         &      ! Convergence file name
       fil_dynin_nsi,         &      ! Input dynamic model
       fil_dynou_nsi                 ! Output dynamic model
=cut
#Se recorre fila a fila hasta que encuentra el fin de la sentencia(ausencia del &)
sub getDeclarationVars {
        my($lines, $lineIdx) = @_;
	
	my $idx = $lineIdx;
	my $finish = 0;
	#Medida de seguridad, como máximo declaraciones de 5000 lineas
	my $max = $lineIdx + 5000;
	#array donde se van metiendo las variables encontradas
	my @locvars;
	while (!$finish && $idx < $max) {
		$line = $lines->[$idx];
		#Si la linea no es un comentario la procesamos
		if ($line!~/^\!/) {
			#quitamos toda la parte de la izquierda a partir de ::
			$line =~ s/.*:://;
			#quitamos los posibles comentarios del final
			$line =~ s/\!.*//;
			#quitamos los parentesis
			$line =~ s/\(.*?\)//g;
			#partimos las variables por la ","
			@vars = split(",", $line);
			for $i ( 0 .. $#vars ) {
				#procesamos cada variable
				$var = $vars[$i];
				#le quitamos el texto que tenga a la derecha
				#remove spaces
				$var =~ s/^\s+//;
				$var =~ s/\s+$//;
				#remove post data
				if($var =~ m/^(\w+)\(*/) {
					#guardamos la variable en el array
					push @locvars, $1;
				}
			}
		}
		#Detección de final de sentencia
		#La sentencia acaba cuando no hay "&" al final, detras puede haber o no comentarios
		if ($line !~ m/(\&\s*(\!\s*.*)?)$/) {
			$finish = 1;
		}
		$idx = $idx + 1;
	}
	#print "Vars:" . join(", ", @locvars). "\n";
	return (\@locvars, $idx);	
}

sub excluded {
  my($data, $module, $var) = @_;
  my $excluded = 0;
  
  #Iterate over each module
	$modulesList = $data->{module};
	@moduleList = $modulesList;
	if (ref $modulesList eq 'ARRAY') {
		@moduleList = @{$modulesList};
	}
	#iterate over each test
	foreach $moduleVars (@moduleList) {
		#check if there is an error
		if ($moduleVars->{moduleName} eq $module) {
        #Iterate over each var
        $varsList = $moduleVars->{excludedVar};
        @varList = $varsList;
        if (ref $varsList eq 'ARRAY') {
          @varList = @{$varsList};
        }
        #iterate over each test
        foreach $excVar (@varList) {
          #check if there is an error
          if ($excVar eq $var) {
            $excluded = 1;
            last;
          }				
        }
        last;			
		}				
	}
  return $excluded;
}

sub checkDefVars {
	my($modulesFolder) = @_;

	opendir(DIR, $modulesFolder) or die $!;
	#Iteramos por cada módulo
	my $checkDefVarsResult = {};
	while (my $file = readdir(DIR)) {
		if (-d "$modulesFolder/$file" && $file ne "." && $file ne ".." && $file ne ".svn") {
			#recorremos los ficheros del módulo buscando los def_ y guardando sus variables
			opendir(DIR_MODULE, "$modulesFolder/$file") or die $!;
			#1.Get the variables defined in the def_ f90 files
			my $defVarsResult = {};
			while (my $fileModule = readdir(DIR_MODULE)) {
				if (-f "$modulesFolder/$file/$fileModule" && $fileModule =~ m/^def_.+\.f90$/) {
					open FILE, "<$modulesFolder/$file/$fileModule";
					my $content = do { local $/; <FILE> };
					close FILE;
					# Obtain the vars
					my @lines = split /\n/, $content;
					my $linesLength = scalar (@lines);
					my $lineIdx = 0;
					my $exit = 0;
					my $begin = 0;
					# Recorremos todas las filas
					while ($lineIdx < $linesLength && !$exit) {
						my $line = $lines[$lineIdx];
						#Miramos a ver si ya hemos llegado al punto de inicio
						if ($line =~ m/^\t*\s*\!--BEGIN REA GROUP/) {
							$begin = 1;
						}
						if ($line =~ m/^\t*\s*\!--END REA GROUP/) {
							$exit = 1;
						}
						#entra si la linea es un comienzo de declaración de tipo variable y no es un parámetro
						if ($line =~ m/^\t*\s*(real|integer|type|character)/ && !($line =~ m/parameter/) 
						    && $begin && !$exit) {
							#obtenemos la lista de variables declaradas con ese tipo
							my ($localVars, $idxAux) = getDeclarationVars(\@lines, $lineIdx);
							#ponemos la lista de variables obtenidas en una hash
							foreach my $localVar (@$localVars) {
								if (!exists $defVarsResult->{$localVar}) {
									my $varResult = {};
									$varResult->{found} = "false";
									$varResult->{file} = [];
									$defVarsResult->{$localVar} = $varResult;
								}
								else {
									print "ERROR DUPLICATE VAR IN DEF FILE, MODULE:$file, VAR:$localVar, FILE:$fileModule\n";
								}
							}
							#saltamos a la linea que contiene la siguiente sentencia a la declaración analizada
							$lineIdx = $idxAux;
						}
						else {
							$lineIdx = $lineIdx + 1;
						}						
					}
				}
			}
			closedir(DIR_MODULE);

			#volvemos a recorrer los ficheros del módulo buscando los _sendat o _parall,
			#para cada uno de ellos validamos si la variable está en uno u otro
			opendir(DIR_MODULE, "$modulesFolder/$file") or die $!;
			#2.Check if variables are in _sendat or _parall
			while (my $fileModule = readdir(DIR_MODULE)) {
				if (-f "$modulesFolder/$file/$fileModule" && 
                                    ($fileModule =~ m/_sendat\.f90$/ || $fileModule =~ m/_parall\.f90$/) ) {
					open FILE, "<$modulesFolder/$file/$fileModule";
					my $content = do { local $/; <FILE> };
					close FILE;
					#Nos recorremos las variables y miramos si están en el fichero
					for (keys %$defVarsResult)
					{	my $var = $_;			
    					        my $found = 0;
						while($content =~ m/\W$var\W/g) {
							$found = 1;
						}
						if ($found) {
							my $varResult = $defVarsResult->{$var};
							$varResult->{found} = "true";
							push(@{$varResult->{file}}, $fileModule);
							$defVarsResult->{$var} = $varResult;
						}

				    	}
				}
			}
			closedir(DIR_MODULE);
			$checkDefVarsResult->{$file} = $defVarsResult;
		}
    	}
    closedir(DIR);
    
    #Load the list of excluded vars
    $xml = new XML::Simple;
    $data = $xml->XMLin("defVarsExcluded.xml");
   
    #print the result:
    print "\n-------------------------------------------------------------\n";
    print "-------------------------------------------------------------\n";
    print "---------------Check Vars Result-----------------------------\n";
    print "-------------------------------------------------------------\n";
    print "-------------------------------------------------------------\n\n";
 


   print "----Not foun vars in _sendat and in _parall------------------\n";
	for (keys %$checkDefVarsResult)
	{	my $module = $_;
		print "\n---------------\n";
		print "---------MODULE:$module\n";
		print "----------------\n";
		my $defVarsResult = $checkDefVarsResult->{$module};
		my $notFoundVars;
		for (keys %$defVarsResult)
		{	my $var = $_;
			$result = $defVarsResult->{$var};
			if ($result->{found} ne "true" && !excluded($data, $module, $var)) {
				$notFoundVars = $notFoundVars . $var . " ,";
			}
		}
		print $notFoundVars . "\n";
	}

=pod
	for (keys %$checkDefVarsResult)
	{	my $module = $_;
		print "\n---------------\n";
		print "---------MODULE:$module\n";
		print "----------------\n";
		my $defVarsResult = $checkDefVarsResult->{$module};
		for (keys %$defVarsResult)
		{	my $var = $_;
			my $message = "-----VAR:$var -->";
			$result = $defVarsResult->{$var};
			if ($result->{found} eq "true") {
				$message = $message . "found in:";
				foreach my $file (@{$result->{file}}) {
					$message = $message . $file . ",";
				}
				$message = $message . "\n";
				print $message;
			}
		}
		for (keys %$defVarsResult)
		{	my $var = $_;
			my $message = "-----VAR:$var -->";
			$result = $defVarsResult->{$var};
			if ($result->{found} ne "true") {
				$message = $message . "NOT FOUND!!!";
				$message = $message . "\n";
				print $message;
			}
		}
	}
=cut
}

sub getSubroutineParams {
	my($lines, $lineIdx) = @_;
  	
	my $warnings;
	my $finish = 0;
	my $finded = 0;
	my $idx = $lineIdx;
	my $max = $lineIdx + 20;
	my @parameters;
	# get the parameters list
        # iterate over each parameter subroutine line
	while (!$finish && $idx < $max) {
		$line = $lines->[$idx];
		while($line =~ m/(\w+)\s*,/g) {
			push @parameters, $1;
		}
		if ($line =~ m/(\w+)\s*\)/) {
			push @parameters, $1;
		}
		if ($line =~ m/(\w+)\s*\&/) {
			push @parameters, $1;
		}
		if ($line =~ m/\)/) {
			$finish = 1;
		}
		if (!$finish) {
			$idx = $idx + 1;
		}
	}

	return @parameters;
}

#Devuelve en numero de parametros opcionales
sub getOptionalParams {
	my($lines, $lineIdx, $subrout) = @_;
	my $subroutinePart;
	my $numLines = scalar(@$lines);
	my $idx = $lineIdx;
	my $finish = 0;
	my $optionalArgs = 0;
	while (!$finish && $idx < $numLines) {
		$line = $lines->[$idx];

		#si contiene optional, contamos las variables que son de este tipo
		if ($line =~ m/\Woptional\W/) {
			my @vars = getVars($line);
			$optionalArgs = $optionalArgs + scalar(@vars);
		}

		#si es final de subrutina salimos
		if ($line =~ m/end subroutine $subrout/) {
			$finish = 1;
		}

		if (!$finish) {
			$idx = $idx + 1;
		}		
	}

	return $optionalArgs;
}

sub getCallNumParams {
	my($lines, $lineIdx, $subrName) = @_;

	my $warnings;
	my $finish = 0;
	my $finded = 0;
	my $idx = $lineIdx;
	my $max = $lineIdx + 20;
	my @parameters;
	#si no tiene parentesis no hay ningún parametro
        if (($lines->[$idx] =~ m/call\s+$subrName\s*\(/) &&
            !($lines->[$idx] =~ m/call\s+$subrName\s*\(\)/)) {
		# get the parameters list
		# iterate over each parameter subroutine line
		while (!$finish && $idx < $max) {
			$line = $lines->[$idx];
			#print "line:$line\n";
			#quitamos el parentesis inicial
			$line =~ s/call\s+(\w+)\s*\(//;
			#print "line2:$line\n";
			#quitamos las comillas simples
			$line =~ s/\'.*?'/a/g;
			#print "line3:$line\n";
			#quitamos los posibles comentarios del final
			$line =~ s/\!.*//;
			#print "line4:$line\n";
			#quitamos los parentesis y su contenido
			while ($line =~ m/\([^\(]*?\)/) {
				$line =~ s/\([^\(]*?\)//g;
			}
			#print "line5:$line\n";
			#quitamos los simbolos raros menos las comas, parentesis, espacio por a
			$line =~ s/[^,)\s\&]/a/g;
			#print "line6:$line\n";
			while($line =~ m/(\w+)\s*,/g) {
				#print "push1:$1\n";
				push @parameters, $1;
			}
			if ($line =~ m/(\w+)\s*\)/) {
				#print "push2:$1\n";
				push @parameters, $1;
			}
			#if ($line =~ m/(\w+)\s*\&/) {
			#	print "line7:$1\n";
			#	push @parameters, $1;
			#}
			if ($line =~ m/\)/) {
				$finish = 1;
			}
			if (!$finish) {
				$idx = $idx + 1;
			}
		}
	}

	return @parameters;
	
}

sub checkNumArguments{
	my($sourceFolder) = @_;
        my @files;
	find(\&print_name_if_dir, $sourceFolder);

	sub print_name_if_dir
	{
		my $file = $File::Find::name;
		
		if (-f $file && $file =~ m/\.f90$/) {
			push @files, $file;
		}
	}
	#Introducimos en el subroutData la lista de subroutinas con su numero de parámetros
	my $subroutData = {};
	#print "Obtain subroutines num of arguments------------------------\n";
	foreach my $file (@files) {
		#print "PROCESSING $file FILE------------------------\n";
		#open file
		open FILE, $file;
		my $content = do { local $/; <FILE> };
		close FILE;

		my @lines = split /\n/, $content;
		my $lineIdx = 0;
		foreach my $line (@lines) {
			if ($line =~ m/^\t*\s*subroutine\s*(\w+)\s*\(/) {
				my $subrout = $1;
				my @parameters = getSubroutineParams(\@lines, $lineIdx);
				my $optPars = getOptionalParams(\@lines, $lineIdx, $subrout);
				my $data = {};
				$data->{numArgs} = scalar(@parameters);
				$data->{optArgs} = $optPars;
				if (exists $subroutData->{$subrout}) {
					$data->{duplicated} = "true";
				}
				$subroutData->{$subrout} = $data;
			}
			$lineIdx = $lineIdx + 1;
		}
	}
	#Nos recorremos los ficheros y vamos buscando los calls que tiene
	my $checkNumArgsResult = {};
	#print "Check subroutine call parameters------------------------\n";
	foreach my $file (@files) {
		#print "PROCESSING $file FILE------------------------\n";
		#open file
		open FILE, $file;
		my $content = do { local $/; <FILE> };
		close FILE;

		my @lines = split /\n/, $content;
		my $lineIdx = 0;
		my @fileCalls;
		foreach my $line (@lines) {
			if ($line =~ m/call\s+(\w+)/ && $line !~/^\s*\!/) {
				my $callName = $1;
				my @parameters = getCallNumParams(\@lines, $lineIdx, $callName);
				my $dataAux = $subroutData->{$callName};
				if ($dataAux && $callName ne 'runend' && !exists ($dataAux->{duplicated})) {
					my $maxArgs = $dataAux->{numArgs};
					my $minArgs = $dataAux->{numArgs} - $dataAux->{optArgs};
					my $callArgs = scalar(@parameters);
					if (!($maxArgs >= $callArgs && $callArgs >= $minArgs)) {
						$badArgs = {};
						$badArgs->{name} = $callName;
						$badArgs->{maxArgs} = $maxArgs;
						$badArgs->{minArgs} = $minArgs;
						$badArgs->{callArgs} = $callArgs;
						push @fileCalls, $badArgs;
					}
				}
				#else {
				#	print "subroutine $callName not processed\n";
				#}
			}
			$lineIdx = $lineIdx + 1;
		}
		if (scalar(@fileCalls) > 0) {
			$checkNumArgsResult->{$file} = \@fileCalls;
		}
	}

	#print the results
	print "\n-------------------------------------------------------------\n";
	print "-------------------------------------------------------------\n";
	print "---------------Check Num Arguments Result--------------------\n";
	print "-------------------------------------------------------------\n";
	print "-------------------------------------------------------------\n\n";
	for (keys %$checkNumArgsResult)
	{	my $file = $_;
		print "\n---------------\n";
		print "---------File:$file\n";
		print "----------------\n";
		my $fileCalls = $checkNumArgsResult->{$file};
		foreach my $call (@{$fileCalls}) {
			$name = $call->{name};
			$maxArgs = $call->{maxArgs};
			$minArgs = $call->{minArgs};
			$callArgs = $call->{callArgs};
			print "call $name bat number of arguments, max:$maxArgs, min:$minArgs, call:$callArgs\n";
		}		
	}

}

sub getMemchk {
	my($lines, $lineIdx) = @_;
  	
	my $warnings;
	my $finish = 0;
	my $finded = 0;
	my $idx = $lineIdx + 1;
	my $max = $lineIdx + 20;
	# get the parameters list
        # iterate over each parameter subroutine line
	while (!$finish && $idx < $max) {
		$line = $lines->[$idx];
		if ($line =~ m/w+/) {
			$finish = 1;
		}
		if ($line =~ m/call memchk/) {
			$finded = 1;
		}
		$idx = $idx + 1;
	}

	return $finded;
}

sub checkAllocates{
	my($sourceFolder) = @_;
        my @files;
	find(\&print_name_if_dir, $sourceFolder);

	sub print_name_if_dir
	{
		my $file = $File::Find::name;
		
		if (-f $file && $file =~ m/\.f90$/) {
			push @files, $file;
		}
	}
	#Introducimos en el subroutData la lista de subroutinas con su numero de parámetros
	my $subroutData = {};
	my $checkAllocateResult = {};
	#print "Obtain subroutines num of arguments------------------------\n";
	foreach my $file (@files) {
		my @lineError;
		#print "PROCESSING $file FILE------------------------\n";
		#open file
		open FILE, $file;
		my $content = do { local $/; <FILE> };
		close FILE;

		my @lines = split /\n/, $content;
		my $lineIdx = 0;
		foreach my $line (@lines) {
			if ($line =~ m/^\t*\s*allocate\s*\(/ && $line =~ m/stat=istat/) {
				my $subrout = $1;
				my $finded = getMemchk(\@lines, $lineIdx);
				if (!$finded) {
					push @lineError, ($lineIdx + 1);
				}
			}
			$lineIdx = $lineIdx + 1;
		}
		if (scalar(@lineError) > 0) {
			$checkAllocateResult->{$file} = \@lineError;
		}
	}	

	#print the results
	print "\n-------------------------------------------------------------\n";
	print "-------------------------------------------------------------\n";
	print "---------------Check Allocate checkins-----------------------\n";
	print "-------------------------------------------------------------\n";
	print "-------------------------------------------------------------\n\n";
	for (keys %$checkAllocateResult)
	{	my $file = $_;
		print "\n---------------\n";
		print "---------File:$file\n";
		print "----------------\n";
		my $lineError = $checkAllocateResult->{$file};
		foreach my $err (@{$lineError}) {
			print "At line:$err\n";
		}		
	}

}


#------------------------------------------------
#---------------main-----------------------------
#------------------------------------------------
checkDefVars("../../../Sources/modules");
checkNumArguments("../../../Sources");
checkAllocates("../../../Sources");
=pod
#--------Individual files checking
my $changes = `$svnlook changed -r "$rev" "$repos"`;
my @files = split(/\n/, $changes);
my @report;
#process files and make validations
foreach $file (@files) {
	my $file = substr $file, 4;
	#File content
	my $content = `$svnlook cat -r "$rev" "$repos" "$file"`;
	if ($content && ($file =~ m/.f90/)) {
		my $fileWarn = {};
		my @warnList;
		$fileWarn->{file} = $file;
		while($content =~ m/(rexcha\s*\([i|j|k|l|m|n|o|p|q])/g) {
			my $warnings = "It contains a call to rexcha with an integer parameter:$1";
			push @warnList, $warnings;
		}
		while($content =~ m/(iexcha\s*\([a|b|c|d|e|f|g|h|r|s|t|u|v|w|x|y|z])/g) {
			my $warnings = "It contains a call to rexcha with an integer parameter:$1";
			push @warnList, $warnings;
		}
		while($content =~ m/(getint\s*\(.*\))/g) {
			my $sentence = $1;
			#get the second parameter
			if (($sentence =~ m/,\s*(.+)\s*,/)) {
				$aux = $1;
				if (!($1 =~ m/_ip/)) {
					my $warnings = "It contains a call to getint without a _ip parameter:$sentence";
					push @warnList, $warnings;
				}				
			}
			else {
				my $warnings = "It contains a call to getint without a correct second parameter:$sentence";
				push @warnList, $warnings;
			}											
		}
		while($content =~ m/(getrea\s*\(.*\))/g) {
			my $sentence = $1;
			#get the second parameter
			if (($sentence =~ m/,\s*(.+)\s*,/)) {
				$aux = $1;
				if (!($1 =~ m/_rp/)) {
					my $warnings = "It contains a call to getrea without a _rp parameter:$sentence";
					push @warnList, $warnings;
				}				
			}
			else {
				my $warnings = "It contains a call to getrea without a correct second parameter:$sentence";
				push @warnList, $warnings;
			}											
		}		
		if (!($file =~ m/def_/)) {
			my @lines = split /\n/, $content;
			my $lineIdx = 0;
			foreach my $line (@lines) {
				#nullify
				if ($line =~ m/pointer/ && !($line =~ m/intent/) && $line =~ m/real|integer|type/) {
					#Its a pointer declaration, get the vars
					@vars = getVars($line);
					foreach my $var (@vars) {
						#search if have a nullify
						if(!($content =~ m/nullify\s*\(.*$var.*/)) {
							my $warnings = "The pointer $var shoud be nullified: nullify($var);";
							push @warnList, $warnings;
						}
					}
				}
				#doxygen subroutine declaration
				if ($line =~ m/^\t*\s*subroutine\s*(\w+)\s*\(/) {
					my $subrout = $1;
					#check if has header documentation
					my $warnings = haveDocumentation(\@lines, $lineIdx - 1);					
					#check if each subroutine parameter is documented
					my @parameters = getSubroutineParams(\@lines, $lineIdx);
					foreach my $parameter (@parameters) {
						#search if have a comment
						if(!($content =~ m/::\s*$parameter\s*(\(.*\))?\s*\!\<\s*\w+/)) {
							$warnings = $warnings . "The \"$parameter\" parameter has no doxygen comments\n";
						}
					}
					if ($warnings) {
						$warnings = "The \"$subrout\" subroutine has no doxygen comments:\n" . $warnings ;
						push @warnList, $warnings;
					}
				}
				$lineIdx = $lineIdx + 1;
			}	
		}
		#doxygen subroutines documentation
		if (scalar(@warnList) > 0) {
			$fileWarn->{warnings} = [ @warnList ];
			push @report, $fileWarn;
		}
	}
}

#create validations report
$warnings;
if (scalar(@report) > 0) {
	$warnings = "-----------------------------------------------------------------------\n";
	$warnings = $warnings . "-------------------------------WARNINGS--------------------------------\n";
	$warnings = $warnings . "-----------------------------------------------------------------------\n\n";
}
for $i ( 0 .. $#report ) {
     my $file = $report[$i]->{file};
     $warnings = $warnings . "------------------------File:$file\n";
     foreach $j ( 0 .. $#{ $report[$i]->{warnings} } ) {
	$message = $report[$i]->{warnings}[$j];
        $warnings = $warnings . "----WARNING!! $j: $message\n";
     }
     $warnings = $warnings . "-----------------------------------------------------------------------\n";
}

if ($warnings) {
	print STDERR $warnings;
	exit(1);
}
=cut
exit(0);
