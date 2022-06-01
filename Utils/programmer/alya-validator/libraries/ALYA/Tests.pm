package ALYA::Tests;

#---------------------------------------------------------------------
#--Package that make one test for a given module----------------------
#---------------------------------------------------------------------

#usage modules
use FindBin;                 # locate this script
use lib "$FindBin::Bin/.";  # use the parent directory
use File::Copy;
use File::Path;
use File::Compare;
use XML::Simple;
use Data::Dumper;
use Text::Diff;
use POSIX;
use ALYA::ExecFactory;

#global variables
$basePath = "modules/";
$alyaExec = "~/alya/Executables/unix/Alya.x";
$testId = "";
$columnCompared = 6;
$reportDir = "report";

#----------------------------------------------------------------------
#-------------------------------functions------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#--Function compareFiles
#--Description:
#----Makes comparisons between files and creates reports about the results
#--Parameters:
#----comparison: Array that contains the comparisons list to do and his properties
#----reportName: Then name of the report to generate
sub compareFiles {
	local($comparison, $reportName, $result, $processors) = @_;

	#log
	$result->{file} = [];

	#compare each file
	@comparFile = $comparison;
	if (ref $comparison eq 'ARRAY') {
		@comparFile = @{$comparison};
	}
	#iterate over each comparison
	foreach $compar (@comparFile) {
		#log
		$compareLog = {};
		$compareLog->{fileName} = $compar->{file};
		push(@{$result->{file}}, $compareLog);
		#obtain files to compare
		if (-d $basePath . "/base/". $processors ."p") {
			$ficheroBase = $basePath . "/base/". $processors ."p/" . $testId . $compar->{file};
		}
		else {
			$ficheroBase = $basePath . "/base/1p/" . $testId . $compar->{file};
		}

		$ficheroComp = $basePath . "/tmp/" . $testId . $compar->{file};	

		#report the result
		if ($compar->{type} eq "diff") {
			compareByDiff($ficheroBase, $ficheroComp, $compareLog, $processors, $testId, $compar);
		}
		elsif ($compar->{type} eq "value") {
			compareByValue($ficheroBase, $ficheroComp, $compareLog, $processors, $testId, $compar);
		}           		
	}
	          	
}

#----------------------------------------------------------------------
#--Function compareFilesByValue
#--Description:
#----Makes comparisons between files and creates reports about the results
#--Parameters:
#----comparison: Array that contains the comparisons list to do and his properties
#----reportName: Then name of the report to generate
sub compareByDiff {
	local($ficheroBase, $ficheroComp, $compareLog, $processors, $testId, $compar) = @_;
	
	if (compare($ficheroBase,$ficheroComp) == 0) {
		$compareLog->{passed} = true;
	}
	else {
		$compareLog->{passed} = false;
		#check if file exists
		unless (-e $ficheroComp) {
			$compareLog->{error} = "File have not been generated";
		}
		else {
			#make the diff between the files
			my $diff = diff $ficheroBase, $ficheroComp, { STYLE => "Context", CONTEXT=> "0" };
			$compareLog->{diff} = $diff;
			#save the file into the report folder
			copy($ficheroComp , $basePath . "/".$reportDir."/" . $processors . "p_".$testId .$compar->{file})
			    or die "copy failed: $!";
		}
	}
}

#----------------------------------------------------------------------
#--Function compareFilesByValue
#--Description:
#----Makes comparisons between files and creates reports about the results
#--Parameters:
#----comparison: Array that contains the comparisons list to do and his properties
#----reportName: Then name of the report to generate
sub compareByValue {
	local($ficheroBase, $ficheroComp, $compareLog, $processors, $testId, $compar) = @_;


	my $section = $compar->{section};
	my $columnToCompare = $compar->{column};
	my $compareListBase = [];
	my $compareListTest = [];

	$compareLog->{section} = $section;
        #----------------FILL SECTION-------------------------------------------
	#fill the comparison structures with the comparison items from the files
	#-----------------------------------------------------------------------
	if ($section eq "FINITE ELEMENT ERRORS") {
		#search for w(n,n) at start line
		#for each line we get each number and make the comparison
		open my $base, $ficheroBase or die "Could not open $ficheroBase: $!";
		while( my $baseLine = <$base>)  { 
			if ($baseLine =~ /^\s*W\(/s) {
				#get each number of the line
				my @compareLine = split(/\s+/,$baseLine);
				my $compRef = [];
				foreach $line (@compareLine) {
					if ($line =~ /^-*\d/s) {
						push(@{$compRef},$line);
					}
				}
				push(@{$compareListBase}, $compRef);
			}
		}

		open my $base, $ficheroComp;
		while( my $baseLine = <$base>)  {    
			if ($baseLine =~ /^\s*W\(/s) {
				#get each number of the line
				my @compareLine = split(/\s+/,$baseLine);
				my $compRef = [];
				foreach $line (@compareLine) {
					if ($line =~ /^-*\d/s) {
						push(@{$compRef},$line);
					}
				}
				push(@{$compareListTest}, $compRef);
			}
		}
		
	}
	elsif ($section eq "ALYA Witness set results") {
		#obtain the last line, get all the columns and compare his values
		open my $base, $ficheroBase or die "Could not open $ficheroBase: $!";
		my $lastBaseLine;
		while( my $baseLine = <$base>)  {
			$lastBaseLine = $baseLine;
		}
		#get each number of the line
		my @compareLine = split(/\s+/,$lastBaseLine);
		my $compRef = [];
		my $idx = 0;
		print "Column:" . $columnToCompare . "\n"; 
		foreach $line (@compareLine) {
			if ($line =~ /^-*\d/s) {
				print "idx:" . $idx . "\n";				
				if (!$columnToCompare || $columnToCompare==$idx) {
					print "Pushed\n";
					push(@{$compRef},$line);
				}
				$idx = $idx + 1;
			}
		}
		push(@{$compareListBase}, $compRef);
		
		open my $base, $ficheroComp;
		my $lastCompLine;
		while( my $baseLine = <$base>)  {    
			$lastCompLine = $baseLine
		}
		#get each number of the line
		my @compareCompLine = split(/\s+/,$lastCompLine);
		my $compCompRef = [];
		$idx = 0;
		foreach $line (@compareCompLine) {
			if ($line =~ /^-*\d/s) {				
				if (!$columnToCompare || $columnToCompare==$idx) {
					push(@{$compCompRef},$line);
				}
				$idx = $idx + 1;
			}
		}
		push(@{$compareListTest}, $compCompRef);
	}

        #----------------END FILL SECTION---------------------------------------
	#fill the comparison structures with the comparison items from the files
	#-----------------------------------------------------------------------

        #----------------COMPARE SECTION----------------------------------------
	#fill the comparison structures with the comparison items from the files
	#-----------------------------------------------------------------------
	
	#check if have the same number of elements
	if (scalar @{$compareListBase} != scalar @{$compareListTest}) {
		$compareLog->{error} = "File doesn't have the same number of lines, base:" . scalar @{$compareListBase} .
                                       " lines, test:" . scalar @{$compareListTest} . " lines";
	}
	else {
		my $line = 0;
		$compareLog->{diffValue} = [];
		foreach $comparLine (@{$compareListBase}) {
			my $column = 0;
			foreach $lineToCompare (@{$comparLine}) {
				$original = $lineToCompare;
				$obtained = $compareListTest->[$line][$column];
				$difference = $original - $obtained;
				if (isDifferent($difference,$compar->{decimals} )) {
					$diffValue = {};
					if ($columnToCompare) {
						$diffValue->{column} = $columnToCompare;
					}
					else {
						$diffValue->{column} = $column;
					}					
					$diffValue->{original} = $original;
					$diffValue->{obtained} = $obtained;
					push(@{$compareLog->{diffValue}},$diffValue);
				}	
				$column = $column + 1;
			}
			$line = $line + 1;	
		}
		if (scalar @{$compareLog->{diffValue}} > 0) {
			$compareLog->{passed} = "false";
		}
		else {
			$compareLog->{passed} = "true";
		}		
	}		
}

sub isDifferent {
	local($difference, $decimals) = @_;
	$isDifferent = 0;
	$differ = $decimals - $difference;
	if (abs($difference) > $decimals) {
		$isDifferent = 1;
	}
	else {
		$isDifferent = 0;
	}
	return $isDifferent;
}


#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#--------------------------------MAIN-----------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
sub doTests {
	local($common, $moduleName, $testName) = @_;
        $basePath = "modules/";
	#xml data
	$xml = new XML::Simple;
	$data = $xml->XMLin("modules/" . $moduleName . "/" . $testName . ".xml");

	$execModule = $common->{execModule};
	if (!$execModule) {
		$execModule = "LocalExec";
	}
	$execution = ALYA::ExecFactory->instantiate($execModule);


	#global initializations
	$basePath = $basePath . $moduleName . "/" . $data->{id};
	$alyaExec = $common->{executable};
	$alya2posExec = $common->{alya2pos};
	$testId = $data->{id};
	rmtree($basePath . "/tmp");

        #report xml initialization
        rmtree( $basePath . "/" . $reportDir );
        mkdir( $basePath . "/" . $reportDir);
	
	my $report = {};
	$report->{testResult} = {};
        $report->{testResult}->{name} = $data->{id};
        $report->{testResult}->{module} = $moduleName;
        $report->{testResult}->{services} = $data->{services};

	$report->{testResult}->{results}->{result} = [];
        my $results = $report->{testResult}->{results};

	#Number of tests
	$processorsTest = $data->{numProcessors}->{numProcessor};
	@processors = $processorsTest;
	if (ref $processorsTest eq 'ARRAY') {
		@processors = @{$processorsTest};
	}
	#iterate over each test
	foreach $numProcessor (@processors) {

		#log
		$result = {};
                $result->{processors} = $numProcessor;
                push(@{$results->{result}}, $result);

		#1.create tmp dir
		mkdir( $basePath . "/tmp") || die "Unable to create directory <$!>\n";

		#exec locally
		#2.copy input files
		opendir DH, $basePath or die "Couldn't open the current directory: $!";
		while ($_ = readdir(DH)) {
			next if $_ eq "." or $_ eq ".." or -d $_;
			if (!($_ eq "base") && !($_ eq $reportDir) && !($_ eq "tmp")) {
				copy($basePath . "/" . $_ , $basePath . "/tmp" . "/" . $_)
				    or die "copy failed: $!";		
			}
		}


		$execOutput = $execution->execTest($basePath, "../../../../",  $numProcessor, $alyaExec, $testId, $data->{timeout});


		#4.check the execution finalization
		print "Exec output:" . $execOutput;
		if ($execOutput =~ m/ALYA  FINISHED NORMALLY/) {	
		   print "finished correctlyyyyyyyyyyyyyyyyyyyyy\n";
		   #if it is required executes alya2pos
		   if ("true" eq $data->{comparisons}->{withAlya2pos}) {
		   	$execution->execWithOpen($basePath, "../../../../",  "$alya2posExec " . $testId, 600);
		   }
		   #compare the files and create the report
		   compareFiles($data->{comparisons}->{comparison}, $reportNameFile, $result, $numProcessor);
		}
		else {
		   print "finished errorssssssssssssssssssssss\n";
		   #log
		   $result->{file} = [];
	           $compareLog = {};
		   $compareLog->{fileName} = "";
		   $compareLog->{error} = "Alya execution finished with errors:\n" . $execOutput;
                   $compareLog->{passed} = "false";
		   push(@{$result->{file}}, $compareLog);
		}
		#5.delete temporal files
		$files_deleted = rmtree($basePath . "/tmp");	

	}

	#write the report
	$xml->XMLout(
	    $report,
	    KeepRoot => 1,
	    NoAttr => 1,
	    OutputFile => $basePath . "/".$reportDir."/" . $data->{id} . ".xml",
	); 

	return $reportNameFile;
}

1;
