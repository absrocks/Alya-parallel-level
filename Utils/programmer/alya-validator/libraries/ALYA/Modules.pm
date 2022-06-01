package ALYA::Modules;

#---------------------------------------------------------------------
#--Package that executes all tests associated to one module-----------
#---------------------------------------------------------------------

#usage modules
use FindBin;                 # locate this script
use lib "$FindBin::Bin/.";  # use the parent directory
use File::Copy;
use File::Path;
use File::Compare;
use File::Copy::Recursive qw(dirmove);
use XML::Simple;
use ALYA::Tests;
use Data::Dumper;
use POSIX;

#global variables
my $basePath = "modules/";
my $reportDir = "report";

#global variables


#----------------------------------------------------------------------
#-------------------------------functions------------------------------
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#--Function checkIfTestPassed
#--Description:
#----Iterates over all the test results and determines if the test has passed
#--Parameters:
#----testId: The name of the test to check
#----result: The xml report variable where to write if test has passed or not
sub checkIfTestPassed {

local($testId, $result) = @_;

	my $xml = new XML::Simple;
	#load the report xml generated, and check if the module passes the test
	$testReport = $xml->XMLin($basePath . "/" . $testId . "/" . $reportDir . "/" . $testId . ".xml");
	my $testPassed = 1;
	
	#Number of tests results
	$resultsReport = $testReport->{results}->{result};
	@resultsRep = $resultsReport;
	if (ref $resultsReport eq 'ARRAY') {
		@resultsRep = @{$resultsReport};
	}

	#iterate over each test
	OUTER_LOOP:
	foreach $resRep (@resultsRep) {
		#check if there is an error
		if ($repRep->{error}) {
			$testPassed = 0;
			last;
		}
		else {
			#iterate over each file			
			$filesReport = $resRep->{file};
			@filesRep = $filesReport;
			if (ref $filesReport eq 'ARRAY') {
				@filesRep = @{$filesReport};
			}
	
			#iterate over each test
			foreach $fileRep (@filesRep) {
				if ($fileRep->{passed} eq 'false') {
					$testPassed = 0;
					last OUTER_LOOP;
				}
				if ($fileRep->{error}) {
					$testPassed = 0;
					last OUTER_LOOP;
				}	
			}
		}
			
	}
	if ($testPassed) {
		$result->{passed} = "true";
	}
	else {
		$result->{passed} = "false";
	}
}

#----------------------------------------------------------------------
#--Function checkIfTestPassed
#--Description:
#----Iterates over all the test results and determines if the test has passed
#--Parameters:
#----testId: The name of the test to check
#----result: The xml report variable where to write if test has passed or not
sub moveReport {

local($testId) = @_;
	rmtree( $basePath . "/" . $reportDir . "/" . $testId);
        mkdir( $basePath . "/" . $reportDir. "/" . $testId);
	dirmove($basePath . "/" . $testId . "/" . $reportDir , $basePath . "/" . $reportDir. "/" . $testId) or die("$!\n");

	#opendir DH, $basePath . "/" . $testId . "/" . $reportDir or die "Couldn't open the current directory: $!";
	#while ($_ = readdir(DH)) {
	#	next if $_ eq "." or $_ eq ".." or -d $_;
	#	copy($basePath . "/" . $testId . "/" . $reportDir . "/" . $_, 
	#	     $basePath . "/" . $reportDir. "/" . $testId . "/" . $_)
	#		    or die "copy failed: $!";		
	#}
	#rmtree( $basePath . "/" . $testId . "/" . $reportDir);
}

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#--------------------------------MAIN-----------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
sub doModules {
	local($common, $moduleName) = @_;
        $basePath = "modules/";

	#xml data
	$xml = new XML::Simple;
	$data = $xml->XMLin($basePath . $moduleName . ".xml");

	$basePath = $basePath . $moduleName;

        #report xml initialization
        rmtree( $basePath . "/" . $reportDir );
        mkdir( $basePath . "/" . $reportDir );
	
	my $report = {};
	$report->{moduleResult} = {};
        $report->{moduleResult}->{name} = $moduleName;
	$report->{moduleResult}->{testResults}->{testResult} = [];
        my $results = $report->{moduleResult}->{testResults};

	#Number of tests
	$modulesTest = $data->{tests}->{test};
	@tests = $modulesTest;
	if (ref $modulesTest eq 'ARRAY') {
		@tests = @{$modulesTest};
	}
	
	#iterate over each test
	foreach $test (@tests) {

		#log
		$result = {};
                $result->{testName} = $test->{testName};
		$result->{description} = $test->{description};
                push(@{$results->{testResult}}, $result);

		#1.Exec the test
		print "Start execution of test $test->{testName}\n";
		ALYA::Tests::doTests($common, $moduleName , $test->{testName});
		
		#2.Checks if test passed
		checkIfTestPassed($test->{testName} , $result);
		
		#3.move the report test to the module
		moveReport($test->{testName});
	}

	#write the report
	$xml->XMLout(
	    $report,
	    KeepRoot => 1,
	    NoAttr => 1,
	    OutputFile => $basePath . "/" . $reportDir . "/" . $moduleName . ".xml",
	);
}

1;
