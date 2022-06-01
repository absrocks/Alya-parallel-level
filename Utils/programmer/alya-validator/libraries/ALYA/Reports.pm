package ALYA::Reports;

#---------------------------------------------------------------------
#--Package that creates the html reports from the xml result----------
#---------------------------------------------------------------------

#usage modules
use FindBin;                 # locate this script
use lib "$FindBin::Bin/.";  # use the parent directory
use File::Copy;
use File::Path;
use File::Compare;
use File::Copy::Recursive qw(dircopy);
use XML::Simple;
use Data::Dumper;
use POSIX;
#xsl test
use XML::XSLT;

#global variables
my $basePath = "modules/";
my $reportDir = "report";

#global variables


#----------------------------------------------------------------------
#-------------------------------functions------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#--Function generateReportName
#--Description:
#----Genereates the report file name, that includes the current time, the module name and the test name.
#--Parameters:
#----moduleName: The module's name that it's tested
#----testName: The test name.
sub generateReportName {
	#local($moduleName, $testName) = @_;
	#get local time in format: MMM_DD_YYYY
	$now_string = strftime "%b_%d_%Y", localtime;

	#Compound the file name
	$fileName = "report_" . $now_string;
	#opendir RH, $basePath . $reportDir or die "Couldn't open the current directory: $!";
        #check if this filename exists
        #$i = 0;
	#while ($_ = readdir(RH)) {
	#	next if $_ eq "." or $_ eq ".." or -d $_;
	#	if ($_ =~ /$fileName/) {
	#		$i++;
	#	}
	#}
        #if filename exists add an aditional identifier to the name
	#if ($i > 0) {
	#	$fileName .= "_" . ($i + 1);
	#}
        #add the extension
	#$fileName .= ".txt";

	return $fileName;
}

#----------------------------------------------------------------------
#--Function doHTMLReport
#--Description:
#----Creates an html report from the testSuite results
#--Parameters:
#----folderName: Folder name where the testSuite results are alocated
#----xmlFile: Name of the main xml testSuite result file
sub doHTMLReport {
	local($folderName, $xmlFile) = @_;
	
	#init the xsl
	my $xsl = "libraries/xsl/testSuite.xsl";
	my $xslt = XML::XSLT->new ($xsl, warnings => 1);

	#process the main xml file.
	$xslt->transform ($folderName . "/" .$xmlFile);
	open REPORT, ">> $folderName/index.html" or die $!;
	print REPORT $xslt->toString;
	close(REPORT);
	$xslt->dispose();

	#iterate over the sub reports
	#xml data
	my $xml = new XML::Simple;
	my $data = $xml->XMLin($folderName . "/" .$xmlFile);

	$resultsReport = $data->{moduleResults}->{moduleResult};
	@resultsRep = $resultsReport;
	if (ref $resultsReport eq 'ARRAY') {
		@resultsRep = @{$resultsReport};
	}
	#iterate over each test
	foreach $resRep (@resultsRep) {
		$name = $resRep->{moduleName};
		$xsl = "libraries/xsl/moduleResultHTML.xsl";
		$xslt = XML::XSLT->new ($xsl, warnings => 1);

		#process the main xml file.
		$xslt->transform ($folderName . "/" . $name . "/" . $name . ".xml" );
		open REPORT, ">> $folderName/$name/$name.html" or die $!;
		print REPORT $xslt->toString;
		close(REPORT);
		$xslt->dispose();

		#iterate over the tests
		my $dataTests = $xml->XMLin($folderName . "/" . $name . "/" . $name . ".xml");
		#Number of tests results
		$resultsReportTest = $dataTests->{testResults}->{testResult};
		@resultsRepTest = $resultsReportTest;
		if (ref $resultsReportTest eq 'ARRAY') {
			@resultsRepTest = @{$resultsReportTest};
		}

		#iterate over each test
		foreach $resRepTest (@resultsRepTest) {
			$testName = $resRepTest->{testName};
			$xsl = "libraries/xsl/testResultHTML.xsl";
			$xslt = XML::XSLT->new ($xsl, warnings => 1);

			#process the main xml file.
			$xslt->transform ($folderName . "/" . $name . "/" . $testName . "/" . $testName . ".xml" );
			open REPORT, ">> $folderName/$name/$testName/$testName.html" or die $!;
			print REPORT $xslt->toString;
			close(REPORT);
			$xslt->dispose();				
		}
	}
	#copy the css styles and images
	mkdir($folderName . "/libs");
	dircopy("libraries/libs", $folderName . "/libs");

}


#----------------------------------------------------------------------
#--Function doResumeReport
#--Description:
#----Creates a basic summary of the testSuite results in text format
#--Parameters:
#----folderName: Folder name where the testSuite results are alocated
#----xmlFile: Name of the main xml testSuite result file
sub doResumeReport {
	local($folderName, $xmlFile) = @_;

	$outputReport = "\n---------------------------------------------------------------------\n";
	$outputReport = $outputReport . "-------------------TestSuite Result----------------------------------\n";
	$outputReport = $outputReport . "---------------------------------------------------------------------\n";
	
	#iterate over the sub reports
	#xml data
	my $xml = new XML::Simple;
	my $data = $xml->XMLin($folderName . "/" .$xmlFile);

	$resultsReport = $data->{moduleResults}->{moduleResult};
	@resultsRep = $resultsReport;
	if (ref $resultsReport eq 'ARRAY') {
		@resultsRep = @{$resultsReport};
	}
	my $testPassed = 1;
	#iterate over each test
	foreach $resRep (@resultsRep) {
		#check if there is an error
		if ($resRep->{error}) {
			$testPassed = 0;
			last;
		}
		if ($resRep->{passed} eq 'false') {
			$testPassed = 0;
			last;
		}					
	}
	if ($testPassed) {
		$outputReport = $outputReport . "RESULT:            ALL TESTS PASSED CORRECTLY\n";
	}
	else {
		$outputReport = $outputReport . "RESULT:            SOME TESTS HAVE NOT PASSED\n";
	}
	$outputReport = $outputReport . "---------------------------------------------------------------------\n";
	$outputReport = $outputReport . "---------------------------------------------------------------------\n";
	$outputReport = $outputReport . "---------------------------------------------------------------------\n";

	#copy the css styles and images
	return $outputReport;
	

}

#----------------------------------------------------------------------
#--Function doEmailReport
#--Description:
#----Creates a basic summary of the testSuite and testSuiteJob results in text format to send to the general email.
#--Parameters:
#----folderName: Folder name where the testSuite results are alocated
#----xmlFile: Name of the main xml testSuite result file
#----jobResultXmlFile: Path to the xml testSuiteJob result
#----sshReport: Data of the ssh configuration information, used to create a link to the remote report
sub doEmailReport {
	local($folderName, $xmlFile, $jobResultXmlFile, $sshReport) = @_;

	$outputReport = "------------testSuite Job Result--------------------------\n";
	if ($sshReport) {
		$host = $sshReport->{sshHost};
		$project = $sshReport->{sshProject};
		$outputReport = $outputReport . "-----REPORT URL: http://$host/$project/index.html \n\n";
	}
	#xml data
	my $xml = new XML::Simple;
	my $data = $xml->XMLin($jobResultXmlFile);

	$outputReport = $outputReport . "\n------------Compilation-----------------------------------\n";
	#Modules
	$outputReport = $outputReport . "---Modules\n";
	$modulesCompResult = $data->{moduleCompilation}->{moduleCompilationResult};
	@modCompResults = $modulesCompResult;
	if (ref $modulesCompResult eq 'ARRAY') {
		@modCompResults = @{$modulesCompResult};
	}	
	foreach $modCompResult (@modCompResults) {
		$passed = "FAIL";
		if ($modCompResult->{passed} eq "true") {
			$passed = "OK"
		}
		$outputReport = $outputReport . "------" . $modCompResult->{moduleName} . " -->" . $passed ."\n";
	}

	#Services
	$outputReport = $outputReport . "---Services\n";
	$servicesCompResult = $data->{serviceCompilation}->{serviceCompilationResult};
	@serCompResults = $servicesCompResult;
	if (ref $servicesCompResult eq 'ARRAY') {
		@serCompResults = @{$servicesCompResult};
	}	
	foreach $serCompResult (@serCompResults) {
		$passed = "FAIL";
		if ($serCompResult->{passed} eq "true") {
			$passed = "OK"
		}
		$outputReport = $outputReport . "------" . $serCompResult->{serviceName} . " -->" . $passed ."\n";
	}

	$outputReport = $outputReport . "\n------------Tests------------------------------------------\n";

	#iterate over the sub reports
	$testResults = "";
	#xml data
	$data = $xml->XMLin($folderName . "/" .$xmlFile);

	$resultsReport = $data->{moduleResults}->{moduleResult};
	@resultsRep = $resultsReport;
	if (ref $resultsReport eq 'ARRAY') {
		@resultsRep = @{$resultsReport};
	}
	my $testPassed = 1;
	$subTestResult = "";
	#iterate over each test
	foreach $resRep (@resultsRep) {
		#check if there is an error
		if ($resRep->{error}) {
			$testPassed = 0;
		}
		if ($resRep->{passed} eq 'false') {
			$testPassed = 0;
		}
		#------------------------------
		#iterate over each subtest		
		$subData = $xml->XMLin($folderName . "/" . $resRep->{moduleName} . "/" . $resRep->{moduleName} . ".xml" );		
		#Number of tests results
		$resultsTest = $subData->{testResults}->{testResult};
		@resultTest = $resultsTest;
		if (ref $resultsTest eq 'ARRAY') {
			@resultTest = @{$resultsTest};
		}

		#iterate over each test
		foreach $resTest (@resultTest) {
			$subPassed = "OK";
			#check if there is an error
			if ($resTest->{error}) {
				$subPassed = "FAIL";
			}
			if ($resTest->{passed} eq 'false') {
				$subPassed = "FAIL";
			}
			$subTestResult = $subTestResult . "-------" . $resTest->{testName} . " -->" . $subPassed . "\n";
		}		
	}
	$testResult = "";
	if ($testPassed) {
		$testResult = "------------ALL TESTS PASSED SUCCESSFULLY------------------\n";
	}
	else {
		$testResult = "------------SOME TESTS HAVE NOT PASSED------------------\n";
	}

	$outputReport = $outputReport . $testResult;
	$outputReport = $outputReport . "-----------------------------------------------------------\n\n";
	$outputReport = $outputReport . $subTestResult;

	#copy the css styles and images
	return $outputReport;
	

}

#----------------------------------------------------------------------
#--Function doEmailReportManager
#--Description:
#----Creates a basic summary of the testSuite and testSuiteJob results in text format to send to the managers
#----The report is personalized for a concrete manager, it only indicates manager modules that fails
#--Parameters:
#----folderName: Folder name where the testSuite results are alocated
#----xmlFile: Name of the main xml testSuite result file
#----jobResultXmlFile: Path to the xml testSuiteJob result
#----manager: Information about the modules managed by the manager
#----allManagerTests: Contains all the tests that are managed by any manager
sub doEmailReportManager {
	local($folderName, $xmlFile, $jobResultXmlFile , $manager, $allManagerTests) = @_;
	$withCompileResults = 0;	

	$outputReport = "------------testSuite Job Result--------------------------\n";
	#xml data
	my $xml = new XML::Simple;
	my $data = $xml->XMLin($jobResultXmlFile);

	$outputReport = $outputReport . "\n------------Compilation-----------------------------------\n";

     
	#---------------------------Modules-----------------------------------

	#put the modules manager into hash managerModules
	$managerModulesList = $manager->{managerModule};
	@managerModuleList = $managerModulesList;
	if (ref $managerModulesList eq 'ARRAY') {
		@managerModuleList = @{$managerModulesList};
	}
	my %managerModules;
	@managerModules{@managerModuleList}=();


        #Report each module result
	$outputReport = $outputReport . "---Modules\n";
	$modulesCompResult = $data->{moduleCompilation}->{moduleCompilationResult};
	@modCompResults = $modulesCompResult;
	if (ref $modulesCompResult eq 'ARRAY') {
		@modCompResults = @{$modulesCompResult};
	}	
	foreach $modCompResult (@modCompResults) {
		$passed = "FAIL";
		if ($modCompResult->{passed} eq "true") {
			$passed = "OK";
		}
                #if module is managed by the manager and it has errors then report the response
		if (exists $managerModules{$modCompResult->{moduleName}} && $passed eq "FAIL") {
			$outputReport = $outputReport . "------" . $modCompResult->{moduleName} . " -->" . $passed ."\n";
			$withCompileResults = 1;
		}
	}

	#--------------------------Services--------------------------------

	#put the services manager into hash managerModules
	$managerServicesList = $manager->{managerService};
	@managerServiceList = $managerServicesList;
	if (ref $managerServicesList eq 'ARRAY') {
		@managerServiceList = @{$managerServicesList};
	}
	my %managerServices;
	@managerServices{@managerServiceList}=();


	#report each service result
	$outputReport = $outputReport . "---Services\n";
	$servicesCompResult = $data->{serviceCompilation}->{serviceCompilationResult};
	@serCompResults = $servicesCompResult;
	if (ref $servicesCompResult eq 'ARRAY') {
		@serCompResults = @{$servicesCompResult};
	}	
	foreach $serCompResult (@serCompResults) {
		$passed = "FAIL";
		if ($serCompResult->{passed} eq "true") {
			$passed = "OK";
		}
		#if module is managed by the manager and it has errors then report the response
		if (exists $managerServices{$serCompResult->{serviceName}} && $passed eq "FAIL") {
			$outputReport = $outputReport . "------" . $serCompResult->{serviceName} . " -->" . $passed ."\n";
			$withCompileResults = 1;
		}
	}

	if (!$withCompileResults) {
		$outputReport = "------------testSuite Job Result--------------------------\n";		
	}

	#-------------------------------------------tests-----------------------------------------
	$withTestResults = 0;
	#put the tests manager into hash managerModules
	$managerTestsList = $manager->{managerTest};
	@managerTestList = $managerTestsList;
	if (ref $managerTestsList eq 'ARRAY') {
		@managerTestList = @{$managerTestsList};
	}
	my %managerTests;
	@managerTests{@managerTestList}=();

	$outputReport = $outputReport . "\n------------Tests------------------------------------------\n";

	#iterate over the sub reports
	$testResults = "";
	#xml data
	$data = $xml->XMLin($folderName . "/" .$xmlFile);

	$resultsReport = $data->{moduleResults}->{moduleResult};
	@resultsRep = $resultsReport;
	if (ref $resultsReport eq 'ARRAY') {
		@resultsRep = @{$resultsReport};
	}
	my $testPassed = 1;
	$subTestResult = "";
	#iterate over each test
	foreach $resRep (@resultsRep) {
		#check if there is an error
		if ($resRep->{error}) {
			$testPassed = 0;
		}
		if ($resRep->{passed} eq 'false') {
			$testPassed = 0;
		}
		#------------------------------
		#iterate over each subtest		
		$subData = $xml->XMLin($folderName . "/" . $resRep->{moduleName} . "/" . $resRep->{moduleName} . ".xml" );		
		#Number of tests results
		$resultsTest = $subData->{testResults}->{testResult};
		@resultTest = $resultsTest;
		if (ref $resultsTest eq 'ARRAY') {
			@resultTest = @{$resultsTest};
		}

		#iterate over each test
		foreach $resTest (@resultTest) {
			$subPassed = "OK";
			#check if there is an error
			if ($resTest->{error}) {
				$subPassed = "FAIL";
			}
			if ($resTest->{passed} eq 'false') {
				$subPassed = "FAIL";
			}
			#if module is managed by the manager and it has errors then report the response
			if (exists $managerTests{$resTest->{testName}} && $subPassed eq "FAIL") {
				$subTestResult = $subTestResult . "-------" . $resTest->{testName} . " -->" . $subPassed . "\n";
				$withTestResults = 1;
			}
			elsif (exists $managerModules{$resRep->{moduleName}} && $subPassed eq "FAIL") {
				print "resTest:" . $resTest->{testName} . "\n";
				if (!exists $allManagerTests->{$resTest->{testName}}) {
					$subTestResult = $subTestResult . "-------" . $resTest->{testName} . " -->" . $subPassed . "\n";
					$withTestResults = 1;	
				}				
			}
		}		
	}
	$testResult = "";
	if ($testPassed) {
		$testResult = "------------ALL TESTS PASSED SUCCESSFULLY------------------\n";
	}
	else {
		$testResult = "------------SOME TESTS HAVE NOT PASSED------------------\n";
	}

	if ($withTestResults) {
		$outputReport = $outputReport . $testResult;
		$outputReport = $outputReport . "-----------------------------------------------------------\n\n";
		$outputReport = $outputReport . $subTestResult;
	}
	

	if ($withCompileResults == 0 && $withTestResults == 0) {
		$outputReport = 0;
	}

	#copy the css styles and images
	return $outputReport;
	

}

#----------------------------------------------------------------------
#--Function doGlobalReport
#--Description:
#----Joins the testSuite html report with testSuiteJob xml result to create an unified html report
#--Parameters:
#----folderName: Folder name where the testSuite results are alocated
#----xmlFile: Name of the main xml testSuite result file
#----repDir: Report folder name
#----jobResultXmlFile: Path to the xml testSuiteJob result
#----sshReport: Data of the ssh configuration information, used to create a link to the remote report
sub doGlobalReport {
	local($xmlFile, $repDir, $jobResultXmlFile, $sshReport, $svnRevision, $metaInfo) = @_;

	#Join the two report xml's
	open FILE, ">$repDir/globalResult.xml" or die $!;

	print FILE "<testSuite>\n";
	#information section;
	print FILE "<information>\n";
	$platform = $metaInfo->{platform};
	print FILE "    <platform>$platform</platform>\n";
	print FILE "    <svnRevision>$svnRevision</svnRevision>\n";
	$compiler = $metaInfo->{compiler};
	print FILE "    <compiler>$compiler</compiler>\n";
	$now_string = strftime "%b_%d_%Y", localtime;
	print FILE "    <date>$now_string</date>\n";
	print FILE "</information>\n";
	open JOBRESULT, "$repDir/$jobResultXmlFile" or die "Couldn't open file: $!"; 
	while (<JOBRESULT>){
		print FILE $_;
	}
	close JOBRESULT;
	open JOBRESULT, "$repDir/$xmlFile"; 
	while (<JOBRESULT>){
		print FILE $_;
	}
	close JOBRESULT;
	print FILE "</testSuite>";
	close FILE;

	#exec the xsl
	#init the xsl
	my $xsl = "libraries/xsl/testSuiteGlobal.xsl";
	my $xslt = XML::XSLT->new ($xsl, warnings => 1);

	#process the main xml file.
	$xslt->transform ("$repDir/globalResult.xml");
	open REPORT, ">$repDir/index.html" or die $!;
	print REPORT $xslt->toString;
	close(REPORT);
	$xslt->dispose();
}


1;
