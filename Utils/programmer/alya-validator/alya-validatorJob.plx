#!/usr/bin/perl

#usage modules
use FindBin;                 # locate this script
use lib "$FindBin::Bin/libraries";  # use the parent directory
use File::Copy;
use File::Path;
use File::Find;
use File::Compare;
use File::Copy::Recursive qw(dirmove);
use XML::Simple;
use Data::Dumper;
use POSIX;
use ALYA::Reports;
use ALYA::Utils;
use Archive::Tar;
use ALYA::ExecFactory;


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
	#get local time in format: MMM_DD_YYYY
	$now_string = strftime "%b_%d_%Y", localtime;

	#Compound the file name
	$fileName = "report_" . $now_string;

	return $fileName;
}


#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#--------------------------------MAIN-----------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
$xml = new XML::Simple;
$data = $xml->XMLin("alya-validatorJob.xml");
my $host = $data->{host};
my $user = $data->{user};
$svnURL = $data->{svnURL};
$password = $data->{password};
my $sshConf = $data->{sshReport};

#1.checkout from svn
print "Checkout ALYA from repository\n";
$reportFolder = "Report";
rmtree($reportFolder);
mkdir($reportFolder);
$output = `svn co $svnURL --password $password $reportFolder`;
print "Start server execution\n";
chdir "$reportFolder/Utils/programmer/alya-validator";
unlink "output.out";
$execOutput = `./alya-validator.plx > output.out`;

my $svnName = generateReportName();
if (-f "output.out") {
        print "Start untar report file\n";
        #$result = `tar -cvzf report.tar.gz output.out`;
	#9.4.General email
	my $message = "------------testSuite Job Result--------------------------\n";
	$host = $sshConf->{sshHost};
	$message = $message . "-----REPORT URL: http://$host/private/output.out \n\n";	
	ALYA::Utils::sendEmail($data->{mail}, "[AlyaValidator]$svnName", "output.out", $message, 0);

} else {
	print "File output.out not found\n";
}

#9.6.sshcopy
$output = `scp -r output.out $sshConf->{sshUser}\@$sshConf->{sshHost}:$sshConf->{sshFolder}`;
chdir "../../../..";
rmtree($reportFolder);

exit 0 
