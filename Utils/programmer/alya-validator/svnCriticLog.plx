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
use XML::XSLT;


#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#--------------------------------MAIN-----------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
$xml = new XML::Simple;
$data = $xml->XMLin("svnCommitsJob.xml");
my $host = $data->{host};
my $user = $data->{user};
$svnURL = $data->{svnURL};
$password = $data->{password};
my $sshConf = $data->{sshReport};


$svnLog = `svn log $svnURL --xml --verbose -r 2170`;
#process the xml
my $xsl = "libraries/xsl/svnCritical.xsl";
my $xslt = XML::XSLT->new ($xsl, warnings => 1);
$xslt->transform ($svnLog);	
my $result = $xslt->toString;
$xslt->dispose();

my $response = "------------------CRITICAL COMMITS SUMMARY---------------------\n";
$response = $response . "---------------------------------------------------------------\n";
my @commits = split /----LOGENTRY/, $result;
foreach my $commit (@commits) {

	if ($commit =~ m/----MESSAGE:TYPE=CRITICAL/) {
		$response = $response . "--------------------------------\n";
		if ($commit =~ m/----REVISION:(.+)/) {
			$response = $response . "----REVISION:$1\n";
		}
		if ($commit =~ m/----AUTHOR:(.+)/) {
			$response = $response . "----AUTHOR:$1\n";
		}
		if ($commit =~ m/----CHANGES:(.+)/) {			
			my $changes = $1;
			my $modules = {};
			my @parts = split /#/, $changes;
			foreach my $part (@parts) {
				if ($part =~ m/\/modules\/(.+?)\//) {
					$modules->{$1} = $1;
				}
				if ($part =~ m/\/services\/(.+?)\//) {
					$modules->{$1} = $1;
				}
				if ($part =~ m/\/kernel\//) {
					$modules->{kernel} = "kernel";
				}
			}
			my @keys = keys %{$modules};
			my $modAux;
			my $total = scalar(@keys);
			my $index = 1;
			foreach my $key (@keys) {
				$modAux = $modAux . $key;
				if ($index != $total)  {
					$modAux = $modAux . ",";
				}
				$index = $index + 1;
			}
			$response = $response . "----CHANGES:$modAux\n";
		}
		if ($commit =~ m/----MESSAGE:([.+|\n])/) {
			$response = $response . "----MESSAGE:$1\n";
		}
	}
}

print $response;


exit 0
