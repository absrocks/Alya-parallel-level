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


#obtain the last revision
open FILE, "lastRevision.txt" or die "Couldn't open file: $!";
binmode FILE;
my $from = <FILE>;
close FILE;
$from =~ s/\R//g;

#obtain the actual revision
$svnInfo = `svn info $svnURL --xml`;
if ($svnInfo =~ m/revision=\"(\d+)\"\>/) {
	my $to = $1;
	if ($to != $from) {
		while ($from < $to) {
			$from = $from + 1;
			#obtain revision info xml
			$svnLog = `svn log $svnURL --xml --verbose -r $from`;
			#process the xml
			my $xsl = "libraries/xsl/svnLog.xsl";
			my $xslt = XML::XSLT->new ($xsl, warnings => 1);
			$xslt->transform ($svnLog);	
			my $result = $xslt->toString;
			$xslt->dispose();
			#get the user
			my $svnUser;
			if ($svnLog =~ m/\<author\>(\w+)\<\/author\>/) {
				$svnUser = $1;
			}
			ALYA::Utils::sendEmail($data->{mail}, "[AlyaSVNCommit]$from-$svnUser", 0, $result, 0);			
		}
		#update the last revision
		open FILE, ">lastRevision.txt" or die $!;
		print FILE $to;
		close FILE;
	}		
}

exit 0
