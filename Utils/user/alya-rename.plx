#!/usr/bin/perl
# configuration
use FindBin;                # locate this script
use lib "$FindBin::Bin/libraries";  # use the parent directory
use File::Copy;
use File::Path;
use File::Compare;
use File::Copy qw(move);
use File::Find;
use Data::Dumper;
use POSIX;


#------------------------------------------------
#---------------main-----------------------------
#------------------------------------------------

#-----Get the arguments------
my $numArgs = $#ARGV + 1;
my $basePath = ".";
my $problemName;
my $newProblemName;
if ($numArgs < 2) {
	print "The usage is: alya-rename.plx PROBLEM_NAME NEW_PROBLEM_NAME [PROBLEM_PATH]";
	exit(0);
}
else {
	$problemName = $ARGV[0];
	$newProblemName = $ARGV[1];
	if ($numArgs == 3) {
		$basePath = $ARGV[2];
	}
}

#openDir
opendir DH, $basePath or die "Couldn't open the current directory: $!";
while ($_ = readdir(DH)) {
	next if $_ eq "." or $_ eq ".." or -d ($basePath . "/" . $_);
	my $fileName = $_;

	#--------------modify the lines where the fileName appears----------------
	my $totalPath = $basePath . "/" . $fileName;
	#Open the file and read data
	#Die with grace if it fails
	open (FILE, "<$totalPath") or die "Can't open $totalPath: $!\n";
	my @lines = <FILE>;
	close FILE;
	#Open same file for writing, reusing STDOUT
	open (FILE, ">$totalPath") or die "Can't open $totalPath: $!\n";
	#Walk through lines, putting into $_, and substitute 2nd away
	for ( @lines ) {
		s/$problemName\./$newProblemName\./;
		print FILE;
	}
	#Finish up
	close FILE;
	#--------------modify file name----------------
	my $newFileName = $fileName;
	$newFileName =~ s/^$problemName\./$newProblemName\./;
	if ($fileName ne $newFileName) {
		move($basePath . "/" . $fileName , $basePath . "/" . $newFileName) or die "copy failed: $!";
	}
}
exit(0);
