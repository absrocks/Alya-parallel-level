package ALYA::LocalExec;

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

#global variables
$execOutput = "";

sub new { 
  my $class = shift; 
  my $self = { }; 
  bless $self, $class;
  return $self;
}

#----------------------------------------------------------------------
#-------------------------------functions------------------------------
#----------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#--------------------------------MAIN-----------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
sub execTest {
	local($self, $basePath, $origPath, $numProcessor, $alyaExec, $testId, $timeout) = @_;
	
	if (!$timeout) {
		$timeout = 0;
	}
	chdir $basePath . "/tmp";
	print "Start execution with $numProcessor processors\n";
	eval {
		local $SIG{ALRM} = sub {die "alarm\n"};
		alarm $timeout;
		$execOutput =  `mpirun -n $numProcessor $alyaExec $testId`;
		alarm 0;
	};

	if ($@) {
		die unless $@ eq "alarm\n";
		#kill the process
		system 'killall', 'mpirun';
		$execOutput = "timed out, it takes more than " . $timeout . " seconds\n" . $execOutput;
	}
	else {
		print "didn't time out\n";
	}
	chdir $origPath;
	return $execOutput;
}

sub execWithOpen {
	local($self, $basePath, $origPath, $command, $timeout) = @_;
	$execOutput = "";
	if (!$timeout) {
		$timeout = 0;
	}
	chdir $basePath . "/tmp";
	eval {
		local $SIG{ALRM} = sub {die "alarm\n"};
		alarm $timeout;
		$pid = open(WRITEME, "| $command");
		print WRITEME "\n";
		close(WRITEME);
		alarm 0;
	};

	if ($@) {
		die unless $@ eq "alarm\n";
		#kill the process
		system 'kill', '-9' , $pid;
		$execOutput = "timed out, it takes more than " . $timeout . " seconds\n";
	}
	else {
		print "didn't time out\n";
	}
	chdir $origPath;

	return $execOutput;
}

sub execCommand {
	local($self, $command, $timeout) = @_;
	print "Executing command:" . $command . "\n";
	$execOutput = "";
	if (!$timeout) {
		$timeout = 0;
	}
	$pid = "";
	eval {
		local $SIG{ALRM} = sub {die "alarm\n"};
		alarm $timeout;
		#$execOutput = `$command 2>error.out`;
		open FILE, ">>output.out" or die $!;
		$pid = open(COMM, "$command 2>&1 |");
		while ( defined(my $line = <COMM> ) )
		{
			chomp($line);
			print FILE $line . "\n";
		}
		close(COMM);
		close FILE;
		alarm 0;
		print $execOutput;
	};
	
	if ($@) {
		die unless $@ eq "alarm\n";		
		#kill the process
		system 'kill', '-9', $pid;
		close(COMM);
		close FILE;
		$execOutput = "timed out, it takes more than " . $timeout . " seconds\n";
	}
	else {
		print "didn't time out\n";
	}

	#loads the output into the var
	open FILE, "<output.out";
	$file_contents = do { local $/; <FILE> };
	close FILE;
	$execOutput = $execOutput . $file_contents;
	unlink "output.out";
	return $execOutput;
}


1;
