package ALYA::LSFExec;

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
	my $jobSubmit;
	my $execOutput;
	eval {
		local $SIG{ALRM} = sub {die "alarm\n"};
		alarm $timeout;
		createJobScript($basePath, $origPath, $numProcessor, $alyaExec, $testId, $timeout, 1);
                #submit the job
		$jobSubmit = `bsub < job.cmd`;
		#get the job id
		if ($jobSubmit =~ m/Job <(\d+)> is submitted/)
      		{
		      select(undef, undef, undef, 6);
		      $jobSubmit = $1;
	        }
		#waits until job finished
		$jobStatus = `bjobs -X`;
		while ($jobStatus =~ m/$jobSubmit/) {
			select(undef, undef, undef, 10);
			$jobStatus = `bjobs -X`;
		}
		alarm 0;
	};

	if ($@) {
		die unless $@ eq "alarm\n";
		#kill the process
		$jobCancel = `bkill $jobSubmit`;
		$execOutput = "timed out, it takes more than " . $timeout . " seconds\n" . $execOutput;
	}
	else {
		print "didn't time out\n";
	}
	#Join all the output in one string, last 30 lines from each file
	my $output =  joinTestOutput();
	$execOutput = $execOutput . "\n" . $output;
	#remove files
	unlink("job.cmd");
	unlink("alya.out");
	unlink("alya.err");
	unlink("out.txt");
	chdir $origPath;

	return $execOutput;
}


sub joinTestOutput {

	#Error file
	open FILE, "<alya.err";
	my $file_contents = do { local $/; <FILE> };
	close FILE;

	my @lines = split(/\n/, $file_contents);
	my $numLines = scalar (@lines);
	$finalLines = $file_contents;
        if ($numLines > 30) {
		$finalLines = join("\n", @lines[($numLines-30) .. ($numLines-1)]);
	}

	my $finalOutput = "-------------ERRORS:\n";
	$finalOutput = $finalOutput . $finalLines . "\n";

	#Job output
	open FILE, "<alya.out";
	$file_contents = do { local $/; <FILE> };
	close FILE;

	@lines = split(/\n/, $file_contents);
	$numLines = scalar (@lines);
	$finalLines = $file_contents;
        if ($numLines > 30) {
		$finalLines = join("\n", @lines[($numLines-30) .. ($numLines-1)]);
	}

	$finalOutput = $finalOutput . "-------------JOB OUTPUT:\n";
	$finalOutput = $finalOutput . $finalLines . "\n";

	#Stdout
	open FILE, "<out.txt";
	$file_contents = do { local $/; <FILE> };
	close FILE;

	@lines = split(/\n/, $file_contents);
	$numLines = scalar (@lines);
	$finalLines = $file_contents;
        if ($numLines > 30) {
		$finalLines = join("\n", @lines[($numLines-30) .. ($numLines-1)]);
	}

	$finalOutput = $finalOutput . "-------------STDOUT:\n";
	$finalOutput = $finalOutput . $finalLines . "\n";
	
	return $finalOutput;
}

sub createJobScript {

   local($basePath, $origPath, $numProcessor, $alyaExec, $testId, $timeout, $withMPI) = @_;

   my $time = convert_time($timeout);
   my $testName = $testId;
   if ($testName eq "") {
	$testName = "alyaCommand"; 
   }

   open FILE, ">job.cmd" or die $!;
   print FILE "#! /bin/bash\n";
   print FILE "\n";

   print FILE "#BSUB -n $numProcessor\n";
   print FILE "#BSUB -R\"span[ptile=2]\"\n";
   print FILE "#BSUB -q debug\n";
   print FILE "#BSUB -o alya.out\n";
   print FILE "#BSUB -e alya.err\n";
   print FILE "#BSUB -J $testName\n";
   print FILE "#BSUB -W $time\n";
   #print FILE "module switch openmpi impi\n";
   print FILE "echo '--| ALYA: JOB SUBMITTED WITH SRUN AT: ' `date`\n";
   if ($withMPI) {
	print FILE "mpirun $alyaExec $testId>out.txt\n";
   }
   else {
	print FILE "$alyaExec $testId>out.txt\n";
   }   
   print FILE "echo '--| ALYA: JOB FINISHED AT: ' `date`\n";   
   close FILE;
   #darle permisos de ejecucion 
   $output = `chmod +x job.cmd`;

}

sub convert_time { 
  my $time = shift;  
  my $hours = int($time / 3600); 
  $time -= ($hours * 3600); 
  my $minutes = int($time / 60);
  
  if ($hours < 1) {
	$hours = '00'
  }
  elsif ($hours < 10) {
	$hours = '0' . $hours;
  }
  if ($minutes < 1) {
	$minutes = '00'
  }
  elsif ($minutes < 10) {
	$minutes = '0' . $minutes;
  }
 
  $time = $hours . ":" . $minutes; 
  return $time; 
}

sub execWithOpen {
	local($self, $basePath, $origPath, $command, $timeout) = @_;
	my $execOutput = "";
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
	my $execOutput = "";
	if (!$timeout) {
		$timeout = 0;
	}
	$pid = "";
	eval {
		local $SIG{ALRM} = sub {die "alarm\n"};
		alarm $timeout;
		createJobScript("", "", 1, $command, "", $timeout, 0);
                #submit the job
		$jobSubmit = `bsub < job.cmd`;
		#get the job id
		if ($jobSubmit =~ m/Job <(\d+)> is submitted/)
      		{
		      $jobSubmit = $1;
		      select(undef, undef, undef, 6);
	        }
		#waits until job finished
		$jobStatus = `bjobs -X`;
		while ($jobStatus =~ m/$jobSubmit/) {
			select(undef, undef, undef, 10);
			$jobStatus = `bjobs -X`;
		}
		alarm 0;		
	};
	
	if ($@) {
		die unless $@ eq "alarm\n";
		#kill the process
		$jobCancel = `mncancel $jobSubmit`;
		$execOutput = "timed out, it takes more than " . $timeout . " seconds\n" . $execOutput;
	}
	else {
		print "didn't time out\n";
	}

	#Join all the output in one string, last 30 lines from each file
	$output =  joinTestOutput();
	$execOutput = $execOutput . "\n" . $output;
	#remove files
	unlink("job.sh");
	unlink("alya.out");
	unlink("alya.err");
	unlink("out.txt");
	return $execOutput;
}

sub execServerJob {
	local($self, $user, $host, $bashJob) = @_;

	print "Launching job on remote server LSF queue\n";
	my $execOutput = `ssh $user\@$host "source /etc/profile; bsub < $bashJob"`;

	print $execOutput . "\n";

	#2.Obtain the job id
	#Submitted batch job 381438
	my $jobId;
	#Job <34288> is submitted to default queue <bsc_case>
	if($execOutput =~ m/Job <(\d+)> is submitted/) {
		$jobId = $1;
		print "jobId:" . $jobId . "\n";
		sleep(6);
	}

	#3.monitorize
	my $jobFinished = 0;
	while(!$jobFinished) {
		$execOutput = `ssh $user\@$host bjobs`;		
		if($execOutput =~ m/$jobId/) {
			print "monitoring:" . $execOutput . "\n";
			sleep(600);
		}
		else {
			$jobFinished = 1;	
		}
	}
	print "job Finished\n";
}


1;
