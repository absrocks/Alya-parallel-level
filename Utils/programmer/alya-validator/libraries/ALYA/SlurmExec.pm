package ALYA::SlurmExec;

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
	$jobSubmit;
	$execOutput;
	eval {
		local $SIG{ALRM} = sub {die "alarm\n"};
		alarm $timeout;
		createJobScript($basePath, $origPath, $numProcessor, $alyaExec, $testId, $timeout);
                #submit the job
		$jobSubmit = `mnsubmit job.sh`;
		#get the job id
		if ($jobSubmit =~ m/Submitted batch job (\d+)/)
      		{
		      $jobSubmit = $1;
	        }
		#waits until job finished
		$jobStatus = `mnq -j $jobSubmit`;
		while ($jobStatus =~ m/$jobSubmit/) {
			select(undef, undef, undef, 10);
			$jobStatus = `mnq -j $jobSubmit`;
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

   local($basePath, $origPath, $numProcessor, $alyaExec, $testId, $timeout) = @_;

   my $time = convert_time($timeout);
   my $testName = $testId;
   if ($testName eq "") {
	$testName = "alyaCommand"; 
   }

   open FILE, ">job.sh" or die $!;
   print FILE "#! /bin/bash\n";
   print FILE "#\n";
   print FILE "#@ job_name         = $testName\n";
   print FILE "#@ output           = alya.out\n";
   print FILE "#@ error            = alya.err\n";
   print FILE "#@ initialdir       = .\n";
   print FILE "#@ total_tasks      = $numProcessor\n";
   print FILE "#@ tasks_per_node   = 1\n";
   #print FILE "#@ class = debug\n";
   print FILE "#@ wall_clock_limit =$time\n";
   print FILE "#@ queue\n";
   print FILE "\n";
   print FILE "echo '--| ALYA: JOB SUBMITTED WITH SRUN AT: ' `date`\n";
   print FILE "time srun $alyaExec $testId>out.txt\n";
   print FILE "echo '--| ALYA: JOB FINISHED AT: ' `date`\n";
   close FILE;
   #darle permisos de ejecucion 
   $output = `chmod +x job.sh`;

}

sub convert_time { 
  my $time = shift;  
  my $hours = int($time / 3600); 
  $time -= ($hours * 3600); 
  my $minutes = int($time / 60); 
  my $seconds = $time % 60;
  
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
  if ($seconds < 1) {
	$seconds = '00'
  }
  elsif ($seconds < 10) {
	$seconds = '0' . $seconds;
  }

  $hours = $hours < 1 ? '00' : $hours .'h '; 
  $time = $hours . ":" . $minutes . ":" . $seconds; 
  return $time; 
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
		createJobScript("", "", 1, $command, "", $timeout);
                #submit the job
		$jobSubmit = `mnsubmit job.sh`;
		#get the job id
		if ($jobSubmit =~ m/Submitted batch job (\d+)/)
      		{
		      $jobSubmit = $1;
	        }
		#waits until job finished
		$jobStatus = `mnq -j $jobSubmit`;
		while ($jobStatus =~ m/$jobSubmit/) {
			select(undef, undef, undef, 10);
			$jobStatus = `mnq -j $jobSubmit`;
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
	print "User:" . $user . " host:" . $host . " bashJob:" . $bashJob . "\n";
	print "Launching job on remote server slurm queue\n";
	$execOutput = `ssh $user\@$host "source /etc/profile; mnsubmit ./$bashJob"`;

	print $execOutput . "\n";

	#2.Obtain the job id
	#Submitted batch job 381438
	my $jobId;
	if($execOutput =~ m/Submitted batch job (\d+)/) {
		$jobId = $1;
		print "jobId:" . $jobId . "\n";
	}


	#3.monitorize
	my $jobFinished = 0;
	while(!$jobFinished) {
		$execOutput = `ssh $user\@$host mnq -j $jobId`;
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
