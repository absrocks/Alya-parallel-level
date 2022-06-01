package ALYA::Utils;

#--------------------------------------------------------------------------------
#--Package that contains common functions used in other packages result----------
#--------------------------------------------------------------------------------

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
use Net::SMTP::SSL;
use MIME::Lite;
use ALYA::ExecFactory;
use ALYA::Reports;


#----------------------------------------------------------------------
#-------------------------------functions------------------------------
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#--Function sendEmail
#--Description:
#----sends the general report email to the people specified in the email xml configuration file
#--Parameters:
#----mailConf: Pointer to a var that contins the email configuration, the from and the to
#----subject: The email subject
#----svnName: The name of the report to send, with this var the function tries to send a .tar.gz file with the report
#----message: the email message
#----sendOutput: Boolean var, if true the testSuiteJob output is send as an attached file
sub sendEmail {
	local($mailConf, $subject, $fileName, $message, $sendOutput) = @_;

	$path = $fileName;

	#Services
	$mailsTo = $mailConf->{mailTo};
	@emailsTo = $mailsTo;
	if (ref $mailsTo eq 'ARRAY') {
		@emailsTo = @{$mailsTo};
	}	
	#creates the configure command
	foreach $mailTo (@emailsTo) {
		print "SENDING EMAIL:" . $mailTo . "\n";
		$msg = MIME::Lite->new(
				 From    => $mailConf->{mailFrom},
				 To      => $mailTo,
				 Subject => $subject,
				 Type    =>'multipart/mixed'
				 );
		$msg->attach(Type     =>'TEXT',
				 Data     =>$message
				 );
		if (-f $path) {
			$msg->attach(Type     =>'application/x-gzip',
				 Path     => $path,
				 Filename => $fileName,
				 Disposition => 'attachment'
				 );
		}
		
		#if ($sendOutput) {
		#	$result = `tar -cvzf output.tar.gz output.out`;
		#	$msg->attach(Type     =>'application/x-gzip',
		#		 Path     => 'output.tar.gz',
		#		 Filename => 'output.tar.gz',
		#		 Disposition => 'attachment'
		#		 );
		#}

		my $smtps = Net::SMTP::SSL->new($mailConf->{mailServer}, 
				               Port => $mailConf->{mailPort},
					       Debug => 1
				               ) or print "$!\n"; 
		# authenticate
		defined ($smtps->auth($mailConf->{mailUser}, $mailConf->{mailPassword}))
		    or die "Can't authenticate: $!\n";
		# send preliminary data
		$smtps->mail("$mailConf->{mailFrom}\n");
		$smtps->to("$mailTo\n");
		$smtps->data();
		#send header
		$smtps->datasend($msg->as_string());
		$smtps->dataend();
		$smtps->quit;
	}
}

#----------------------------------------------------------------------
#--Function senOneEmail
#--Description:
#----Sends only one email to the specified email direction, message and subject
#--Parameters:
#----mailConf: Pointer to a var that contains the email configuration.
#----mailTo: the email direction
#----subject: The email subject
#----message: the email message
sub sendOneEmail {
	local($mailConf, $mailTo, $subject, $message) = @_;

	print "SENDING EMAIL:" . $mailTo . "\n";
	$msg = MIME::Lite->new(
			 From    => $mailConf->{mailFrom},
			 To      => $mailTo,
			 Subject => $subject,
			 Type    =>'multipart/mixed'
			 );
	$msg->attach(Type     =>'TEXT',
		     Data     =>$message
	             );
	
	my $smtps = Net::SMTP::SSL->new($mailConf->{mailServer}, 
				        Port => $mailConf->{mailPort},
					Debug => 1
				        ) or print "$!\n"; 
	# authenticate
	defined ($smtps->auth($mailConf->{mailUser}, $mailConf->{mailPassword}))
	    or die "Can't authenticate: $!\n";
	# send preliminary data
	$smtps->mail("$mailConf->{mailFrom}\n");
	$smtps->to("$mailTo\n");
	$smtps->data();
	#send header
	$smtps->datasend($msg->as_string());
	$smtps->dataend();
	$smtps->quit;
}


#----------------------------------------------------------------------
#--Function senOneEmail
#--Description:
#----Sends only one email to the specified email direction, message and subject
#--Parameters:
#----mailConf: Pointer to a var that contains the email configuration.
#----mailTo: the email direction
#----subject: The email subject
#----message: the email message
sub sendSeveralEmails {
	local($mailConf, $mailTo, $subject, $message, $path) = @_;
        print "Sending emails:" . join(", ", @{$mailTo}) . "\n";
	my $to;
	my $idx = 1;
	my $numMails = scalar(@{$mailTo});
	foreach $mail (@{$mailTo}) {
		$to = $to . $mail;
		if ($idx != $numMails) {
			$to = $to . ", ";
		}
		$idx = $idx + 1;
	}

	$msg = MIME::Lite->new(
			 From    => $mailConf->{mailFrom},
			 To      => $to,
			 Subject => $subject,
			 Type    =>'multipart/mixed'
			 );
	$msg->attach(Type     =>'TEXT',
		     Data     =>$message
	             );
	print "Path:" . $path . "\n";
	if (-f $path) {
			$msg->attach(Type     =>'application/x-gzip',
				 Path     => $path,
				 Filename => $fileName,
				 Disposition => 'attachment'
				 );
	}
	my $smtps = Net::SMTP::SSL->new($mailConf->{mailServer}, 
				        Port => $mailConf->{mailPort},
					Debug => 1
				        ) or print "$!\n"; 
	# authenticate
	defined ($smtps->auth($mailConf->{mailUser}, $mailConf->{mailPassword}))
	    or die "Can't authenticate: $!\n";
	# send preliminary data
	$smtps->mail("$mailConf->{mailFrom}\n");
	$smtps->recipient(@{$mailTo});
	$smtps->data();
	#send header
	$smtps->datasend($msg->as_string());
	$smtps->dataend();
	$smtps->quit;
}

#----------------------------------------------------------------------
#--Function managersEmail
#--Description:
#----Sends emails to the managers according to the manager's xml configuration data
#--Parameters:
#----mail: email configuration data
#----managers: managers configuration data: relationship between modules and manager
#----subject: The email subject
#----svnName: The name of the testSuiteJob report folder
#----tSuiteName: The name of the testSuite report folder
#----reportFile: The name of the testSuiteJob result xml file
sub managersEmail {
	local($mail, $managers, $subject, $reportName, $folderName, $reportFile) = @_;
	
	#get each testSuite manager
	$managersList = $managers->{manager};
	@managerList = $managersList;
	if (ref $managersList eq 'ARRAY') {
		@managerList = @{$managersList};
	}	

	#put all the managerTest in one hash
	my $managerTests = {};
	foreach $manager (@managerList) {
		$managerTestsList = $manager->{managerTest};
		@managerTestList = $managerTestsList;
		if (ref $managerTestsList eq 'ARRAY') {
			@managerTestList = @{$managerTestsList};
		}
		@{$managerTests}{@managerTestList}=();
	}	

	#for each manager creates its summary report
	foreach $manager (@managerList) {
		print "manager name:" . $manager->{managerName} . "\n";
		$summary = ALYA::Reports::doEmailReportManager($folderName, "testSuiteResult.xml", $reportFile, $manager, $managerTests);
		#The summary is only generated if the manager has one of its modules with erros.
		if ($summary) {
			print "sumariing manager name:" . $manager->{managerMail} . "\n";
			sendOneEmail($mail, $manager->{managerMail}, $subject, $summary);
		}
	}

}

sub svnUserChanges {
	local($moduleName, $moduleType, $from, $to, $svnUrl) = @_;
	#open the users file
	my $xmlUsers = new XML::Simple;
	my $dataUsers = $xmlUsers->XMLin("users.xml");

        print "svn log -v $svnUrl/Sources/$moduleType/$moduleName --xml -r $from:$to\n";
	print "TO:" . $to . "\n";
	$output = `svn log -v $svnUrl/Sources/$moduleType/$moduleName --xml -r $from:$to`;
	print "Output:" . $output . "\n";
	@users = ($output =~ m/<author>(\w+)<\/author>/g);
	print "Users:" . join(", ", @users) . "\n";
	#obtain the user's email
	my $emailsTo = {};
	my $withNoEmail = {};
	foreach $user (@users) {
		my $email = 0;
		print "User:" . $user . "\n";
		print "Dumper" . Dumper($dataUsers);
		@usersList = @{$dataUsers->{user}};
		foreach $userData (@usersList) {
			print "UserName:" . $userData->{userName} . "\n";
			if ($userData->{userName} eq $user) {
				$email = $userData->{userEmail};
			}
		}
		if ($email) {
			$emailsTo->{$email} = $email;
		}
		else {
			$withNoEmail->{$user} = $user;	
		}
	}
	@emails = keys %$emailsTo;
	print "Emails:" . join(", ", @emails) . "\n";
	return @emails;
}

sub sendErrorEmail {
	local($moduleName, $moduleType, $from, $to, $svnUrl, $mail, $subject, $folderName, $reportName) = @_;
        my $xml = new XML::Simple;
	if ($moduleType eq "modules" || $moduleType eq "services") {
		#determine which people changed this module
		@emails = svnUserChanges($moduleName, $moduleType, $from, $to, $svnUrl);
		#create the file content
		$report = "------------testSuite Job Result--------------------------";
		$report = $report . "Module $moduleName compilation errors:\n";
		#open the error file
		open FILE, "$folderName/$moduleName.txt" or die "Couldn't open file: $!";
		binmode FILE;
		my $compileError = <FILE>;
		close FILE;
		$report = $report . $compileError;
		#Send the email
		print "Emails2:" . join(", ", @emails) . "\n";
		sendSeveralEmails($mail, \@emails, $subject, $report, "");
	}
	else {
		#test errors
		#determine which people changed this module
		@emails = svnUserChanges($moduleName, "modules", $from, $to, $svnUrl);
		#determine which people changed the kernel
		@emailsKernel = svnUserChanges("", "kernel", $from, $to, $svnUrl);
		#join all the emails
		@userEmails = (@emails, @emailsKernel);
		my %hash   = map { $_, 1 } @userEmails;
		my @emails = keys %hash;
		$report = "------------testSuite Job Result--------------------------";
		$report = $report . "Test module $moduleName errors\n";
		#iterate over each subtest		
		$subData = $xml->XMLin($folderName . "/" . $moduleName . "/" . $moduleName . ".xml" );		
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
			if ($subPassed eq "FAIL") {
				$report = $report . "-------" . $resTest->{testName} . " -->" . $subPassed . "\n";
				if ($resTest->{error}) {
					$report = $report . "---" . $resTest->{error} . "\n";
				}
				#Iterate over the services used in this test and obtains the users that change this services
				$testData = $xml->XMLin("modules/" . $moduleName . "/" . $resTest->{testName} . ".xml" );
				#Number of tests results
				$servicesTest = $testData->{services}->{service};
				@serviceTest = $servicesTest;
				if (ref $servicesTest eq 'ARRAY') {
					@serviceTest = @{$servicesTest};
				}
				#iterate over each test
				foreach $service (@serviceTest) {
					@emailsServices = svnUserChanges($service, "services", $from, $to, $svnUrl);
					#join the emails
					@userEmails = (@userEmails, @emailsServices);
				}						
			}
		}	
		#Send the email
		sendSeveralEmails($mail, \@emails, $subject, $report, $reportName . ".tar.gz");
	}
}

#----------------------------------------------------------------------
#--Function errorsEmail
#--Description:
#----Fore each module, service and test error obtains the users that changed the code from svn, and sends emails
#--Parameters:
#----mail: email configuration data
#----subject: The email subject
#----svnName: The name of the testSuiteJob report folder
#----tSuiteName: The name of the testSuite report folder
#----reportFile: The name of the testSuiteJob result xml file
sub errorsEmail {
	local($mail, $subject, $reportName, $folderName, $reportFile, $to, $svnUrl) = @_;
	print "ReportName:" . $reportName . "\n";
	print "FolderName:" . $folderName . "\n";
	#open the result file
	my $xml = new XML::Simple;
	my $data = $xml->XMLin($reportFile);

	#load the last revision viewed
	open FILE, "lastRevision.txt" or die "Couldn't open file: $!";
  	binmode FILE;
	my $from = <FILE>;
	close FILE;
	$from =~ s/\R//g;

	#module errors
	$modulesCompResult = $data->{testSuiteJobResult}->{moduleCompilation}->{moduleCompilationResult};
	@modCompResults = $modulesCompResult;
	if (ref $modulesCompResult eq 'ARRAY') {
		@modCompResults = @{$modulesCompResult};
	}
	foreach $modCompResult (@modCompResults) {
		if ($modCompResult->{passed} ne "true" && $modCompResult->{moduleName} ne "All") {
			$moduleName = $modCompResult->{moduleName};
			sendErrorEmail($moduleName, "modules" , $from, $to, $svnUrl, $mail, $subject, $folderName, $reportName);
		}                
	}

	#service errors
	$servicesCompResult = $data->{testSuiteJobResult}->{serviceCompilation}->{serviceCompilationResult};
	@serCompResults = $servicesCompResult;
	if (ref $servicesCompResult eq 'ARRAY') {
		@serCompResults = @{$servicesCompResult};
	}
	foreach $serCompResult (@serCompResults) {
		if ($serCompResult->{passed} ne "true" && $serCompResult->{serviceName} ne "All") {
			$serviceName = $serCompResult->{serviceName};
			sendErrorEmail($serviceName, "services" , $from, $to, $svnUrl, $mail, $subject, $folderName, $reportName);
		}                
	}
	
	#test errors
	$modulesResult = $data->{testSuiteResult}->{moduleResults}->{moduleResult};
	@modulesResults = $modulesResult;
	if (ref $modulesResult eq 'ARRAY') {
		@modulesResults = @{$modulesResult};
	}
	foreach $modResult (@modulesResults) {
		if ($modResult->{passed} ne "true") {
			$moduleName = $modResult->{moduleName};
			sendErrorEmail($moduleName, "tests" , $from, $to, $svnUrl, $mail, $subject, $folderName, $reportName);
		}                
	}
		
	#update the last revision
	open FILE, ">lastRevision.txt" or die $!;
	print FILE $to;
	close FILE;
}


#----------------------------------------------------------------------
#--Function sshCopy
#--Description:
#----Function used to copy the report to a remote machine using ssh
#--Parameters:
#----sshConf:Pointer to the ssh configuration xml data.
sub sshCopy {
	local($sshConf, $reportDir) = @_;
	if ($sshConf) {
		if (-d $reportDir) {
		   #$output = `ssh $sshConf->{sshUser}@$sshConf->{sshHost} rm -R $sshConf->{sshFolder}/report_*`;
		   $output = `scp -r $reportDir/* $sshConf->{sshUser}\@$sshConf->{sshHost}:$sshConf->{sshFolder}`;
		}
	}
}


#----------------------------------------------------------------------
#--Function saveError
#--Description:
#----Saves to a file the 30 last lines of a given output
#--Parameters:
#----path:The path where to save the file
#----test:The name of the file to save
#----stdout: The content to save into the file
sub saveError {
	local($path, $test, $stdout) = @_;
		
	my @lines = split(/\n/, $stdout);
	my $numLines = scalar (@lines);
	my $finalLines = join("\n", @lines[($numLines-30) .. ($numLines-1)]);

	#save to disc the output
	open FILE, ">$path/$test.txt" or die $!;
	print FILE $finalLines;
	close FILE;

	return "$test.txt";
}


1;
