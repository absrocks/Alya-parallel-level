#!/usr/bin/perl

############################################################################
#### Checking License Generator
############################################################################

#Generate Valid License
$output = `./generateLicense 20160330 20280425 20180330 3000`;
if(index($output, "Generated Successful") == -1){
	print "Error to generate Alya-License\n";
}
else {
	print "Alya-License generated\t\t OK\n";
}

#Executing the program with the generated License
$output = `./hello Alya-License.lic`;
if(index($output, "Hello World ¬¬") == -1){
	print "Error to Execute the program hello\n";
}
elsif(index($output, "Error") != -1){
	print "License not valid\n";
}
else {
	print "Program Execution\t\t OK\n";
}


############################################################################
#### Checking License Control
############################################################################

#Checking ended license period
$output = `./hello test/expired.lic 2>&1`;
if(index($output, "Error: Your license has expired!") == -1){
	#print $output;
	print "Checking Expired License\t FAIL\n";
}
else {
	print "Checking Expired License\t OK\n";
}

#Checking license out of execution hours
$output = `./hello test/maxHourReached.lic 2>&1`;
if(index($output, "Error: Your license has expired!") == -1){
	print "Checking Max Hour \t\t FAIL\n";
}
else {
	print "Checking Max Hour \t\t OK\n";
}

#Checking Checking Activation Expired
$output = `./hello test/activationExpired.lic 2>&1`;
if(index($output, "Error: Your license has expired!") == -1){
	print "Checking Activation Expired \T FAIL\n";
}
else {
	print "Checking Activation Expired \t OK\n";
}

#Checking Checking Activation Expired
$output = `./hello test/notvalid.lic 2>&1`;
if(index($output, "Error: Your license is not valid!") == -1){
	print "Checking not valid license \t FAIL\n";
}
else {
	print "Checking not valid license \t OK\n";
}

#Test that generate a valid license and execute it many time until
#it reach the max hour.
$output = `./generateLicense 20160330 20280425 20180330 5000`;
for (my $i=0; $i <= 100; $i++) {
	$output = `./hello Alya-License.lic 2>&1`;
}
if(index($output, "Error: Your license has expired!") == -1){
	print "Execution expiration \t\t FAIL\n";
}
else {
	print "Execution expiration \t\t OK\n";
}

$output = `./hello test/cheat.lic 2>&1`;
if(index($output, "Are you trying to cheat me?") == -1){
        print "Cheating test: \t\t\t FAIL\n";
}
else {
        print "Cheating test: \t\t\t OK\n";
}

