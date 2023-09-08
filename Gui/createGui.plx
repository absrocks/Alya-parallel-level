#!/usr/bin/perl

#usage modules
use FindBin;                 # locate this script
use lib "$FindBin::Bin/libraries";  # use the parent directory
use File::Copy;
use File::Copy::Recursive qw(dircopy);
use File::Path;
use File::Find;
use File::Compare;
use Data::Dumper;
use POSIX;


#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#--------------------------------MAIN-----------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
$paraViewPath = $ARGV[0];

if (!$paraViewPath) {
	print "You need to pass the path to your paraView build folder(example: /home/user/paraViewBuild)\n";
}
else {
	#copia de seguridad del AlyaGuiSource
	dircopy("AlyaGuiSource", "AlyaGuiSourceCopy");
	#create the GUI in AlyaGui dir
	chdir "createGui";
	$output = `./createGUI.plx`;
        print $output;
	#configure and compile the GUI
	chdir "../";
	rmtree("AlyaGuiBin");
	mkdir("AlyaGuiBin");
	chdir("AlyaGuiBin");
	#-DParaView_DIR:PATH=/home/tartigues/paraviewPruebas/ParaViewBin
	$output = `cmake -DCMAKE_INSTALL_PREFIX:PATH=. -DParaView_DIR:PATH=$paraViewPath ../AlyaGuiSource && make`;
	#copy the alya binaries and utilities to the GUI	
	#copy the alya binaries and utilities to the GUi	
	mkdir("alyaBin") or die "MKDIR failed: $!";
	$output = `pwd`;
	print "pwd:" . $output;
	copy("../../Executables/unix/Alya.x","alyaBin") or die "Copy failed: $!";
	copy("../../Utils/user/alya2pos/alya2pos.x","alyaBin") or die "Copy failed: $!";
	chmod 0755, "alyaBin/Alya.x", "alyaBin/alya2pos.x";
	print $output; 
	print "----------------------------------------------------------\n";
	print "----------------------------------------------------------\n";
	print "Alya Gui Source code generated in: Gui/AlyaGui\n";
	print "Alya Gui binaries generated in: Gui/AlyaGuiBin\n";
	print "----------------------------------------------------------\n";
	print "----------------------------------------------------------\n";
}

exit 0
