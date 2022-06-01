#!/usr/bin/perl -w
#
# This script builds a makefile for Alya.
#
#--------------------------------------------------------------------------
#
# 1) Check which f90 compiler is available and print header
#
my ($compi,$linki);
my $OS =  `uname`;
chomp($OS);
if ($OS=~/inux/i){    
    print "Linux OS identified. \n";
    print "Using Intel Fortran Compiler ifc. \n";
    $compi = "ifc -quiet -module \$O -c ";
    $linki = "ifc -quiet ";
} elsif ($OS=~/irix/i) {
    print "SGI Irix OS identified. \n";
    print "Using SGI f90 Fortran Compiler. \n";
    #$compi = "f90 -64 -mips4 -r10000 -i8 -r8 -col120 -backslash -c ";
    $compi = 'f90 -64 -mips4 -r10000 -i8 -r8 -col120 -backslash -c $(INCLU)';
    $linki = "f90 -64 -mips4 -r10000 -i8 -r8 -col120 -backslash ";        
} elsif ($OS=~/ygwin/i) {
    print "CygWin under DOS identified. \n";
    print "Using CygWin g95 Fortran Compiler. \n";
    $compi = "g95 -c ";
    $linki = "g95 ";
}

#
#--------------------------------------------------------------------------
#
# 2) Create a list of files
#
my ($file,$srcdir,$dirs,$i,$j,$k,$ncf) ;
my (@zepdirs,@alldirs,@patdirs,@auxdirs);
my (@allfiles,@comfiles,@objfiles);

$srcdir="./Sources";
@zepdirs=($srcdir."/kernel" , $srcdir."/modules" , $srcdir."/services");
@alldirs=();

$i=0;
foreach $dirs (@zepdirs) {
    opendir(ALL_DIRS,$dirs);
    @auxdirs=readdir(ALL_DIRS);
    foreach (@auxdirs) {$patdirs[$i]=$dirs ."/".$_ ;$i++}
    @alldirs=(@alldirs,@auxdirs);
    closedir(ALL_DIRS);
}

$i=0;
foreach $dirs (@alldirs) {
    if (! (($dirs=~/CVS/) || ($dirs=~/\./)) ) {
        $dirs = $patdirs[$k];
        opendir(ALL_FILES,$dirs);
        @allfiles=readdir(ALL_FILES);
        foreach (@allfiles) {
            if (!($_=~/\#/)) {
                if ($_=~/\.f90$/i) {
                    $i++;
                    @comfiles=(@comfiles,$dirs."/".$_);
                    s/f90/o/i;
                    #print "$dirs $_ \n";
                    @objfiles=(@objfiles,$_);
                } elsif ($_=~/\.f$/i) {
                    $i++;
                    @comfiles=(@comfiles,$dirs."/".$_);
                    s/\.f/.o/i;
                    #print "$dirs $_ \n";
                    @objfiles=(@objfiles,$_);
                }
            }
        }
        closedir(ALL_FILES);                    
    }
    $k++;
}

$ncf=$i;
#
#--------------------------------------------------------------------------
#
# 3) Create a hash of files with subroutines, functions and modules as keys
#
my ($interf);

my (@words);
@words=(" " , " ");

my (%subfiles,%funfiles,%modfiles);
%subfiles=();
%funfiles=();
%modfiles=();

$i=0;
$interf=0;
foreach $file (@comfiles) { 
    open (IN,$file);
    $j=0;
    while (<IN>) {
    	$j++;
        if ( ( ($file =~ /\.f90$/) && (!($_=~/^\s*\!/)) ) || ( ($file =~ /\.f$/ ) && (!($_=~/^\s*c/)) ) ) {
            chomp($_);
            s/^\s+//;
            s/\s+$//;
            @words = split(/(\s+|\s*\(\s*|\s*\)\s*|=)/, lc) ;
            if(@words>0) {
                #print " $#words \n";
                #print " @words \n";
                #for ($k=0; $k < @words; $k++) {
                #    print "$words[$k]"."\n";
                #}
                #print "---------------------------\n";
                for ($k=0; $k < @words; $k++) {
                    if($words[$k]=~/\!/) {
                    # The rest of the sentence is a comment => break loop
                        $k=@words;
                    } elsif($words[$k]=~/\'/) {
                    # This is a character constant => break loop
                        $k=@words;
                    # Subroutines:
                    } elsif( ($words[$k]=~/^interface$/i) ) {
                    	if($k+2 > @words ) { 
                            if( $k<2 ) {
                                #print "Interface in file $file, at line $j \n";
                                #print "@words \n" ;print "$_ \n";
                                $interf=1;
                            } elsif(! ($words[$k-2]=~/^end$/i) ) {
                                #print "Interface (k>1) in file $file, at line $j \n";
                                #print "@words \n" ;print "$_ \n";
                                $interf=1;
                            } elsif( ($words[$k-2]=~/^end$/i) && ($interf==1) ) {
                                #print "End interface in file $file, at line: $j \n";
                                #print "@words \n" ;print "$_ \n";
                                $interf=0;
                            }
                        } else {
                            if( $k<2 ) {
                                %subfiles = (%subfiles,$words[$k+2],$i);
                            } elsif(! ($words[$k-2]=~/^end$/i) ) {
                                %subfiles = (%subfiles,$words[$k+2],$i);
                            }
                        }
                    } elsif( ($words[$k]=~/^subroutine$/i) && ($interf==0) ) {
                    	if($k+2 > @words ) { 
                            if( $k<2 ) {
                                #print "Error in file $file, at line $j \n";
                                #print "@words \n" ;print "$_ \n";
                            } elsif(! ($words[$k-2]=~/^end$/i) ) {
                                #print "Error in file $file, at line $j \n";
                                #print "@words \n" ;print "$_ \n";
                            }
                        } else {
                            if( $k<2 ) {
                                %subfiles = (%subfiles,$words[$k+2],$i);
                            } elsif(! ($words[$k-2]=~/^end$/i) ) {
                                %subfiles = (%subfiles,$words[$k+2],$i);
                            }
                        }
                    # Modules:
                    } elsif ( ($words[$k]=~/^module$/i) && ($interf==0) ) {
                    	if($k+2 > @words) { 
                            if( $k<2 ) {
                                #print "Error in file $file, at line $j \n";
                                #print "@words \n" ;print "$_ \n";
                            } elsif(! ($words[$k-2]=~/^end$/i) ) {
                                #print "Error in file $file, at line $j \n";
                                #print "@words \n" ;print "$_ \n";
                            }
                        } else {
                            if( $k<2 && ( !($words[$k+2]=~/^procedure$/i)) ) {
                                %modfiles = (%modfiles,$words[$k+2],$i);
                            } elsif(! ($words[$k-2]=~/^end$/i) ) {
                                %modfiles = (%modfiles,$words[$k+2],$i);
                            }
                        }
                    # Functions:
                    } elsif ( ($words[$k]=~/^function$/i) && ($interf==0) ) {
                    	if($k+2 > @words) { 
                            if( $k<2 ) {
                                #print "Error in file $file, at line $j \n";
                                #print "@words \n" ;print "$_ \n";
                            } elsif(! ($words[$k-2]=~/^end$/i) ) {
                                #print "Error in file $file, at line $j \n";
                                #print "@words \n" ;print "$_ \n";
                            }
                        } else {
                            if( $k<2 ) {
                                %funfiles = (%funfiles,$words[$k+2],$i);
                            } elsif(! ($words[$k-2]=~/^end$/) ) {
                                %funfiles = (%funfiles,$words[$k+2],$i);
                            }
                        }
                    }
                }
            }
        } 
    }
    close (IN);
    $i++;
}

my (@funs);
@funs = keys %funfiles;
#
# To Debug:
#
#my (@subs,@sfil);
#my (@mods,@mfil);
#my (@funs,@ffil);
#@subs = keys %subfiles;
#@sfil = values %subfiles;
#@mods = keys %modfiles;
#@mfil = values %modfiles;
#@funs = keys %funfiles;
#@ffil = values %funfiles;
#
#print $subfiles{"inirun"};

#
#--------------------------------------------------------------------------
#
# 4) Loop on the procedures tree (it is actually a graph) to decide which 
#    files need to be compiled.
#

my ($mlevel,$klevel,$godown,@level);
$mlevel=1;
for ($i=0; $i <= $ncf; $i++) {
    $level[$i]=0;
}
$i=0;
foreach $file (@objfiles) { 
    if($file=~/Alya/) {
        $level[$i]=1;
    }
    $i++;
}

$klevel=0;
while ($klevel<$mlevel) {
    $i=0;
    $klevel++;
    print "Tree level: $klevel \n";
    #print "Files: \n";
    $godown=0;
    foreach $file (@comfiles) { 
        if($level[$i]==$klevel) {
            #print "    $file \n";
            open (IN,$file);
            $j=0;
            while (<IN>) {
                $j++;
                if ( ( ($file =~ /\.f90$/) && (!($_=~/^\s*\!/)) ) || ( ($file =~ /\.f$/ ) && (!($_=~/^\s*c/)) ) ) {
                    chomp($_);
                    s/^\s+//;
                    s/\s+$//;
                    @words = split(/(\s+|\s*\(\s*|\s*\)\s*|=|\*|\+)/, lc ) ;
                    if(@words>0) {
                        #print " $#words \n";
                        #print " @words \n";
                        #for ($k=0; $k < @words; $k++) {
                        #    print "$words[$k]"."\n";
                        #}
                        #print "---------------------------\n";
                        $k=0;
                        while($k<@words) {
                            if($words[$k]=~/\!/) {
                                # the rest of the sentence is a comment => break loop
                                $k=@words;
                            } elsif($words[$k]=~/\'/) {
                                # this is a character constant => break loop
                                $k=@words;
                            } elsif($words[$k]=~/^call$/i) {
                                # touch file that contains this routine (if it has no been touch)
                                # print "$words[$k+2] \n";
                                if(exists $subfiles{$words[$k+2]}) {
                                    if($level[$subfiles{$words[$k+2]}]==0) {
                                        $level[$subfiles{$words[$k+2]}]=$klevel+1;
                                        $godown=1;
                                    }
                                } else {
                                    #print "Routine in library: $words[$k+2] \n";
                                    #print "File $file \n";
                                    #print "line: @words \n";
                                }
                            } elsif($words[$k]=~/^use$/i) {
                                # touch file that contains this module
                                if(exists $modfiles{$words[$k+2]}) {
                                    if($level[$modfiles{$words[$k+2]}]==0) {
                                        $level[$modfiles{$words[$k+2]}]=$klevel+1;
                                        $godown=1;
                                    }
                                }
                            } elsif( ($words[$k]=~/\(/) ) {
                                # this could be a function
                                if( $k > 1) {
                                    # look for words[$k-1] on the list of functions
                                    for $ff (@funs) {
                                        if($ff eq $words[$k-1]) {
                                            #print "Using function $ff in file $file line $j \n "
                                            if(exists $funfiles{$words[$k-1]}) {
                                                if($level[$funfiles{$words[$k-1]}]==0) {
                                                    $level[$funfiles{$words[$k-1]}]=$klevel+1;
                                                    $godown=1;
                                               }
                                           }
                                       }
                                    }
                                }
                            }
                            $k++;
                        }
                    }
                }
            }
        }
        $i++;
    }
    if($godown==1) {$mlevel++}
    #print "------------------------------------------------------------\n";
}

print "Unused files: \n";
$i=0;
$k=0;
$nlast=0;
foreach $file (@comfiles) { 
    if($level[$i]==0) {
        print "$file \n";
        $k++;
    } else {
        $nlast=$i;
    }
    $i++;
}
print "Total number of unused files: $k \n";

