#!/usr/bin/perl

$mak    =  'makefile';
$dot    = 'makefile.dot';
# 
# Check colors at: http://wingraphviz.sourceforge.net/wingraphviz/language/colorname.htm
#
my @colsub  = ();
$colsub[ 0] = 'red';
$colsub[ 1] = 'blue';
$colsub[ 2] = 'green';
$colsub[ 3] = 'yellow';
$colsub[ 4] = 'purple';
$colsub[ 5] = 'orange';
$colsub[ 6] = 'lightsalmon';
$colsub[ 7] = 'lightblue';
$colsub[ 8] = 'orchid';
$colsub[ 9] = 'lightyellow';
$colsub[10] = 'seagreen1';
$colsub[11] = 'wheat';
$colsub[12] = 'thistle2';
$colsub[13] = 'yellow3';
$colsub[14] = 'yellowgreen';
$colsub[15] = 'seagreen';
$colsub[16] = 'royalblue';
#
# Usage
#
print "\n";
print "Output results in dot format in file $dot \n";
print "Visualize the resulrs e.g. with Graphviz\n";
print "\n";

open ($MAKDOT, '>', $dot);
print $MAKDOT "digraph graphname { \n";
#
# Loop over arguments
#
my $n=0;
foreach $subru (@ARGV) {

    $outfile= 'outfile-'. $subru ;    
    open (MAK,$mak);
    #
    # Loop over deps
    #
    @depes = ();
    $i=0;
    while (<MAK>) {
	@rgmod = split;
	if ($rgmod[0]=~/$subru.o\:/) {
	    $i=1;
	}
	if ($i==1) {
	    if ($_=~/FCOPT/) {
		$i=3;
	    }
	    if ($i==1) {
		if (!($_ =~/ources/)) {
		    chop;
		    s/\.o//;
		    s/\$O//;
		    s/\///;
		    s/\\//;
		    tr/ //d;
		    @depes=(@depes,$_);
		}
	    }
	}
    }
    #
    # Eliminate duplicates
    #
     sub uniq {
	my %seen;
	grep !$seen{$_}++, @_;
    }
    my @unique = uniq(@depes);    
    #%hash   = map { $_ => 1 } @depes;
    #@unique = keys %hash;
    #
    # Print output
    #
    print $MAKDOT "$subru  [style=filled, fillcolor=$colsub[$n]]\n";    
    foreach (@unique) {
	print $MAKDOT "$subru -> $_\; \n";
    }
    #
    # Close
    #    
    ++$n;
    close (MAK,$mak);
}

print $MAKDOT "} \n";

exit;
 
 
