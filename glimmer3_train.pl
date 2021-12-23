#!/usr/bin/env perl

#  BAGEL3 program 
# 	Anne de Jong
#	University of Groningen
#	the Netherlands
#	anne.de.jong@rug.nl
#
#  2012-october
#
#  This routine makes a training model for Glimmer3
#	Problem: From relative small fragments the training model will be inaccurate 
#	Solution: Around 2Mb are samples from the scaffold un which the model will be generated

use strict ;
use warnings ;
use lib "/data/bagel3/lib";
use genomics ;

my $sessiondir = "./" ;
my $queryfile ;
my $trainingfasta = 'glimmer.trainingset.fasta';
my $outputfile = "glimmer.icm" ;

my $usage = "./glimmer3_train.pl
		-s Sessiondir [default=current folder]
		-i queryfile  [this is the multiple entries fasta file]
		-t DNA fasta file containing the trainingset  [default=glimmer.trainingset.fasta]
		-o output file [default=glimmer.icm]\n
	e.g.  ./glimmer3_train.pl -i NC_009004.fna 		
" ;
&parseparam();

my $max = 2000000 ;  # DNA sample size
my %conf = genomics::read_conf('/data/bagel3/bagel3.conf');


# -----------------------------------------------------------------------------------  start main pipeline ---------------------------------------------------------------------------

# 1. Extract around 2Mb (or change the $max) if possible
	my @lines = genomics::read_lines("$queryfile");
	chomp @lines ;
	open FILE ,">$sessiondir/glimmer.trainingset.fasta" or die ("Could not write to $sessiondir/$trainingfasta\n") ;
	print FILE ">glimmer fasta for building training model containing $max bases\n";
	my $count = 0 ;
	foreach my $line (@lines) {
		if ($count < $max) { 
		print FILE $line."\n" if ($line !~ m/^\>/g ) ; 
		$count += length($line) ;
		}	
	}	
	close FILE ;

# 2. build training model

	# step1:
	# Find long, non-overlapping orfs to use as a training set
	my $tmp = "$conf{glimmer}/long-orfs -n -t 1.15 $sessiondir/$trainingfasta $sessiondir/glimmer.longorfs > /dev/null 2>&1" ;
	system ($tmp) ;

	# step2:
	# Extract the training sequences from the genome file
	$tmp = "$conf{glimmer}/extract -t $sessiondir/$trainingfasta $sessiondir/glimmer.longorfs > $sessiondir/glimmer.train" ;
	system($tmp) ;

	#step3:
	# Build the icm from the training sequences
	$tmp = "$conf{glimmer}/build-icm -r $sessiondir/$outputfile < $sessiondir/glimmer.train" ;
	system($tmp) ;


sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        $sessiondir = shift(@arg) if($var eq '-s') ;
		$queryfile	= shift(@arg) if($var eq '-i') ;
        $outputfile = shift(@arg) if($var eq '-o') ;
    }
    die "No filename found\n$usage" if (!$queryfile) ;
    #remove the last / from sessiondir to make it universal
    $sessiondir =~ s/\/$// ;
}
