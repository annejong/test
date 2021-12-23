#!/usr/bin/env perl

#   Get the data from the bacteriocin database and write it to fasta  
# 	
#	University of Groningen
#	the Netherlands
#	anne.de.jong@rug.nl
#
#   
#  2012 July, Extract data from the Bacteriocin MySQL database
#		4 July 2012 start		
# 

use strict;
use warnings;
use Storable ;
use DBI;



# -----------------------------------------------------  parameters  ---------------------------------------
my $sessiondir = '.';
my $extract ;
my $db = 'BAGEL3_classI';
my $outputfilename = "bacteriocinI.db";
my $usage = "bagel2_db_2_fasta.pl
	-s sessiondir and output folder [default=current folder]
	-db database [default=BAGEL3_classI]
	-o output file name [default=bacteriocinI.db] will be stored in the sessiondir

e.g.  ./bagel2_db_2_fasta.pl -db BAGEL3_classII -o bacteriocinII.db
";

&parseparam();
mkdir $sessiondir ;
&get_data_from_database();



# ---------------------------------------------------------------------------------------- functions ----------------------------------------------------------------------------------------

sub sessionstop {
	open (FILE,'>', "$sessiondir/sessionstop") or die( "Cannot write $sessiondir/sessionstop \n" );;
	print FILE 'DONE!';
	close FILE ;
}


sub get_data_from_database {
	# connect to MySQL database
	my $dbh = DBI->connect("DBI:mysql:dbname=bagel2;host=localhost", "bagel2", "lanti4bac1", {'RaiseError' => 1});
	# Get information on the genome of the input:genomename   e.g. NC00904
	my $sth = $dbh->prepare("SELECT Name, Organism, Subclass, Sequence FROM $db ") ;
	$sth->execute();
	open FILE , '>', $outputfilename or die ("Could not write $outputfilename/n") ;
	while (my @vector = $sth->fetchrow) {
		print FILE ">".$vector[0]." |".$vector[1]." |".$vector[2]."\n".$vector[3]."\n";
		print $vector[0]." |".$vector[1]." |".$vector[2]."\n".$vector[3]."\n";
	}
	$sth->finish;
	$dbh->disconnect();
	close FILE ;
}

sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        $sessiondir				= shift(@arg) if($var eq '-s') ;
		$db						= shift(@arg) if($var eq '-db') ;
		$outputfilename			= shift(@arg) if($var eq '-o') ;
    }
	$outputfilename = "$sessiondir/$outputfilename";
}
