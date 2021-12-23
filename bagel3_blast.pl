#!/usr/bin/perl

use strict ;
use Data::Dumper ;


my $sessiondir = './';
my $genome ;
my $query ;
my $program ;
my $outfmt = 0 ;
my $threads = 6 ;
my $max_targets = 1 ;
my $evalue = 0.00001 ;
my $outputfile = 'blast.results';

my $usage = "./webblast
				-s Sessiondir [default=current folder]
				-i query
				-g genome
				-m max_targets [default=1]
				-outfmt  outputformat. default=0 means ' as alignmment' 
				-p program [blastp or blastn]
				-o output file [default=blast.results]\n
		e.g.  ./bagel3_webblast.pl -i query -g genome -p blastp		
\n - Anne de Jong - \n" ;

my $debug = " example :
>Ericin A_AAL15567.1|Bacillus subtilis|I|A|
MSKFDDFDLDVVKVSKQDSKITPQVLSKSLCTPGCITGPLQTCYLCFPTFAKC
>Ericin S_AAL15569.1|Bacillus subtilis|I|A|
MSKFDDFDLDVVKVSKQDSKITPQWKSESVCTPGCVTGVLQTCFLQTITCNCHISK
>Mutacin B-Ny266 _AAL73241.1|Streptococcus mutans|I|A|
FKSWSFCTPGCAKTGSFNSYCC

" ;

&parseparam() ;


# copy the genbank to sessiondir
	my $tmp = "cp $genome $sessiondir".'/genome.fasta'  ;
	print "==>".$tmp."\n";
	system($tmp) ;

# format the database for blast > version 2.2
	$tmp =	'formatdb -i '.$sessiondir."/genome.fasta -p T -l ".$sessiondir."/formatdb.log" if ($program eq 'blastp'); 
	$tmp =	'formatdb -i '.$sessiondir."/genome.fasta -p F -l ".$sessiondir."/formatdb.log" if ($program eq 'blastn'); 
	print "==>".$tmp."\n";
	system($tmp) ;

# perform the blast
	$tmp = "$program -db $sessiondir/genome.fasta -query $sessiondir/$query -out $sessiondir/$outputfile -evalue $evalue -outfmt $outfmt -max_target_seqs $max_targets -num_threads $threads" ;
	print "==>".$tmp."\n";
	system($tmp) ;

# reformat results
	&parse_results();


open FILE, ">$sessiondir/sessionstop" or die ("Could not write $sessiondir/sessionstop\n");
print FILE 'Session ready'  ;
close FILE ;

 
sub parse_results {
	# 1. read the blast resultfile
	open FILE, "<$sessiondir/$outputfile" or die ("Could not read $sessiondir/$outputfile\n");
	my @lines = <FILE>;	
	chomp @lines ;
	close FILE ;

	# 2. take the interresing lines
	my $printline = 0 ;
	my @result ;
	foreach my $line (@lines) {
		$printline = 1 if ($line =~ m/Query=/g) ; 
		$printline = 0 if ($line =~ m/^Lambda/g) ; 
		push @result, $line if ($printline) ;	
	}
	
	# 3. write the filtered lines as HTML format
	my $url = "http://www.uniprot.org/uniprot";
	my $cys = "<b><font color=red>C</font></b>";
	open FILE, ">$sessiondir/results.html" or die ("Could not read $sessiondir/results.html\n");
	print FILE "<FONT FACE='Courier New'>\n";
	foreach my $line (@result) {
		$line =~ s/\ /&nbsp/g ;
		$line =~ s/C/$cys/g if ($line =~ m/^Query&nbsp/g or $line =~ m/^Sbjct&nbsp/g) ;
		if ($line =~ m/Query=/g) { 
			print FILE "<hr><font color=blue><b>$line</b></font><br>\n" ; 
		} elsif ($line =~ m/^\>/g) {
			if ($line =~ m/^\>(.*)\|.*/g) {
				print FILE "<b><a href=$url/$1 target=_blank>$line</a> </b><br>\n";
			} else {
				print FILE "<b>$line</b><br>\n";
			}	
		} else { 
			print FILE $line."<br>\n";
		}	
	}
	print FILE "</FONT>\n";
	close FILE ;
}
 
  
sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        $sessiondir = shift(@arg) if($var eq '-s') ;
		$query	= shift(@arg) if($var eq '-i') ;
		$genome	= shift(@arg) if($var eq '-g') ;
		$max_targets	= shift(@arg) if($var eq '-m') ;
		$outfmt	= shift(@arg) if($var eq '-outfmt') ;
		$program	= shift(@arg) if($var eq '-p') ;
        $outputfile = shift(@arg) if($var eq '-o') ;
    }
    die "$usage" if (!$genome or !$program or !$query) ;
	$sessiondir =~ s/\/$// ;  # remove the last /
}