#!/usr/bin/env perl

#  Wrapper for BAGEL3 program 

# 	Anne de Jong
#	University of Groningen
#	the Netherlands
#	anne.de.jong@rug.nl
#
#  2012-november
#
#	This is the wrapper for the main pipeline of BAGEL3
#	this routine optional unzip the file and subsequently scans all the files and folders
#	From the webserver it just scans and analyse all the files, 
#	but as stand-alone the folder content can be filtered using a regular expression for the file names
#
# 2016-July: Add php code injection security issue for webserver



use strict ;
use warnings ; 
use Data::Dumper ;
use lib "/data/bagel3/lib";
use genomics ;



my $program_dir = '/data/bagel3';
my $make_glimmer_model = 1;
my $sessiondir = "./my_session" ;
my $queryfolder ;
my $webserver = 0 ;
my $regex = ".*";
my $usage = "./bagel3
				-s Sessiondir [default=./my_session]
				-i queryfolder or zip_file
				-webserver [default=0; Meaning that the analysis will de done on the folder from -i, else a query from the webserver is expected] 
				-glimmer make glimmer model [default=1; means model will be made] Usefull for multiple entry FastA files with relative small DNA fragments 
				-r regular expression for identifing files [default=.* meaning all files ]
		e.g.  ./bagel3_wrapper.pl -i examples -r \.fna\$ 		
" ;

&parseparam() ;
my $reference_genome = $queryfolder ; # this is only used when a NCBI genome is analyzed
my $min_scaffold_size = 2000 ;
my $tmp;
my $num_genomes = 0;



# -----------------------------------------------------------------------------------  start main pipeline ---------------------------------------------------------------------------

sub check_php_code_injection {
	# will return true is PHP code is found in the upload zip file
	my $folder = shift ;
	my $result = 0 ;
	my $filenames = &searchpath( $folder,  '.*') ;
	foreach my $file (@$filenames) {
		print "Checking for PHP code: $file\n";
		open(FILE,$file);
		my $content = <FILE>;
		print "==>$content\n";
		if ($content =~ /\<\?php/) {  # check if the file contains php code
			$result = 1 ;
			last ;
		}	
		close(FILE);
	}
	return $result ;
}

mkdir $sessiondir ;

# 1. Prepare and unpacking the query file from the webserver
	if ($webserver) {
		my $unpackfolder = "$sessiondir/query";
		#print "unpackfolder=$unpackfolder\n";
		mkdir $unpackfolder ;

		my $tmp = "unzip -j $sessiondir/$queryfolder -d $sessiondir/query" ;
		$tmp = "unzip -j $sessiondir/$queryfolder -d $sessiondir/query"  		if ($queryfolder =~ m/.*zip$/g ); # this looks double but needed in perl 
		$tmp = "unrar x $sessiondir/$queryfolder $sessiondir/query" 		if ($queryfolder =~ m/.*rar$/g );
		$tmp = "tar -xvf $sessiondir/$queryfolder -C $sessiondir/query" 	if ($queryfolder =~ m/.*tar$/g );
		$tmp = "cp $sessiondir/$queryfolder $sessiondir/query/" 			if ($queryfolder =~ m/.*fasta$/g ); # this is the default if any other extension than zip, rar or tar is found

		if (check_php_code_injection("$sessiondir/query")) {
			&genomics::write_log("$sessiondir/resulttable2.html","\n============ PHP code in upload is not allowed =============\n",'true');	
			&genomics::write_log("$sessiondir/sessionstop","\n============ PHP code in upload is not allowed =============\n",'true');	
		}
		
		if ($queryfolder =~ m/(.*)\/.*gbk$/g ) { # this is a file selected from the organims list on the webserver
			$queryfolder = $1 ;
			$regex = '.fna';
			print "QUERYFOLDER=$queryfolder\n";
		} else {
			$queryfolder = $unpackfolder ;
			print $tmp."\n";
			system($tmp) ;	
		}


	}	
	&genomics::write_log("$sessiondir/Bagel3_wrapper.log","Scanning folder $queryfolder for files: $regex WEBSERVER=$webserver",'true');
	my $filenames = &searchpath( $queryfolder,  $regex) ;
	my $numtotal = scalar(@$filenames) ;


# 2. Scan the files and analyse them using bagel3.pl 
	foreach my $file (@$filenames) {
		my %query ;
		my $queryname ;
		my $key ;
		$num_genomes++ ;
		#print "Analyzing nr $num_genomes/$numtotal ==> $file\n" ;
		&genomics::write_log("$sessiondir/Bagel3_wrapper.log","Analyzing nr $num_genomes/$numtotal ==> $file",'true');
		# Check each individual entry
		my @lines = &genomics::read_lines( $file ) ;
		if ($file =~ m/.*\/(.+)$/) {$queryname = $1 ; } else { $queryname = $file ; } 
		foreach my $line (@lines) {
			if ($line =~ m/^\>(.*)( |)/g) {
				$key=$queryname."___".$1;
				$key = &genomics::properfilename($key) ;
				#print $key."\n";
			} else {
				$query{$key} .= $line;
			}	
		}

		# start screening each fasta dna file
		&genomics::write_log("$sessiondir/stats.txt","QUERY\tDNA length\tprocessed",'false');
		foreach my $key (keys %query) {
			my $dnalength = length($query{$key}) ;
			if ( $dnalength > $min_scaffold_size) { 
				my $queryfile = genomics::properfilename($key);
				open FILE ,">$sessiondir/$queryfile" or die ("could not write to $sessiondir/$queryfile") ;
				print FILE ">$queryfile\n" ;
				print FILE $query{$key} ;
				close FILE ;
				if ($make_glimmer_model) { # make a model for glimmer
					$tmp = "$program_dir/glimmer3_train.pl -s $sessiondir -i $sessiondir/$queryfile"  ;
					#print "$tmp\n"; exit();
					system($tmp) if ($make_glimmer_model);
				}	
				$tmp = "$program_dir/bagel3.pl -s $sessiondir -i '$sessiondir/$queryfile' -r $reference_genome -glimmer $make_glimmer_model"; 
				system($tmp) ;
				unlink ("$sessiondir/$queryfile") ;  # remove the file when done
				# Collect all the results in ALL_resulttable.html
				#$tmp ="cat $sessiondir/resulttable.html >>$sessiondir/ALL_resulttable.html";

				&genomics::write_log("$sessiondir/stats.txt","$key\t$dnalength\tyes",'false');
			} else {
				&genomics::write_log("$sessiondir/Bagel3_wrapper.log","============ DNA fragment to small of $key =============",'true');
				&genomics::write_log("$sessiondir/stats.txt","$key\t$dnalength\tno",'false');
			}	
		}

	}
	$tmp ="cat $sessiondir/resulttable2.html >>$sessiondir/ALL_resulttable.html";
	system($tmp) ;
	$tmp ="cat $sessiondir/AOI.scaffolds.fna >>$sessiondir/ALL_AOI.fna";
	system($tmp) ;
	

&genomics::write_log("$sessiondir/sessionstop","\n============ Analysis done =============\n",'true');



# -----------------------------------------------------------------------------------  functions ---------------------------------------------------------------------------



sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        $sessiondir	= shift(@arg) if($var eq '-s') ;
		$queryfolder= shift(@arg) if($var eq '-i') ;
		$webserver	= shift(@arg) if($var eq '-webserver') ;
        $make_glimmer_model		= shift(@arg) if($var eq '-glimmer') ;
        $regex		= shift(@arg) if($var eq '-r') ;
    }
    die $usage if (!$queryfolder) ;
	#remove the last / from sessiondir to make it universal
	$sessiondir =~ s/\/$// ;
}

sub filesearch {														# use an anonimous function to limit code complexity
    my ($regexp, $path, $files) = @_ ;									# parse the parameters							
    $path .= '/' if ($path !~ /\/$/) ;									# add a slash to the end of the path if it is not provided
    opendir( my $ha, $path ) ;											# open the directory
    my @all = sort map {$path.$_} grep {$_ !~ /^\./} readdir( $ha ) ;	# get all the files and subdirectories in the directory 
    closedir( $ha )	;													# close the directory
    my @tocheck = grep {-d $_ } @all ;									# get the subdirectories which need checking	
    my @toadd   = grep {$_ =~ /$regexp/} @all ;							# get the files to add 	
    foreach my $fn (@toadd ) {
		push @$files, $fn ;												# add the files to the store
    }		
    foreach my $sdir ( @tocheck ) {										# check the subdirectories	
		($regexp, $path, $files) = &filesearch($regexp, $sdir, $files) ;
    }
    return($regexp, $path, $files) ;									# return the results
} 


sub searchpath {	# searches the path and its sub directories for files matching with $regexp
    my ($path, $regexp) = @_ ;									# parse the commandline parameters	
    my $files = [] ;
    ($regexp, $path, $files) = &filesearch($regexp, $path, $files) ;
    return( $files ) ;											# return the previous info
}