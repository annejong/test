#!/usr/bin/perl

# when you move the BAGEL3 results to another location, use this tool to change the url's of results files

use strict ;
use Data::Dumper ;
use lib "/data/molgentools/lib";
use anne_files ;

# /usr/bagel3/bagel3_replacetool.pl


my $sessiondir = '.';
my $old_url = "http://bagel.molgenrug.nl/bagel3results/home/anne/Wolfgang/results" ;
my $new_url = "http://ngs.molgenrug.nl/wolfgang/results" ;
my $regex	= "html\$";


$old_url = 'http://genome2d.molgenrug.nl/genome2d_results/GSEA_pro/129.125.142.9552l9m65vuh5ujb0p8gr48g0u66967' ;
$new_url = 'http://genome2d.molgenrug.nl/GSEA_Pro';



my $usage = "/usr/bagel3/bagel3_replacetool.pl
				-s Sessiondir [default=current folder]
				-old old url
				-new new url
				";

my $filenames = anne_files::get_files_from_subdirs( $sessiondir, $regex) ;
foreach my $file (@$filenames) {
	print "$file\n";
	my @lines = anne_files::read_lines($file);
	my @result ;
	foreach my $line (@lines) {
		$line =~ s/$old_url/$new_url/g ;
		push @result, $line ;
	}
	anne_files::write_lines($file, @result) ;	
}
			
				
sub parseparam {
    my $var ;
    my @arg = @ARGV ;

    while(@arg) {
        $var = shift(@arg) ;
        $sessiondir = shift(@arg) if($var eq '-s') ;
        $old_url	= shift(@arg) if($var eq '-old') ;
        $new_url	= shift(@arg) if($var eq '-new') ;
    }
}