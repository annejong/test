#!/usr/bin/env perl

##
##  Make a genome table from Progresql:molgenedb
##  compatible with jDataTables
##
##

use warnings;
use strict;
use DBI;

my $organism_id = $ARGV[0] ;

# connect to ProgreSQL for genome data
	my $dbh = DBI->connect("DBI:Pg:dbname=molgenedb;host=localhost", "molgen", "!ggdimme!", {'RaiseError' => 1});
# execute SELECT query
	my $sth = $dbh->prepare("SELECT organismname FROM organism WHERE organismid=$organism_id");
	$sth->execute();
	my @list = $sth->fetchrow ;
	chomp @list ;
	my $organismname=$list[0];
	#print "$organismname\n";
	$sth = $dbh->prepare("SELECT file FROM dnafeature WHERE organism='$organismname'");
	$sth->execute();
	@list = $sth->fetchrow ;
	chomp @list ;
	my $file=$list[0];
	print "$file\n";
# close db
	$sth->finish;
	$dbh->disconnect();