#!/usr/bin/env perl

##
##  Make a genome table from Progresql:molgenedb
##  compatible with jDataTables
##
##

use warnings;
use strict;
use DBI;


# connect to ProgreSQL for genome data
	my $dbh = DBI->connect("DBI:Pg:dbname=molgenedb;host=localhost", "molgen", "!ggdimme!", {'RaiseError' => 1});
# execute SELECT query
	my $sth = $dbh->prepare("SELECT organismname, organismid FROM organism");
	$sth->execute();
# print the select list	
	my %organisms ;
	while (my @vetor = $sth->fetchrow) { $organisms{$vetor[0]} = $vetor[1] ;	}
	print "<select name=organism_id size=1>\n\t<option value=empty>-- Select a genome from this list ---</option>\n" ;
	foreach my $key (sort keys %organisms) {	print "\t<option value=$organisms{$key}>$key</option>\n";	}
	print "</select>\n";
# close db
	$sth->finish;
	$dbh->disconnect();