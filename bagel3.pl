#!/usr/bin/env perl

#  BAGEL3 program 
# 	Anne de Jong and Auke van Heel
#	University of Groningen
#	the Netherlands
#	anne.de.jong@rug.nl
#
# 	v0.1 2012-october	- first full replacement of BAGEL2 by BAGEL3
#	v1.0 2013 april 10	- major update and first deplay version
#
#  This is the main pipeline of BAGEL3
# 


use strict ;
use warnings ;
use lib "/data/bagel3/lib";
use bagel3_graphics ;
use genomics ;

my %conf = read_conf('/data/bagel3/bagel3.conf') ;


# -----------------------------------------------------------------------------------  global vars ---------------------------------------------------------------------------
my $sessiondir = "./my_session/" ;
my $queryfile ;
my $queryname ;
my $reference_genome ;
my $outputfile = "results.txt" ;
my $is_circular = 1 ;
my $glimmer_train = 0 ;
my $make_graphics = 1 ;
my %hmmcontext ;
my $usage = "./bagel3
				-s Sessiondir [default=current folder]
				-i queryfile
				-r NCBI reference genome
				-is_circular [default=1 for yes]
				-glimmer [default=0, 1 means that the training model is present in the sessiondir; normally made by bagel3_wrapper.pl]
				-graphics [default=1, 1 means that a graphic is drawn of the AOI]
				-o output file [default=results.txt]\n
		e.g.  ./bagel3.pl -s /tmp/BAGEL3WRAPPER -i genomes/NC_009004.fna 		
" ;

&parseparam() ;


# -----------------------------------------------------------------------------------  start main pipeline ---------------------------------------------------------------------------


# initialize session
	# removing old files
	my @files = ('hmm_search.log','hmm_primary_search.log','AOI_tmp_orfs.fasta','PROT_small_orfs','PROT_orfs','AOI.all.faa','all_proteins.faa','AOI.putative_small_orfs.faa','AOI.combined_pfam.context.fasta','AOI.proteins.fasta','AOI.combined_pfam.filtered.fasta','lantibiotics.fasta');
	foreach my $file (@files) { unlink "$sessiondir/$file"; }
	# extract a queryname from the query filename
	if ($queryfile =~ m/.*\/(.+)$/) {$queryname = $1 ; } else { $queryname = $queryfile ; } 
	mkdir $sessiondir ;
	my $tmp ;
	my $query_DNA = read_DNA($queryfile);
	my $queryLen = length($query_DNA);
	my $logfile = "$sessiondir/$conf{logfile}" ;
	my $sessionprogress = "$sessiondir/sessionprogress" ;
	if ($conf{bacteriocinI_db} eq "SQL") {
		# Retrieve the lantibiotics from the MySQL database
		$tmp = "$conf{bagel3dir}/bagel3_db_2_fasta.pl -s $sessiondir" ;
		system($tmp) ;
		$tmp =	"formatdb -i $sessiondir/bacteriocinI.db -p T -l $sessiondir/formatdb.log";
		system($tmp) ;
		$conf{bacteriocinI_db} = "$sessiondir/bacteriocinI.db";
	} 	
	genomics::write_log($logfile,"\n-------------------------------------------------------------------- Start analysing $queryname -----------------------------------------",'true') ;
	genomics::write_log($sessionprogress,"Start analysing $queryname | $queryLen bp<br>",'false');
	genomics::write_log($logfile,"Genome size is $queryLen bp",'true') ;

	
# 1. Call all orfs for primary context screening
	genomics::write_log($logfile,"\n------------ step 1: ORF calling ------------\n",'true');
	genomics::write_log($sessionprogress,"ORF calling<br>",'false');
	$tmp = "$conf{bagel3dir}/orf_calling.pl -s $sessiondir -i $queryfile -min $conf{protein_minlen} -max $conf{protein_maxlen} -o $conf{orfs}" ;
	system($tmp) ; 
	if (!-e "$sessiondir/PROT_$conf{orfs}") { &genomics::write_log($logfile, "\nNo ORFS found in this sequence\n",'true') ; exit; }
	

# 2. Blast all orfs against bacteriocin database I, II and III on the basis of a glimmer orf calling
	genomics::write_log($logfile,"\n------------ step 2 Blast all ORFs against BacteriocinII and III database------------\n",'true');
	&genomics::write_log($sessionprogress,"Blast all ORFs against Bacteriocin I, II and III database<br>",'false');
	&genomics::write_string("$sessiondir/query.fna",">QUERY start=1 end=$queryLen\n$query_DNA\n") ;
	my %all_orfs = &orf_calling_glimmer("$sessiondir/query.fna", "all_orfs.ffn", "all_proteins.faa","rbs=off");
	my %blast_bacteriocinI  = blast_all_proteins("$conf{bagel3dir}/$conf{bacteriocinI_db}","$sessiondir/all_proteins.faa", $conf{blast_evalue_bacteriocinI});
	my $count = scalar (keys %blast_bacteriocinI) ;
	&genomics::write_log($logfile,"Number of bacteriocins found on the basis of blast against BacteriocinI database: $count",'true');	
	my %blast_bacteriocinII = blast_all_proteins("$conf{bagel3dir}/$conf{bacteriocinII_db}","$sessiondir/all_proteins.faa", $conf{blast_evalue_bacteriocinII});
	my $count = scalar (keys %blast_bacteriocinII) ;
	&genomics::write_log($logfile,"Number of bacteriocins found on the basis of blast against BacteriocinII database: $count",'true');	
	my %blast_bacteriocinIII = blast_all_proteins("$conf{bagel3dir}/$conf{bacteriocinIII_db}","$sessiondir/all_proteins.faa", $conf{blast_evalue_bacteriocinIII});
	$count = scalar (keys %blast_bacteriocinIII) ;
	&genomics::write_log($logfile,"Number of bacteriocins found on the basis of blast against BacteriocinIII database: $count",'true');	


	
# 3. Search primary PFAMS in context genes
	genomics::write_log($logfile,"\n------------ step 3 Search primary context genes------------\n",'true');
	&genomics::write_log($sessionprogress,"Search primary context genes<br>",'false');
	my %AOI_identified_primaryPFAMrules = search_primary_PFAM_rules() ;


	
# 4. Select Area Of Interrest on the basis of primary context genes
	genomics::write_log($logfile,"\n------------ step 4: Finding Area Of Interrest [AOI]------------",'true');
	&genomics::write_log($sessionprogress,"Finding Area Of Interrest [AOI]<br>",'false');
	my %AOI ;
	my $AOI_count 	= Select_AOI();
	my $primary_pfams_count = keys %AOI_identified_primaryPFAMrules;
	if (scalar %AOI eq 0) {
		&genomics::write_log($logfile,"===> No AOI found <===\n",'true') ;
		&genomics::write_log("$sessiondir/No_AOI_found.html","No hits found in $queryname<br>\n",'false') ;
		exit ;
	}	


# 5. Orf calling on AOI using glimmer or prodigal
	genomics::write_log($logfile, "\n------------ step 5 ORF calling on AOI------------",'true');
	&genomics::write_log($sessionprogress,"ORF calling on AOI<br>",'false');
	if ($conf{orfcalling_prog} eq 'glimmer') {
		my %glimmer = &orf_calling_glimmer("$sessiondir/AOI.scaffolds.fna", "AOI.$conf{orfcalling_prog}.ffn", "AOI.$conf{orfcalling_prog}.faa","rbs=on");
	} else {
		my %prodigal = &orf_calling_prodigal() ;
	}	

	
# 6. smallORFs calling on AOI using bagel3 orf calling program
	genomics::write_log($logfile, "\n------------ step 6 Finding small orfs in AOI------------",'true');
	&genomics::write_log($sessionprogress,"Finding small orfs in AOI<br>",'false');
	my %small_orfs = &putative_small_orfs_in_AOI('AOI.putative_small_orfs.faa') ;

	
# 7. Filter the small orf on the basis of the rules in table table_small_orf_requirement.txt	
	genomics::write_log($logfile, "\n------------ step 7 Identifying interesting small orfs in AOI------------",'true');
	&genomics::write_log($sessionprogress,"Identifying interesting small orfs in AOI<br>",'false');	
	&filter_small_orf_on_class_table()  if ( $primary_pfams_count>0 ) ;
	

# 8. combine the two files in one FastA file: AOI.all.faa
	$tmp = "cat $sessiondir/AOI.$conf{orfcalling_prog}.faa  $sessiondir/AOI.putative_small_orfs.faa >$sessiondir/AOI.all.faa"; 
	system($tmp);
	&filter_file_for_duplicate_orfs("$sessiondir/AOI.all.faa");
	

	
# 9. Annotate proteins in AOI on the basis of PFAM domains
	&genomics::write_log($logfile,"\n------------ step 9 Annotating proteins ------------",'true');
	&genomics::write_log($sessionprogress,"Annotating proteins<br>",'false');

	print "AOI_publicPFAM\n";
	my %AOI_publicPFAM = search_public_PFAMS("$sessiondir/AOI.all.faa", $conf{pfam_db}) ;
	my $AOI_publicPFAMcount = scalar keys %AOI_publicPFAM ;	
	&genomics::write_log($logfile,"Number of ORFs having a public PFAM domain: $AOI_publicPFAMcount",'true');

	print "AOI_bacteriocinPFAM\n";
	my %AOI_bacteriocinPFAM = search_PFAMS("$sessiondir/AOI.all.faa", "$conf{bagel3dir}/tables/bacteriocin_hmm.txt") ;
	my $AOI_bacteriocinPFAMcount = scalar keys %AOI_bacteriocinPFAM ;	
	&genomics::write_log($logfile,"Number of ORFs having a bacteriocin PFAM domain: $AOI_bacteriocinPFAMcount",'true');

	print "AOI_contextPFAM\n";
	my %AOI_contextPFAM = search_PFAMS("$sessiondir/AOI.all.faa", "$conf{bagel3dir}/tables/context_hmm.txt") ;
	my $AOI_contextPFAMcount = scalar keys %AOI_contextPFAM ;	
	&genomics::write_log($logfile,"Number of ORFs having a context PFAM domain: $AOI_contextPFAMcount",'true');

	print "AOI_primaryPFAM\n";
	my %AOI_primaryPFAM = search_PFAMS("$sessiondir/AOI.all.faa", "$conf{bagel3dir}/tables/primary_hmm.txt") ;
	my $AOI_primaryPFAMcount = scalar keys %AOI_primaryPFAM ;	
	&genomics::write_log($logfile,"Number of ORFs having a primary PFAM domain: $AOI_primaryPFAMcount",'true');

	
# 10. Combine all orfs in AOI.all.faa with the PFAM: AOI.all.faa and AOI.PFAM.results and write to: AOI.combined_pfam.tab
	&genomics::write_log($logfile,"\n------------ step 10 Combining ORFs ------------",'true');
	&genomics::write_log($sessionprogress,"Combining ORFs<br>",'false');
	my %result_table = combine_orf_pfam();
	# headers of the result_table ;
	my @headers = genomics::get_second_hashkey(%result_table) ;

	
# 11. blast the proteins to identify bacteriocins
	&genomics::write_log($logfile,"\n------------ step 11 Blast the proteins to identify bacteriocins ------------",'true');
	&genomics::write_log($sessionprogress,"Blast against BacteriocinI database<br>",'false');
	&genomics::write_log($logfile,"BacteriocinI blast hits:",'true');
	my $AOI_blast_bacteriocinI_count = blast_proteins("$conf{bagel3dir}/$conf{bacteriocinI_db}",'bacteriocinI',"$sessiondir/AOI.proteins.fasta");
	&genomics::write_log($logfile,"Number of bacteriocins found on the basis of blast against BacteriocinI database: $AOI_blast_bacteriocinI_count",'true');	
	&genomics::write_log($sessionprogress,"Blast against BacteriocinII database<br>",'false');
	&genomics::write_log($logfile,"BacteriocinII blast hits:",'true');
	my $AOI_blast_bacteriocinII_count = blast_proteins("$conf{bagel3dir}/$conf{bacteriocinII_db}",'bacteriocinII',"$sessiondir/AOI.proteins.fasta");
	&genomics::write_log($logfile,"Number of bacteriocins found on the basis of blast against BacteriocinII database: $AOI_blast_bacteriocinII_count",'true');		
	&genomics::write_log($sessionprogress,"BacteriocinII blast hits:",'false');
	my $AOI_blast_bacteriocinIII_count = blast_proteins("$conf{bagel3dir}/$conf{bacteriocinIII_db}",'bacteriocinIII',"$sessiondir/AOI.proteins.fasta");
	&genomics::write_log($logfile,"Number of bacteriocins found on the basis of blast against BacteriocinIII database: $AOI_blast_bacteriocinIII_count",'true');	
	&genomics::write_log($logfile,"Context ORFS:",'true');	
	my $AOI_blast_contextORFScount = blast_contextORFs();
	#my %AOI_unirefBLAST = blast_uniref();
	
	
# 12. Remove ovelapping small_orfs OR small_orfs without PFAM or Blast hit in AOI in file AOI.combined_pfam.tab 
	&genomics::write_log($logfile, "\n------------ step 12 Filter small orfs ------------",'true');
	&genomics::write_log($sessionprogress,"Filter small orfs<br>",'false');
	&filter_small_orfs_hits_or_overlap();


# 13. Add blast results 
	&genomics::write_log($logfile,"\n------------ step 13 Add blast results ------------",'true');
	&genomics::write_log($sessionprogress,"Add blast results<br>",'false');
	&add_blast_results();
	

# 14. # remove orfs with similar end position or orfs using the a wrong RBS 
	&genomics::write_log($logfile,"\n------------ step 14 Remove duplicate orfs ------------",'true');
	&genomics::write_log($sessionprogress,"Remove duplicate orfs<br>",'false');
	&remove_duplicate_orfs();

	
# 15. Add physical properties 
	&genomics::write_log($logfile,"\n------------ step 15 Add physical properties ------------",'true');
	&genomics::write_log($sessionprogress,"Add physical properties <br>",'false');
	&add_properties();

# 16. Add small orf classes 
	&genomics::write_log($logfile,"\n------------ step 16 Adding small orf classes ------------",'true');
	&genomics::write_log($sessionprogress,"Adding small orf classes<br>",'false');
	&add_smallorf_class() ;
	&search_leaders();	

# 17. Blast reference genome	
	&genomics::write_log($logfile,"\n------------ step 17 Add NCBI references ------------",'true');
	&genomics::write_log($sessionprogress,"Add NCBI references<br>",'false');
	&blast_reference_genome();

	
# 18. APPEND results to ALL_results.container  (txt and html)
	&genomics::write_log($logfile,"\n------------ step 18 Add results to ALL_results.container.txt ------------",'true');
	&genomics::write_log($logfile,"$queryname analysis data files are stored in folder $sessiondir",'true');
	@headers = ('smallorf_class','queryname','AOI','AOIstart','orf','Gene_start','Gene_end','strand','blast_context','reference','bacteriocinI','bacteriocinII','bacteriocinIII','blast_uniref50',
				'pfam','pfamname','pfam_bacteriocin','pfam_context','dist_proximity_PFAM_bp','CysThrSer','charged_aa','pI','leader_name','processing_site',
				'rbs','length','pfam_primary','find_processing_site','leader_seq','sequence') ;
	&write_result_table("$sessiondir/AOI.result.filtered.tab", @headers);
	$tmp = "cat $sessiondir/AOI.result.filtered.tab >>$sessiondir/ALL_result.container.txt"; 
	system($tmp);
	&genomics::table2html("$sessiondir/ALL_result.container.txt","$sessiondir/ALL_result.container.html") ;

	
# 19. Write Web Table with links to results	
	&genomics::write_log($logfile,"\n------------ step 19 Create HTML table ------------",'true');
	&html_link_table() ;	

	
# 20. Draw graphics	
	&genomics::write_log($logfile,"\n------------ step 20 Adding graphics ------------",'true');
	$tmp = "cp $conf{bagel3dir}/Legend_BAGEL3_graphics.png $sessiondir/Legend_BAGEL3_graphics.png"; 
	system($tmp);	
	&make_image();


# 21. done
	&genomics::write_log($logfile, "\n====================================================================== done =========================================================\n",'true') ;	
	&genomics::write_log($sessionprogress,"============ done =============<br><br>",'false');	

	
	
	
	
	
# -----------------------------------------------------------------------------------  functions ---------------------------------------------------------------------------


sub blast_reference_genome {
	my $genome = $reference_genome.".faa" ;  # this should be the file of the genome2D database
	&genomics::write_log($logfile, "Copying reference genome $genome",'true') ;	
	foreach my $key (keys %result_table) { $result_table{$key}{reference} = ''; }
	if (-e $genome) {	# if the reference genome is in the database, copy this to the sessionfolder
		my $ref_genome = "$sessiondir/ref_genome.faa" ;
		my $tmp = "cp $genome $ref_genome";
		system($tmp) ;
		$tmp = "formatdb -p T -i $ref_genome";
		system($tmp) ;
		$tmp = my ($db, $fasta, $evalue) ;
		my %blast = blast_all_proteins($ref_genome, "$sessiondir/AOI.all.faa", $conf{blast_evalue_context});
		my $url = 'http://www.ncbi.nlm.nih.gov/protein/' ;
		foreach my $key (keys %blast) {
			my @items = split /\|/, $blast{$key}{Subject} ;
			$result_table{$key}{reference} = "<a href=$url$items[2]>$items[1]</a>" if (defined($result_table{$key}{reference}));
		}
	}
	&clean_table();
}


sub html_link_table {
	my %table ;
	my @rowcolor = ("CCE5FF","99CCFF","3399FF");
	my $odd = 0 ;
	my $url = "http://bagel.molgenrug.nl/bagel3results".$sessiondir;
	$url =~ s/\/tmp\/BAGEL3WRAPPER//g ;
	my @wrapperlines = "<b>Table 1, Putative bacteriocin(s) or modified peptide(s)</b><br>
						<table>\n\t<tr bgcolor=$rowcolor[2]><td colspan=4>$queryname</td></tr>";
	my @wrapperdetailed_report ;
	# Make a summary Table 1 of putative bacteriocins
	my @lines = "<hr><a href=$url/resulttable.html target=_blank>BOOKMARK your session here!</a><br><br>";
	push @lines, "<b>Table 1, Putative bacteriocin(s) or modified peptide(s)</b><br>";
	push @lines, "<table>";
	push @lines, "\t<tr bgcolor=$rowcolor[2]><td>Area Of Interest</td><td>Protein ID</td><td>Bacteriocin type</td><td>Homology</td><td>Detailed report</td></tr>";
	push @wrapperlines, "\t<tr bgcolor=$rowcolor[1]>
										<td>Area Of Interest</td>
										<td>Class</td>
										<td>Identified on the basis of</td>
										<td></td>
									</tr>" ;
	foreach my $key (sort { $result_table{$a}{AOI} cmp $result_table{$b}{AOI} || $result_table{$a}{Gene_start} <=> $result_table{$b}{Gene_start} } keys %result_table ) {
		$result_table{$key}{class} = '';
		$result_table{$key}{class} = 'ClassIII' if ( $result_table{$key}{bacteriocinIII} ne '' ) ; 
		$result_table{$key}{class} = 'ClassII' 	if ( $result_table{$key}{bacteriocinII} ne '' ) ;
		$result_table{$key}{class} = 'ClassI' 	if ( $result_table{$key}{bacteriocinI} ne '' ) ;
		$result_table{$key}{class} = $result_table{$key}{smallorf_class} if ( $result_table{$key}{smallorf_class} ne '' ) ; 
		my $homology = $result_table{$key}{bacteriocinI}.' '.$result_table{$key}{bacteriocinII}.' '.$result_table{$key}{bacteriocinIII};
		my $image = "$result_table{$key}{AOI}"."_$queryname.image.html" ;
		if ($result_table{$key}{class} ne '') {
			if ($odd) { $odd=0; } else { $odd=1; }
			push @wrapperdetailed_report, "\t<tr bgcolor=$rowcolor[0]>
												<td>$result_table{$key}{AOI}</td>
												<td>$AOI{$result_table{$key}{AOI}}{class}</td>
												<td>$AOI{$result_table{$key}{AOI}}{aoi_select}</td>
												<td><a href=$url/$image target=_blank>Detailed report</a></td>
											</tr>" ;
			push @lines, "\t<tr bgcolor=$rowcolor[$odd]>
								<td>$result_table{$key}{AOI}</td>
								<td>$key</td><td>$result_table{$key}{class}</td>
								<td>$homology</td><td>
								<a href=$url/$image target=_blank>Detailed report</a></td>
							</tr>";
		}
	}
	push @wrapperlines, sort(genomics::unique_array(@wrapperdetailed_report)) ; # add one time each AOI
	push @lines, "</table>\n";
	push @lines, "<hr>";

	# Make a link Table 2 to additional data
	my @table2 = "<b>Table 2, Additional data</b><br>";
	push @table2, "<table>";
	#push @table2, "\t<tr bgcolor=$rowcolor[$odd]><td><a href=$url/images.html target=_blank>IMAGEs of the Areas Of Interest (AOIs)</a></td></tr>" ; 
	#if ($odd) { $odd=0; } else { $odd=1; }
	push @table2, "\t<tr bgcolor=$rowcolor[$odd]><td><a href=$url/ALL_result.container.html target=_blank>Results in HTML table</a></td></tr>" ; 
	if ($odd) { $odd=0; } else { $odd=1; }
	push @table2, "\t<tr bgcolor=$rowcolor[$odd]><td><a href=$url/ALL_result.container.txt target=_blank>Result file as tab delimited text file </a></td></tr>" ; 
	if ($odd) { $odd=0; } else { $odd=1; }
	push @table2, "\t<tr bgcolor=$rowcolor[$odd]><td><a href=$url/Bagel3.log target=_blank>Detailed log of the analysis</a></td></tr>" ; 
	if ($odd) { $odd=0; } else { $odd=1; }
	push @table2, 		"\t<tr bgcolor=$rowcolor[$odd]><td><a href=$url/Bagel3_wrapper.log target=_blank>Files that have been analysed</a></td></tr>" ; 
	#push @wrapperlines, "\t<tr bgcolor=$rowcolor[1]><td><a href=$url/Bagel3_wrapper.log target=_blank>Files that have been analysed</a></td></tr>" ; 
	if ($odd) { $odd=0; } else { $odd=1; }
	push @table2, "\t<tr bgcolor=$rowcolor[$odd]><td><a href=$url/stats.txt target=_blank>Summary of the DNA fragments analysed</a></td></tr>" ; 
	if ($odd) { $odd=0; } else { $odd=1; }
	#push @table2, "\t<tr bgcolor=$rowcolor[$odd]><td><a href=$url/AOI.scaffolds.fna target=_blank>DNA sequences of the Areas Of Interest in FastA format</a></td></tr>" ; 
	push @table2, "\t<tr bgcolor=$rowcolor[$odd]><td><a href=$url/ALL_AOI.fna target=_blank>DNA sequences of the Areas Of Interest in FastA format</a></td></tr>" ; 
	push @table2, "</table>\n";
	push @table2, "<hr><a href=$url/ALL_resulttable.html target=_blank>BOOKMARK your session here!</a><br><br>";

	push @wrapperlines, "<tr bgcolor=white><td></td></tr></table>\n";
	
	&genomics::write_lines("$sessiondir/resulttable1.html", @wrapperlines) ;
	&genomics::write_lines("$sessiondir/resulttable2.html", @table2) ;
	&genomics::write_lines("$sessiondir/resulttable.html", @lines) ;
	$tmp ="cat $sessiondir/resulttable1.html >>$sessiondir/ALL_resulttable.html";
	system($tmp) ;
}


sub make_image {
	# make an image of each AOI 
	my @lines ;
	my @rowcolor = ("CCE5FF","99CCFF","3399FF");
	my $odd = 0 ;
	my $font = "$conf{bagel3dir}/$conf{font}";
	&add_context_color_and_name() ;
	my $url = "http://bagel.molgenrug.nl/bagel3results".$sessiondir;
	$url =~ s/\/tmp\/BAGEL3WRAPPER//g ;
	foreach my $AOIkey (sort keys %AOI) {
		my @AOIlines ;
		my @besthitblast ;
		my $filename = "$AOIkey"."_$queryname.png" ;
		my %tmptable = %result_table ;
		# %tmptable will only contain entries of the current $AOIkey
		foreach my $key (keys %tmptable) {
			if ($key !~ m/$AOIkey/g) { delete $tmptable{$key} ; }
		}
		# Add colors to genes
		foreach my $key (keys %tmptable) {
			if ($tmptable{$key}{AOI} !~ m/$AOIkey/g) { delete $tmptable{$key} ; }
		}		
		&genomics::write_log($logfile,"Image $filename of $AOIkey saved",'true');
		push @AOIlines, "<a href=$url/$filename target=_blank>IMAGE of $AOIkey of $queryname</a><br>" ; 			# link to the AOI image
		push @AOIlines, "<img src=$url/$filename alt='IMAGE of $AOIkey' width=1100 height=300>" ; 				# the AOI image
		push @AOIlines, "<img src=$url/Legend_BAGEL3_graphics.png alt='legend' width=200 height=200><br>" ; 	# the legend
		&bagel3_graphics::draw_AOI("$sessiondir/$filename", $font, $queryLen, %tmptable);
		# decide which is the bacteriocin gene
		my @Table3 ;
		push @AOIlines, "Table 1. Identified putative bacteriocin(s) or modified peptide(s)<br>" ;
		push @AOIlines, "<table>\n" ;
		push @AOIlines, "<tr bgcolor=$rowcolor[2]><td>Type</td><td>Protein ID</td><td>protein sequence</td></tr>\n" ; 	
		foreach my $key (sort  keys %tmptable) {
			my $seq = $tmptable{$key}{sequence} ;
			if ($result_table{$key}{leader_seq} ne '') {
				my $leader = substr $seq, 0, length($result_table{$key}{leader_seq});
				my $mature = substr $seq, length($result_table{$key}{leader_seq}) ; 
				$seq = "<u>$leader</u>$mature" ;
			}
#			if ($AOI{$AOIkey}{aoi_select} eq 'blast_bacteriocinI' or $result_table{$key}{smallorf_true} eq 'TRUE' or $result_table{$key}{smallorf_rules_true} eq 'TRUE') {
			if ($result_table{$key}{smallorf_true} eq 'TRUE' or $result_table{$key}{smallorf_rules_true} eq 'TRUE') {
				$seq =~ s/C/\<b\>\<font color\=red\>C\<\/font\>\<\/b\>/g ;
				$seq =~ s/S/\<b\>\<font color\=green\>S\<\/font\><\/b\>/g ;
				$seq =~ s/T/\<b\>\<font color\=green\>T\<\/font\><\/b\>/g ;
				if ($odd) { $odd=0; } else { $odd=1; }
				if ($result_table{$key}{smallorf_tophit_true} eq 'TRUE') {
					push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>$tmptable{$key}{class}</td><td>$key</td><td>$seq</td><td>Best Hit</td></tr>\n" ;
					@besthitblast = blast_one_protein("$conf{bagel3dir}/$conf{bacteriocinI_db}", $key, $result_table{$key}{sequence} );
				} else {	
					push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>$tmptable{$key}{class}</td><td>$key</td><td>$seq</td></tr>\n" ;
				}	
				push @Table3, $key ;
			}
			if ($tmptable{$key}{bacteriocinII} ne '' and $AOI{$AOIkey}{aoi_select} eq 'blast_bacteriocinII') {
				if ($odd) { $odd=0; } else { $odd=1; }
				push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>$tmptable{$key}{class}</td><td>$key</td><td>$seq</td></tr>\n" ; 	
				push @Table3, $key ;
			}
			if ($tmptable{$key}{bacteriocinIII} ne '' and $AOI{$AOIkey}{aoi_select} eq 'blast_bacteriocinIII') {
				if ($odd) { $odd=0; } else { $odd=1; }
				push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>$tmptable{$key}{class}</td><td>$key</td><td>$seq</td></tr>\n" ; 	
				push @Table3, $key ;
			}
			if ($tmptable{$key}{smallorf_class} ne '') {
				push @Table3, $key ;
			}
		}
		push @AOIlines, "</table><br>\n" ;
		
		push @AOIlines, @besthitblast ; 	
		
		# Detailed report table
		my $href_pfam = $conf{href_pfam} ;
		push @AOIlines, "<br>Table 2. Annotation of the context<br>" ;
		my @headers = ('Gene_start','Gene_end','strand','pfam','pfamname','length','blast_context','reference') ;
		push @AOIlines, "<table><tr bgcolor=$rowcolor[2]><td>ProteinID</td><td>".join("</td><td>",@headers)."</td></tr>\n" ; 	
		foreach my $key (sort { $tmptable{$a}{Gene_start} <=> $tmptable{$b}{Gene_start} } keys %tmptable) {
			if ($odd) { $odd=0; } else { $odd=1; }
			my @row ;
			foreach my $header (@headers) {	
				if ($tmptable{$key}{$header} =~ m /(PF\d+)\..*/) { # add href to PFAM domains
					push @row, "<a href=$conf{href_pfam}/".$1." target=_blank\>".$tmptable{$key}{$header}."</a>" ; 
				} else {
					push @row, $tmptable{$key}{$header} ; 
				}
			}
			push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>$key</td><td>".join("</td><td>",@row)."</td></tr>\n" ; 	
		}
		push @AOIlines, "</table><br>\n" ;
		
		# add detailed information og the bacteriocin from the result table
		push @AOIlines, "Table 3. Detailed information of putative bactericon(s) or modified peptide(s)<br>" ;
		@Table3 = genomics::unique_array(@Table3) ; 
		push @AOIlines, "<table>\n" ;
		foreach my $key (@Table3) {
			push @AOIlines, "<tr bgcolor=$rowcolor[2]><td>Protein ID</td><td>$key</td></tr>\n" ; 	
			push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>Leader sequence</td><td>$result_table{$key}{leader_seq}</td></tr>\n" if ($result_table{$key}{leader_seq} ne ''); 	
			if ($odd) { $odd=0; } else { $odd=1; }
			push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>Area of Interrest ID (AOI)</td><td>$result_table{$key}{AOI}</td></tr>\n" ; 	
			if ($odd) { $odd=0; } else { $odd=1; }
			push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>AOI start position</td><td>$result_table{$key}{AOIstart}</td></tr>\n" ; 	
			if ($odd) { $odd=0; } else { $odd=1; }
			push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>Postition of the gene in AOI</td><td>$result_table{$key}{Gene_start} .. $result_table{$key}{Gene_end}</td></tr>\n" ; 	
			if ($odd) { $odd=0; } else { $odd=1; }
			push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>Strand</td><td>$result_table{$key}{strand}</td></tr>\n" ; 	
			if ($odd) { $odd=0; } else { $odd=1; }
			push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>length</td><td>$result_table{$key}{length}</td></tr>\n" ; 	
			if ($odd) { $odd=0; } else { $odd=1; }
			push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>Blast hit in bacteriocin I database</td><td>$result_table{$key}{bacteriocinI}</td></tr>\n" if ($result_table{$key}{bacteriocinI} ne ''); 	
			if ($odd) { $odd=0; } else { $odd=1; }
			push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>Blast hit in bacteriocin II database</td><td>$result_table{$key}{bacteriocinII}</td></tr>\n" if ($result_table{$key}{bacteriocinII} ne ''); 	
			if ($odd) { $odd=0; } else { $odd=1; }
			push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>Blast hit in bacteriocin III database</td><td>$result_table{$key}{bacteriocinIII}</td></tr>\n" if ($result_table{$key}{bacteriocinIII} ne ''); 	
			if ($odd) { $odd=0; } else { $odd=1; }
			push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>Lantibiotic PFAM domain</td><td>$result_table{$key}{pfam_bacteriocin}</td></tr>\n" if ($result_table{$key}{pfam_bacteriocin} ne ''); 	
			if ($odd) { $odd=0; } else { $odd=1; }
			push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>Number of Cys and Thr</td><td>$result_table{$key}{CysThrSer}</td></tr>\n" if ($result_table{$key}{CysThrSer} ne ''); 	
			if ($odd) { $odd=0; } else { $odd=1; }
			push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>pI</td><td>$result_table{$key}{pI}</td></tr>\n" if ($result_table{$key}{pI} ne ''); 	
			if ($odd) { $odd=0; } else { $odd=1; }
			push @AOIlines, "<tr bgcolor=$rowcolor[$odd]><td>RBS</td><td>$result_table{$key}{rbs}</td></tr>\n" if ($result_table{$key}{rbs} ne ''); 	
		}
		push @AOIlines, "</table><br>\n" ;
		
		push @AOIlines, "<a href=$url/ALL_resulttable.html target=_blank>Return to Results table</a><br><hr>";

		
		&genomics::write_lines("$sessiondir/$AOIkey"."_$queryname.image.html", @AOIlines) ;
		push @lines, @AOIlines ;
	}
}




sub add_context_color_and_name {
	my %table = genomics::read_table_to_hash("$conf{bagel3dir}/tables/context_annotation_and_color_rules.txt") ;
	# show the color scheme table and change the rule to a regular expression
	genomics::write_log($logfile,"Color scheme for context genes",'true') ;
	my $regexAND = '.*)(?=.*';
	foreach my $key (sort { $table{$a}{order} <=> $table{$b}{order}  } keys %table )  {
		my @line ;
		foreach my $key2 ( keys %{$table{$key}}) {	push @line, "$key2=$table{$key}{$key2}"; }
		my $row = "($key)\t".join("\t",@line) ;
		genomics::write_log($logfile, $row,'true') ;
		# change the rule to a regular expression rule
		$table{$key}{rule} =~ s/ AND /$regexAND/g ;
		$table{$key}{rule} = "(?=.*$table{$key}{rule}).*" ;
	}
	
	foreach my $key ( keys %result_table ) {
		my $pfam = $result_table{$key}{pfam_context} ;
		$result_table{$key}{contextname} = '' ;
		$result_table{$key}{color} = 'gray' if (!defined($result_table{$key}{color}));
		$result_table{$key}{color} = 'greenblue' if ($result_table{$key}{pfam_context} ne '');
		# Apply the rules from the context_annotation_and_color_rules.txt
		foreach my $key_rules (sort { $table{$a}{order} <=> $table{$b}{order}  } keys %table )  {
				if ( $pfam =~ m/$table{$key_rules}{rule}/ ) {
					if ($result_table{$key}{length} < $table{$key_rules}{max_len} ) {  
						genomics::write_log($logfile, "Color applied to: $key\t$table{$key_rules}{name}\t$table{$key_rules}{color}\t$table{$key_rules}{rule}" ,'true') ;
						$result_table{$key}{contextname}  = $table{$key_rules}{name} ; 
						$result_table{$key}{color} = $table{$key_rules}{color} ; 
					}	
				}
		}
		#$result_table{$key}{color} = 'greenIII' if ($result_table{$key}{smallorf_class} ne '') ;
		$result_table{$key}{color} = 'greenI' if ($result_table{$key}{bacteriocinII} ne '') ;
		$result_table{$key}{color} = 'greenI' if ($result_table{$key}{bacteriocinIII} ne '') ;
		$result_table{$key}{color} = 'greenII' if ($result_table{$key}{smallorf_true} eq 'TRUE' or $result_table{$key}{smallorf_rules_true} eq 'TRUE') ;
		$result_table{$key}{color} = 'greenI' if ($result_table{$key}{smallorf_tophit_true} eq 'TRUE') ;
	}
}


sub search_leaders {
	# scan for leader seqs on the basis of leaderHMMs
	genomics::write_log($logfile,"Search leader HMMs",'false');
	# the pfams are defined in 
	my @motifs = ( 'leaderLanBC', 'leaderLanM');
	my @lines ;
	my @files ;
	# write the lanti proteins and initialize values
	foreach my $key (keys %result_table) { 	
		#push (@lines, ">$key\n$result_table{$key}{sequence}\n") if ($result_table{$key}{lantibiotic} eq 'TRUE') ; 
		push (@lines, ">$key\n$result_table{$key}{sequence}\n") if ($result_table{$key}{smallorf_class} ne '') ; 
		$result_table{$key}{leader_name} = ''; 
		$result_table{$key}{leader_evalue} = 1 ;
		$result_table{$key}{leader_from} = 1 ;
		$result_table{$key}{leader_to} = 1 ;
		$result_table{$key}{processing_site} = '';
		$result_table{$key}{leader_seq} = '';
		$result_table{$key}{find_processing_site} = '';
	}
	my $proteins = "$sessiondir/lantibiotics.fasta" ;
	genomics::write_lines($proteins, @lines) ;
	my $resultdomfile ;
	foreach my $motif ( @motifs ) {
		genomics::write_log($logfile,"Search for leader $motif",'false');
		$resultdomfile = "$sessiondir/$motif.leaderdom.tmp.result";
		push (@files, $resultdomfile) ;
		$tmp = "$conf{hmmsearch} --noali --cpu $conf{hmmsearch_processors} --domtblout $resultdomfile -E $conf{hmm_leader_evalue} $conf{hmm_models_dir}$motif.hmm $proteins >>$logfile" ;
		system($tmp) ;
	}

	# merge all results to one file
	my $resultfile = "$sessiondir/all.leaderdom.result";
	$tmp = "cat $sessiondir/*.leaderdom.tmp.result > $resultfile";
	system($tmp) ;

	# parse the best result to the result_table
	@lines = genomics::read_lines($resultfile);
	foreach my $line (@lines) {
		if ($line !~ m/^\#/g) {
			my @items = split " ", $line ;
			my $key = $items[0] ;
			my $besthit = 1 ;
			if ($result_table{$key}{leader_name} ne '') {
				$besthit = 0 if ($result_table{$key}{leader_evalue} < $items[11]) ;
			} 
			if ($besthit) {	
				$result_table{$key}{leader_name} = $items[3];
				$result_table{$key}{leader_evalue} = $items[11] ;
				$result_table{$key}{leader_from} = $items[17] ;
				$result_table{$key}{leader_to} = $items[18] ;
				$result_table{$key}{processing_site} = $result_table{$key}{leader_from} + 1;
				$result_table{$key}{processing_site} = $result_table{$key}{leader_to} if ($result_table{$key}{leader_name} =~ m/LanM/g );
				$result_table{$key}{leader_seq} = substr $result_table{$key}{sequence}, 0, $result_table{$key}{processing_site} ; 
				$result_table{$key}{find_processing_site} = find_processing_site($result_table{$key}{sequence});
				genomics::write_log("$sessiondir/$conf{logfile}", "$key\t$result_table{$key}{leader_name}\t$result_table{$key}{leader_evalue}\t$result_table{$key}{leader_from}\t$result_table{$key}{leader_to}",'true');
			}

		}
	}
	unlink @files ;  # remove all tmp files
}


sub find_processing_site {
# routine inhireted from BAGEL2 to compare with the new PFAM based prediction
    # return the leader sequence the cut site must be between position 19 and 45
    my ($sequence) = @_;
    my $site = substr($sequence, 17, 35) ;
	my $hit = 0 ;
	if ($site =~ m/(SGA[ES]PR)/g and (length (substr($sequence, 17+pos($site))) > 15) ) { 
		$hit=1; 
	} else { 
		if ($site =~ m/(G[GAS])/g  and (length (substr($sequence, 17+pos($site))) > 15)) { 
			$hit=1; 
		} else { 
			if ($site =~ m/(P[RQS])/g  and (length (substr($sequence, 17+pos($site))) > 15)) { $hit=1; }
		}
	}	
	if ($hit eq 1) {
		return substr($sequence, 0, 17+pos($site))."^".substr($sequence, 17+pos($site)); 
	} else { 
		return '' ; 
	}
}




sub clean_table {
	# sometimes perl does not delete the key properly, therefore this small function is added
	foreach my $key (keys %result_table ) {
		delete $result_table{$key} if (!defined($result_table{$key}{AOI})) ;
	}
}

sub write_result_table {
	my($filename, @headers) = @_;
	open FILE, ">$filename" or die ("Could not write $filename\n") ;
	print FILE "ProteinID\t".join("\t",@headers)."\n" ; 	# write header
	foreach my $key (sort { $result_table{$a}{AOI} cmp $result_table{$b}{AOI} || $result_table{$a}{Gene_start} <=> $result_table{$b}{Gene_start} } keys %result_table ) {
		my @row ;
		push @row, $key ;
		if (defined($result_table{$key}{Gene_start})) {
			foreach my $header (@headers) { push @row, $result_table{$key}{$header} ;}
			print FILE join("\t",@row)."\n" ; 
		}	
	}
	close FILE ;
	&genomics::write_log($logfile,"Intermediate result table created: $filename",'true');
}

sub add_smallorf_class {
	# adding small orf class on the basis of the class properties file 
	# $result_table{$key}{smallorf_true} = 'TRUE' for bacteriocinI / II and lanti_pfam
	# $result_table{$key}{smallorf_rules_true} = 'TRUE' if no smallorf_true
	
	my ($rules, $pfams) = read_table_class_properties() ;
	&add_gene_order() ;
	my %rules = %$rules ;
	my %pfams = %$pfams ;
	my @regexpressions ;
	my @headers = genomics::get_second_hashkey(%rules) ;
	foreach my $key (@headers) {
		if ($key =~ m/regex=(.*)/) { push @regexpressions, $1 ; }
	}
	my $regex = join " and ", @regexpressions ;
	genomics::write_log($logfile,"Regular expressions to check: $regex", 'true');
	# $AOI{$AOIkey}{class} = e.g. Lanthipeptide
	my @truekeys ;
	my @rulekeys ;
	foreach my $key (keys %result_table) {
		$result_table{$key}{smallorf_class} = '';
		$result_table{$key}{smallorf_true} = '';
		$result_table{$key}{smallorf_rules_true} = '' ;	
		$result_table{$key}{smallorf_tophit_true} = '';
		$result_table{$key}{dist_proximity_PFAM_bp} = '';
		$result_table{$key}{dist_proximity_PFAM_order} = '';
		$result_table{$key}{sameori_proximity_PFAM} = '';
		$result_table{$key}{upstream_proximity_PFAM} = '';
		$result_table{$key}{pvalueII} = 1;
		$result_table{$key}{pvaluePFAM} = 1;
		$result_table{$key}{pvalueI} = 1;
		my $AOIkey = $result_table{$key}{AOI} ;
		if ($AOI{$AOIkey}{aoi_select} eq 'primary_pfams') {  # apply the rules only for AOI with primay_pfam 
			if ( ($result_table{$key}{bacteriocinI} ne '' or $result_table{$key}{bacteriocinII} ne '' or $result_table{$key}{pfam_bacteriocin} ne '') and $result_table{$key}{length} <=120) {
				$result_table{$key}{smallorf_true} = 'TRUE' ;
				push @truekeys, $key ;
				genomics::write_log($logfile, "$key in AOI=$AOIkey bacteriocinI=$result_table{$key}{bacteriocinI} II=$result_table{$key}{bacteriocinII} pfam=$result_table{$key}{pfam_bacteriocin}", 'true');
				if ($result_table{$key}{bacteriocinII}	=~ m/.*\[(.+)\].*/) { $result_table{$key}{pvalueII} = $1; } # get the pvalue if present
				if ($result_table{$key}{pfam_bacteriocin} 	=~ m/.*\[(.+)\].*/) { $result_table{$key}{pvaluePFAM} = $1; } # get the pvalue if present
				if ($result_table{$key}{bacteriocinI} 	=~ m/.*\[(.+)\].*/) { $result_table{$key}{pvalueI} = $1; } # get the pvalue if present
			}
			my $ruleskey = $AOI{$AOIkey}{class} ;
			my $score = 0 ;
			my $yes_count = 1 ; # 1 because besides yes it also needs to be within min en max_len
			my $len = length($result_table{$key}{sequence}) ;
			$score++ if ($len >= $rules{$ruleskey}{min_len} and $len <= $rules{$ruleskey}{max_len} ) ;
			foreach my $header (@regexpressions) {
				if ($rules{$ruleskey}{"regex=$header"} eq 'yes') {
					$yes_count++ ;
					$score++ if ($result_table{$key}{sequence} =~ m/$header.*/) ;
				}	
			}
			if ($score >= $yes_count) {
				$result_table{$key}{smallorf_class} = $ruleskey  ;
				my $proximitykey = get_proximity_PFAM($key, $AOIkey, $rules{$ruleskey}{proximity} ) ;
				$result_table{$key}{dist_proximity_PFAM_bp}		= abs($result_table{$proximitykey}{Gene_start} - $result_table{$key}{Gene_start} ) ;
				$result_table{$key}{dist_proximity_PFAM_order}	= abs($result_table{$key}{gene_order} - $result_table{$proximitykey}{gene_order} ) ;
				$result_table{$key}{sameori_proximity_PFAM} 	= 1;
				$result_table{$key}{sameori_proximity_PFAM} 	= 0 if ($result_table{$proximitykey}{strand} ne $result_table{$key}{strand}) ;
				$result_table{$key}{upstream_proximity_PFAM} 	= 0 ;
				$result_table{$key}{upstream_proximity_PFAM} 	= 1 if ($result_table{$key}{gene_order} < $result_table{$proximitykey}{gene_order}  and $result_table{$proximitykey}{strand} eq '+') ;
				$result_table{$key}{upstream_proximity_PFAM} 	= 1 if ($result_table{$key}{gene_order} > $result_table{$proximitykey}{gene_order}  and $result_table{$proximitykey}{strand} eq '-') ;
				
				push @rulekeys, $key ;
				genomics::write_log($logfile, "$key in AOI=$AOIkey fits rules for $ruleskey; yes=$yes_count score=$score proximityrule=$rules{$ruleskey}{proximity}  distancePFAM=$result_table{$key}{dist_proximity_PFAM_bp}", 'true');
			}
		}
	}
	
	my %truehash = map { $_ => $result_table{$_} } @truekeys;
	my @besttruekeys ;
	foreach my $AOIkey (sort %AOI) {
		foreach my $key (sort { $truehash{$a}{pvalueI} <=> $truehash{$b}{pvalueI} ||
							$truehash{$a}{pvaluePFAM} <=> $truehash{$b}{pvaluePFAM} || 
							$truehash{$a}{pvalueII} <=> $truehash{$b}{pvalueII} 	} keys %truehash ) {
			if ( $result_table{$key}{AOI} eq $AOIkey) {
				push @besttruekeys, $key ; 
				$result_table{$key}{smallorf_tophit_true} = 'TRUE';
				last ;
			}
		}
	}		

	
	# only keep the closest to the proximity PFAM
	foreach my $AOIkey (sort %AOI) {
		genomics::write_log($logfile, "Checking AOI=$AOIkey", 'true');
		my $hits = 0 ;
		# check if the AOI contains a true key for blast or pfam hit
		my $blast_pfam = 'FALSE' ;
		foreach my $key (@besttruekeys) { 
			if ( $result_table{$key}{AOI} eq $AOIkey) {
				$blast_pfam = 'TRUE' if ( $result_table{$key}{AOI} eq $AOIkey) ;
				genomics::write_log($logfile, "truekey=$key in AOI=$AOIkey", 'true');
	
			}	
		}
		# If there is no bacteriocinI or pfam_bacteriocin, then add orfs on the basis of proximity
		if ($blast_pfam eq 'FALSE') { 
			
			my %hash = map { $_ => $result_table{$_} } @rulekeys;
			foreach my $key (sort { $hash{$b}{sameori_proximity_PFAM} <=> $hash{$a}{sameori_proximity_PFAM} ||
									$hash{$b}{upstream_proximity_PFAM} <=> $hash{$a}{upstream_proximity_PFAM} || 
									$hash{$a}{dist_proximity_PFAM_order} <=> $hash{$b}{dist_proximity_PFAM_order} 	} keys %hash ) {
				if ( $result_table{$key}{AOI} eq $AOIkey) {
					$result_table{$key}{smallorf_rules_true} = 'TRUE' ;	# add this to the result table
					genomics::write_log($logfile, "$key is selected as most putative candidate; sameori_proximity_PFAM=$hash{$key}{sameori_proximity_PFAM} upstream=$hash{$key}{upstream_proximity_PFAM} distance=$hash{$key}{dist_proximity_PFAM_order} ", 'true');
					last ; # only take the first one
				}	
			}	
		}
	}
}

sub get_proximity_PFAM {
	# returns the closest distance of the orf to an orf with the pfam
	my ($smallorfkey, $AOIkey, $pfam) = @_ ;
	my $result = 100000 ;
	my $resultkey ;
	foreach my $key (keys %result_table) {
		if ($result_table{$key}{AOI} eq $AOIkey and $result_table{$key}{pfam_context} =~ m/$pfam.*/) {
			my $distance = abs($result_table{$smallorfkey}{Gene_start} - $result_table{$key}{Gene_start} ) ;
			if ( $distance < $result ) {
				$result = $distance ;
				$resultkey = $key ;
			}	
		}
	}
	return $resultkey ;
}

sub add_gene_order {
	# add the order of the genes on the basis of the gene_start, to calculate the distances between genes
	my $count = 0 ;
	foreach my $key (sort { $result_table{$a}{AOI} cmp $result_table{$b}{AOI} || $result_table{$a}{Gene_start} <=> $result_table{$b}{Gene_start} } keys %result_table ) {
		$result_table{$key}{gene_order} = $count++ ;
	}
}

sub add_properties {
	foreach my $key ( keys %result_table ) {
		my $AOI = $result_table{$key}{AOI} ;
		$result_table{$key}{charged_aa} = '';
		$result_table{$key}{CysThrSer} = '';
		$result_table{$key}{pI} = ''; 
		if ($result_table{$key}{length} <= 120 ) {
			my $Cys    = genomics::count_char($result_table{$key}{sequence},'C') ;
			my $SerThr = genomics::count_char($result_table{$key}{sequence},'(S|T)') ;
			$result_table{$key}{CysThrSer} = "Cys=$Cys SerThr= $SerThr" if ($Cys and $SerThr) ;
			$result_table{$key}{charged_aa} = genomics::count_char($result_table{$key}{sequence},'(K|R|H|D|E|C|Y)') ;
			$result_table{$key}{pI}	= genomics::calculate_pI($result_table{$key}{sequence}) ;
		}
	}
}



sub remove_duplicate_orfs {
	# remove orfs with similar end position
	my $prev_key ;
	my @delkeys ;
	foreach my $key (sort { $result_table{$a}{AOI} cmp $result_table{$b}{AOI} ||$result_table{$a}{strand} cmp $result_table{$b}{strand} || $result_table{$a}{Gene_start} <=> $result_table{$b}{Gene_start} } keys %result_table ) {
		if (defined($prev_key)) {
			if (defined($result_table{$prev_key}{Gene_end})) {
				my $remove_key = 0;
				if ($result_table{$key}{strand} eq '+') { 
					$remove_key = 1 if ($result_table{$prev_key}{Gene_end} == $result_table{$key}{Gene_end}) ;
				} else { # strand is -
					$remove_key = 1 if ($result_table{$prev_key}{Gene_start} == $result_table{$key}{Gene_start}) ;
				}
				if ($remove_key) { # duplicate found
					# choose for the glimmer over the smallORF 
					if (($result_table{$prev_key}{orfcall} eq 'glimmer') and ($result_table{$key}{orfcall} eq 'smallORF')) { 
						push @delkeys, $key ; }
					elsif (($result_table{$prev_key}{orfcall} eq 'smallORF') and ($result_table{$key}{orfcall} eq 'glimmer')) { 
						push @delkeys, $prev_key ;}
					# choose for the best rbs if both are smallORFs
					elsif (($result_table{$prev_key}{orfcall} eq 'smallORF') and ($result_table{$key}{orfcall} eq 'smallORF')) {
						if ( genomics::RBS_score($result_table{$prev_key}{rbs}) > genomics::RBS_score($result_table{$key}{rbs})) { 
							push @delkeys, $key ; }	else { push @delkeys, $prev_key ; }
					}
				}
			}	
		}
		$prev_key = $key ;
	}
	foreach my $key (@delkeys) { 
		&genomics::write_log($logfile, "$key removed",'true') ;
		delete ($result_table{$key}) ; 
	}
	clean_table();
}

sub filter_file_for_duplicate_orfs {
	my $filename = shift ;
	my $key ="";
	# 1. Read the proteins
	my %proteins ;
	my @lines = &genomics::read_lines( $filename ) ;
	chomp @lines ;
	foreach my $line (@lines) {
		if (length($line)>5) {  # so, this is no an empty line
			if ( $line =~ m/\>(.+).+strand(.).*start=(\d+).*end=(\d+).*RBS=(.*) /g ) {
				$key=$1 ;
				$proteins{$key}{header}=$line ;
				$proteins{$key}{strand}=$2 ;
				$proteins{$key}{Gene_start}=$3 ;
				$proteins{$key}{Gene_end}=$4 ;
				$proteins{$key}{rbs}=$5 ;
				my @items = split /\;/, $key ;  # this is split for sorting purposes
				$proteins{$key}{AOI}=$items[0] ;
				$proteins{$key}{orf}=$items[1] ;
			} else {
				$proteins{$key}{sequence}=$line ;
			}
		}	
	}
	
	my $prev_key = 'NULL';
	my @delkeys ;
	foreach my $key (sort { $proteins{$a}{AOI} cmp $proteins{$b}{AOI} || $proteins{$a}{Gene_start} <=> $proteins{$b}{Gene_start} } keys %proteins ) {
		if ($prev_key ne 'NULL') {
			my $remove_key = 0;
			if ($proteins{$prev_key}{Gene_start} == $proteins{$key}{Gene_start} and $proteins{$prev_key}{Gene_end} == $proteins{$key}{Gene_end} ) {
				# choose for the glimmer over the smallORF 
				if ( $key =~ m/smallORF/g ) { 
					push @delkeys, $key ;
				} else {
					push @delkeys, $prev_key ;
				}	
			}	
		}
		$prev_key = $key ;
	}
	foreach my $delkey (@delkeys) { delete $proteins{$delkey} ; }
	@lines = ();
	foreach my $key ( sort  { $proteins{$a}{AOI} cmp $proteins{$b}{AOI} || $proteins{$a}{Gene_start} <=> $proteins{$b}{Gene_start} } keys %proteins ) {
		push @lines , "$proteins{$key}{header}\n$proteins{$key}{sequence}";
	}
	&genomics::write_lines($filename,@lines) ;

}



sub blast_all_proteins {
	my ($db, $fasta, $evalue) = @_ ;  # $resultlable is the label in the result_table has file
	my $blastresultfile = "$sessiondir/blast.all_proteins.results" ;
	my $tmp = "blastp -outfmt 6 -db $db -query $fasta -max_target_seqs 1 -num_threads $conf{blast_processors} -evalue $evalue -out $blastresultfile" ;
	print "$tmp\n";
	system($tmp) ;
	&genomics::add_header($blastresultfile, "Query\tSubject\tpercent\talign_len\tmistmatches\tgap\tq_start\tq.end\ts_start\ts.end\tevalue\tbit_score");
	my %blast = genomics::read_table_to_hash($blastresultfile) ;
	return %blast ;
}


sub blast_one_protein {
	# just blast one protein and return a html file of the alignment
	my $db = shift ;
	my $seqname = shift ;
	my $fasta = shift ;
	my $seqfile = "$sessiondir/oneprotein.fasta" ;
	&genomics::write_string($seqfile, ">$seqname\n$fasta\n") ;
	my $blastresultfile = "$sessiondir/oneproteinblast.results" ;
	my $tmp = "blastp -outfmt 0 -db $db -query $seqfile -max_target_seqs 1 -num_threads 6 -evalue 1E-05 -out $blastresultfile" ;
	system($tmp) ;
	my $cys = "<b><font color=red>C</font></b>";
	my @lines = genomics::read_lines($blastresultfile);
	chomp @lines ;
	my @results = "<FONT FACE='Courier New'>" ;
	my $printon = 0 ;
	foreach my $line (@lines) {
		$line =~ s/\ /&nbsp/g ; # space to html space
		$line =~ s/C/$cys/g if ($line !~ m/Query=.*/) ; # make Cys red
		$printon = 1 if ($line =~ m/Query=.*/) ;
		$printon = 0 if ($line =~ m/Lambda.*/) ;
		push (@results, $line."<br>") if ($printon and $line ne '');
	}
	push @results, "</FONT>";
	&genomics::write_lines("$sessiondir/oneproteinblast.html",@results);
	return @results ;
}

sub blast_proteins {
	# read AOI, write all proteins to AOI.proteins.fasta and blast them against the $db
	my ($db, $resultlable, $fasta) = @_ ;  # $resultlable is the label in the result_table has file
	#my $fasta = "$sessiondir/AOI.proteins.fasta" ; 
	open FILE, ">$fasta" or die ("Could not write $fasta\n") ;
	foreach my $key ( keys %result_table ) {
		print FILE ">$key\n".$result_table{$key}{sequence}."\n" ; 
	}
	close FILE ;
	my $blastresultfile = "$sessiondir/blast.$resultlable.results" ;
	my $tmp = "blastp -outfmt 6 -db $db -query $fasta -max_target_seqs 1 -num_threads $conf{blast_processors} -evalue $conf{blast_evalue_bacteriocinI} -out $blastresultfile" ;
	system($tmp) ;
	&genomics::add_header($blastresultfile, "Query\tSubject\tpercent\talign_len\tmistmatches\tgap\tq_start\tq.end\ts_start\ts.end\tevalue\tbit_score");
	my %blast = genomics::read_table_to_hash($blastresultfile) ;
	foreach my $key (keys %blast) {
		$result_table{$key}{$resultlable} = $blast{$key}{Subject}."[".$blast{$key}{evalue}."]";
		&genomics::write_log($logfile,$key."\t".$result_table{$key}{$resultlable},'true');
	}
	my $count = 0;
	foreach my $key (keys %result_table) { 
		if (defined($result_table{$key}{$resultlable})) {
			$count ++ ;
		} else {
			$result_table{$key}{$resultlable}=''  ;
		}
	}
	return $count ;
}

sub blast_contextORFs {
	# read AOI and write to context.fasta file
	my $blastcontexturl = "http://www.uniprot.org/uniprot/";
	my $fasta = "$sessiondir/AOI.combined_pfam.context.fasta" ; 
	open FILE, ">$fasta" or die ("Could not write $fasta\n") ;
	foreach my $key ( keys %result_table ) {
		print FILE ">$key\n".$result_table{$key}{sequence}."\n" if ($key !~ m/.*smallORF.*/g);
	}
	close FILE ;
	my $blastresultfile = "$sessiondir/blast.context.results" ;
	my $database = "$conf{bagel3dir}/$conf{context_db}" ;
	$tmp = "blastp -outfmt 6 -db $database -query $fasta -max_target_seqs 1 -num_threads $conf{blast_processors} -evalue $conf{blast_evalue_context} -out $blastresultfile" ;
	system($tmp)  ;
	&genomics::add_header($blastresultfile, "Query\tSubject\tpercent\talign_len\tmistmatches\tgap\tq_start\tq.end\ts_start\ts.end\tevalue\tbit_score");
	my %blast = genomics::read_table_to_hash($blastresultfile) ;
	foreach my $key (keys %blast) {
		$result_table{$key}{blast_context} = "<a href=$blastcontexturl$blast{$key}{Subject}>$blast{$key}{Subject}</a>" ;
		&genomics::write_log($logfile,$key."\t".$result_table{$key}{blast_context},'true');
	}
	my $count = 0;
	foreach my $key (keys %result_table) { 
		$result_table{$key}{blast_context}='' if (!defined($result_table{$key}{blast_context})) ;
		$count ++ ;
	}
	return $count ;
}


sub add_blast_results {
	my $fasta = "$sessiondir/AOI.combined_pfam.filtered.fasta" ; 
	open FILE, ">$fasta" or die ("Could not write $fasta\n") ;
	foreach my $key ( keys %result_table ) {
		print FILE ">$key\n" ;
		print FILE $result_table{$key}{sequence}."\n"; ;
	}
	close FILE ;
	

	# blast uniref50
	if ($conf{blast_uniref50}) {
		&genomics::write_log($logfile, "Blast to uniref50 database",'true') ;	
		my $blastresultfile = "$sessiondir/blast.uniref50.results" ;
		$tmp = "blastp -outfmt 6 -db $conf{uniref_db} -query $fasta -max_target_seqs 1 -num_threads $conf{blast_processors} -evalue $conf{blast_evalue_context} -out $blastresultfile" ;
		system($tmp)  ;
		&genomics::add_header($blastresultfile, "Query\tSubject\tpercent\talign_len\tmistmatches\tgap\tq_start\tq.end\ts_start\ts.end\tevalue\tbit_score");
		my %blast = genomics::read_table_to_hash($blastresultfile) ;
		foreach my $key (keys %blast) {
			$result_table{$key}{blast_uniref50} = $blast{$key}{Subject}."[".$blast{$key}{evalue}."]" ;
		}	
	}	
	
	# Add blast results to the file AOI.result.tab
	my @headers = ('queryname','AOI','AOIstart','orf','Gene_start','Gene_end','strand','blast_context','bacteriocinI','bacteriocinII','blast_uniref50','pfam','pfamname','pfam_bacteriocin','pfam_context','length','rbs','sequence') ;
	open FILE, ">$sessiondir/AOI.result.tab" or die ("Could not write $sessiondir/AOI.result.tab\n") ;
	# write header
	print FILE join("ProteinID\t",@headers)."\n" ;
	foreach my $key (sort { $result_table{$a}{AOI} cmp $result_table{$b}{AOI} || $result_table{$a}{Gene_start} <=> $result_table{$b}{Gene_start} } keys %result_table ) {
		if (defined($small_orfs{$key}{rbs})) { 
			$result_table{$key}{rbs} = $small_orfs{$key}{rbs} ; 
		} else { 
			$result_table{$key}{rbs} = "" if ($result_table{$key}{rbs} eq ''); 
		}
		$result_table{$key}{blast_uniref50} = '' if (!defined($result_table{$key}{blast_uniref50})) ;
		my @row ;
		push @row, $key ;
		foreach my $header (@headers) { push @row, $result_table{$key}{$header} ; }
		print FILE join("\t",@row)."\n" ;
	}
	close FILE ;
	&genomics::write_log($logfile, "Combined results written to AOI.result.tab",'true') ;	
}


sub filter_small_orfs_hits_or_overlap {
	# write all orfs but only the smallORF if it contains a PFAM or a blast
	my @delkeys ;
	foreach my $AOIkey (sort keys %AOI) {
		my $remove_smallorf = 0;
		foreach my $key ( keys %result_table ) {
			if ( $AOIkey eq $result_table{$key}{AOI} ) {
				$remove_smallorf = 1 if ($result_table{$key}{pfam_bacteriocin} ne '' or $result_table{$key}{bacteriocinI} ne '') ;
			}
		}
		if ($remove_smallorf) {
			genomics::write_log($logfile,"$AOIkey contains small orf putative bacteriocin, Removing other small orfs",'false') ;
			foreach my $key ( keys %result_table ) {
				if ( $AOIkey eq $result_table{$key}{AOI} ) {
					push (@delkeys, $key) if ( ($result_table{$key}{pfam_bacteriocin} eq '') and ($result_table{$key}{bacteriocinI} eq '') and ($key =~ m/.*smallORF.*/g)) ; 
				}
			}	
		} else {
			genomics::write_log($logfile,"$AOIkey NO small orf putative bacteriocin, Removing overlapping small orfs",'false') ;
			push @delkeys, remove_ovelapping_orfs($AOIkey) ;
		}
		
	}
	
	# remove the selected keys
	@delkeys = genomics::unique_array(@delkeys) ;
	foreach my $key (@delkeys ) {
		genomics::write_log($logfile,"Removing ".$key."\t".$result_table{$key}{sequence},'false') ;
		delete $result_table{$key} ;
	}	
	
	&remove_identical_orfs() ;
	
	# write file
	my @headers = ('queryname','AOI','AOIstart','orf','Gene_start','Gene_end','strand','blast_context','bacteriocinI','bacteriocinII','pfam','pfamname','pfam_bacteriocin','pfam_context','length','rbs','sequence') ;
	open FILE, ">$sessiondir/AOI.combined_pfam.filtered.tab" or die ("Could not write $sessiondir/AOI.combined_pfam.filtered.tab\n") ;
	print FILE join("\t",@headers)."\n" ; # write header
	foreach my $key (sort { $result_table{$a}{AOI} cmp $result_table{$b}{AOI} || $result_table{$a}{Gene_start} <=> $result_table{$b}{Gene_start} } keys %result_table ) {
		foreach my $header (@headers) { 
			if (defined($result_table{$key}{$header})) { print FILE $result_table{$key}{$header}."\t" ; } else {print FILE "\t";}
		}
		print FILE "\n";
	}
	close FILE ;
}

sub remove_identical_orfs {
	# this routine is made to remove identical small orfs that are both called by glimmer and the bagel_orf_calling method
	my $prev_key ;
	my $firstkey=1;
	my @delkeys ;
	foreach my $key (sort { $result_table{$a}{AOI} cmp $result_table{$b}{AOI} || $result_table{$a}{Gene_start} <=> $result_table{$b}{Gene_start} } keys %result_table ) {
		if ($firstkey) { 
			$prev_key = $key ;
			$firstkey = 0 ;
		} else {
			if ( $result_table{$key}{Gene_start} eq  $result_table{$prev_key}{Gene_start} and  $result_table{$key}{Gene_end} eq  $result_table{$prev_key}{Gene_end} ){
				if ($key =~ m/smallORF/g) { 
					push @delkeys, $key ; 
				} else {
					push @delkeys, $prev_key ; 
				}
			}
		}
	}
	
	# remove the selected keys
	@delkeys = genomics::unique_array(@delkeys) ;
	foreach my $key (@delkeys ) {
		genomics::write_log($logfile,"Removing identical small ORF ".$key."\t".$result_table{$key}{sequence},'false') ;
		delete $result_table{$key} ;
	}	
}

sub remove_ovelapping_orfs {

	my $AOIkey = shift ;
	my $max_overlap = 30 ; # number of bases allowed for overlap
	my @results ;
	my %table2 ;
	# Table2 will contain only the glimmerORFS
	genomics::write_log($logfile,"Checking $AOIkey", 'true') ; 
	foreach my $key ( keys %result_table ) { 
		if ($key !~ m/smallORF.*/g and $result_table{$key}{AOI} eq $AOIkey) { 
			genomics::write_log($logfile,"Added $key to table2", 'true') ; 
			$table2{$key}{Gene_start} = $result_table{$key}{Gene_start} ;
			$table2{$key}{Gene_end}   = $result_table{$key}{Gene_end} ;
			$table2{$key}{AOI}   = $result_table{$key}{AOI} ;
		}
	}
	
	
	# check if the start or end of a smallORF is not ovelapping a long orf, in the same AOI
	foreach my $key (sort { $result_table{$a}{AOI} cmp $result_table{$b}{AOI} || $result_table{$a}{Gene_start} <=> $result_table{$b}{Gene_start} } keys %result_table ) {
		if ($key =~ m/smallORF.*/g ) {
			genomics::write_log($logfile,"key-$key $result_table{$key}{AOI} eq $AOIkey", 'true') ; 
			if ($result_table{$key}{AOI} eq $AOIkey) {
				genomics::write_log($logfile,"Checking key=$key orf=$result_table{$key}{orf}", 'true') ; 
				my $overlap = 0 ;
				# Check if the start or the stop is not in other gene regions within the same AOI
				foreach my $key_mask (keys %table2) {
					if ( $result_table{$key}{Gene_start} > $table2{$key_mask}{Gene_start} and  
						 $result_table{$key}{Gene_start} < $table2{$key_mask}{Gene_end}-$max_overlap ) { $overlap = 1 ; }
					if ( $result_table{$key}{Gene_end}   > $table2{$key_mask}{Gene_start}+$max_overlap and  
						 $result_table{$key}{Gene_end}   < $table2{$key_mask}{Gene_end}-$max_overlap ) { $overlap = 1 ; }
					if ( $result_table{$key}{strand} eq '+' and $result_table{$key}{Gene_end}   == $table2{$key_mask}{Gene_end})   { $overlap = 1 ; }
					if ( $result_table{$key}{strand} eq '-' and $result_table{$key}{Gene_start} == $table2{$key_mask}{Gene_start}) { $overlap = 1 ; }
				}
				genomics::write_log($logfile,"$key from $AOIkey removed due overlap\t".$result_table{$key}{sequence},'true') if ($overlap) ;
				push @results, $key if ($overlap) ;
			}	
		}
	}
	return @results ;
}


sub combine_orf_pfam {
	# this function creates the first results in a hash table
	my %table ;
	my $key ="";
	# 1. Add protein information
	my @lines = &genomics::read_lines("$sessiondir/AOI.all.faa" ) ;
	chomp @lines ;
	foreach my $line (@lines) {
		if (length($line)>5) {  # so, this is no an empty line
			if ( $line =~ m/\>(.+).+strand(.).*start=(\d+).*end=(\d+).*RBS=(.*) /g ) {
				$key=$1 ;
				$table{$key}{strand}=$2 ;
				$table{$key}{Gene_start}=$3 ;
				$table{$key}{Gene_end}=$4 ;
				$table{$key}{rbs}=$5 ;
				my @items = split /\;/, $key ;  # this is split for sorting purposes
				$table{$key}{AOI}=$items[0] ;
				$table{$key}{orf}=$items[1] ;
			} else {
				$table{$key}{sequence}=$line ;
			}
		}	
	}
	# 1. Initialize values
	my @new_fields = ('pfam','pfamname','pfam_bacteriocin','pfamname_bacteriocin','pfam_context','pfamname_context','pfamclass_context','pfam_primary');
	foreach my $key ( keys %table ) {
		foreach my $header (@new_fields) {	
			$table{$key}{$header}	= '';
		}	
	}

	# 2. Add public PFAM domains 
	foreach my $key (keys %AOI_primaryPFAM) {
		$table{$key}{pfam_primary} = $AOI_primaryPFAM{$key}{pfam} ;
	}
	
	# 2. Add public PFAM domains 
	foreach my $key (keys %AOI_publicPFAM) {
		$table{$key}{pfam} = $AOI_publicPFAM{$key}{pfam} ;
		$table{$key}{pfamname} = $AOI_publicPFAM{$key}{pfamname} ;
	}
	# 3. Add lanti PFAM domains 
	foreach my $key (keys %AOI_bacteriocinPFAM) {
		$table{$key}{pfam_bacteriocin} = $AOI_bacteriocinPFAM{$key}{pfam} ;
		$table{$key}{pfamname_bacteriocin} = $AOI_bacteriocinPFAM{$key}{pfamname} ;
	}
	# 4. Add context PFAM domains 
	foreach my $key (keys %AOI_contextPFAM) {
		$table{$key}{pfam_context} = $AOI_contextPFAM{$key}{pfam} ;
		$table{$key}{pfamname_context} = $AOI_contextPFAM{$key}{pfamname} ;
	}

	# write combined table to AOI.combined_pfam.tab
	my @tableheader = ('ProteinID','Organism','AOI','AOI_start','ORF','Gene_start','Gene_end','Strand','PFAM','PFAM_name','LantiPFAM','ContextPFAM','Length','sequence') ;
	open (FILE,">$sessiondir/AOI.combined_pfam.tab") ; 
	print FILE join("\t",@tableheader)."\n" ; 	# write header
	foreach my $key (sort { $table{$a}{AOI} cmp $table{$b}{AOI} || $table{$a}{Gene_start} <=> $table{$b}{Gene_start} } keys %table ) {
		$table{$key}{queryname}	= $queryname ;
		$table{$key}{AOIstart}	= $AOI{$table{$key}{AOI}}{start} ;
		$table{$key}{length}  = length ($table{$key}{sequence}) ;
		print FILE "$key\t$queryname\t$table{$key}{AOI}\t$AOI{$table{$key}{AOI}}{start}\t$table{$key}{orf}\t$table{$key}{Gene_start}\t$table{$key}{Gene_end}\t$table{$key}{strand}\t";
		print FILE "$table{$key}{pfam}\t$table{$key}{pfamname}\t$table{$key}{pfam_bacteriocin}\t$table{$key}{pfam_context}\t$table{$key}{length}\t$table{$key}{sequence}\n" ;
	}
	close FILE ;
	genomics::write_log($logfile, "Combined file for PFAM annotation written to $sessiondir/AOI.combined_pfam.tab\n",'true') ;
	return %table ;
}


sub putative_small_orfs_in_AOI {
	# Call all small orfs in each AOI, the orf calling routine will write 2 file; DNA_small_orfs  and PROT_small_orfs 
	my $filename = shift ; 
	my $startcodons = 'ATG|GTG|TTG';
	$startcodons = 'ATG' ;
	my %table ;
	my $tmp_dnafilename = "$sessiondir/tmp_dna" ;
	open (FILE,">$sessiondir/$filename") ; 
	foreach my $key (sort keys %AOI) {
		if ($AOI{$key}{aoi_select} eq 'primary_pfams') {
			&genomics::write_string( $tmp_dnafilename, $AOI{$key}{dna} ) ;
			genomics::write_log($logfile,"Scanning $key for putative small orfs",'true');
			$tmp = "$conf{bagel3dir}/orf_calling.pl -s $sessiondir -i $tmp_dnafilename -min $conf{lantibiotic_minlen} -max $conf{lantibiotic_maxlen} -startcodons $startcodons -o small_orfs" ;
			system($tmp); 
			my @lines = &genomics::read_lines( "$sessiondir/PROT_small_orfs" ) ;
			chomp @lines ;
			my $tmp_start ;
			my $AOI_ID = '';
			foreach my $line (@lines) { 
				$line =~ s/\>/\>$key;small/ ; # the AOI_ID to the header
				print FILE $line."\n" ; 
				# load this data in the hash table
				if ($line =~ m/\>(.*) strand(.).*start=(\d+).end=(\d+).RBS=(.+) GENE.*/g) {
					$AOI_ID = $1 ;
					$table{$AOI_ID}{AOI} = $key;
					$table{$AOI_ID}{orfcall} = 'smallORF';
					$table{$AOI_ID}{strand} = $2;
					$table{$AOI_ID}{start} = $3; 
					$table{$AOI_ID}{end} = $4 ;
					$table{$AOI_ID}{rbs} = $5;
					$table{$AOI_ID}{rbs_score} = genomics::RBS_score($5) ;
				} else {
					$table{$AOI_ID}{sequence} = $line ;
				}	
			}
		}	
	}	
	close FILE ;
	genomics::write_log($logfile,"Putative smallORFs in AOIs written to $filename",'true') ;
	return %table ;
}


sub filter_small_orf_on_class_table {
	my ($rules, $pfams) = read_table_class_properties() ;
	my %rules = %$rules ;
	my %pfams = %$pfams ;
	my @regexpressions ;
	my @headers = genomics::get_second_hashkey(%rules) ;
	foreach my $key (@headers) {
		if ($key =~ m/regex=(.*)/) { push @regexpressions, $1 ; }
	}
	
	# $AOI{$AOIkey}{class} = e.g. Lanthipeptide
	my @orf_keys ;
	my @orf_delkey ;
	foreach my $orf (keys %small_orfs) {
		my $AOIkey = $small_orfs{$orf}{AOI} ;
		my $ruleskey = $AOI{$AOIkey}{class} ;
		my $score = 0 ;
		my $yes_count = 1 ; # 1 because besides yes it also needs to be within min en max_len
		my $len = length($small_orfs{$orf}{sequence}) ;
		$score++ if ($len >= $rules{$ruleskey}{min_len} and $len <= $rules{$ruleskey}{max_len} ) ;
		foreach my $header (@regexpressions) {
			if ($rules{$ruleskey}{"regex=$header"} eq 'yes') {
				$yes_count++ ;
				$score++ if ($small_orfs{$orf}{sequence} =~ m/$header.*/) ;
			}	
		}
		if ($score >= $yes_count) {
			push (@orf_keys, $orf) ;
		} else {
			push (@orf_delkey, $orf) ;
		}
		genomics::write_log($logfile,"$small_orfs{$orf}{sequence} with SCORE=$score of $yes_count needed, ", 'false');
	}
	
	foreach my $key (@orf_delkey) {	delete $small_orfs{$key} ; }
	
	# write the small ORFs  
	my @orfs ;
	foreach my $key (@orf_keys) { 
		genomics::write_log($logfile,"SAVED: $key\t$small_orfs{$key}{sequence}",'true') ;
		push @orfs, ">$key strand=$small_orfs{$key}{strand} start=$small_orfs{$key}{start} end=$small_orfs{$key}{end} RBS=$small_orfs{$key}{rbs} GENE\n$small_orfs{$key}{sequence}\n" ;
	}
	&genomics::write_lines( "$sessiondir/PROT_small_orfs" , @orfs ) ;
}


sub search_PFAMS {
	my $proteins 	= shift ;
	my $pfam_table 	= shift ;
	my $logfile 	= "$sessiondir/hmm_search.log" ;
	my $table_out ;
	my @files ;
	genomics::write_log($logfile,"Search PFAMs in $proteins. PFAMS based on table $pfam_table",'true');
	# the pfams are defined in a special bagel3 formatted table
	my %pfams = genomics::read_table_to_hash($pfam_table) ;
	foreach my $motif ( keys %pfams ) {
		#genomics::write_log($logfile,"Search for motif:$motif",'true');
		$table_out = "$sessiondir/$motif.hmm.tmp.result";
		push (@files, $table_out) ;
		$tmp = "$conf{hmmsearch} --noali --cpu $conf{hmmsearch_processors} --tblout $table_out -E $pfams{$motif}{threshold} $conf{hmm_models_dir}$motif.hmm $proteins >>$logfile" ;
		system($tmp) ;
	}

	# merge all results to one file
	$table_out = "$sessiondir/all.hmm.tmp.result";
	push (@files, $table_out) ;
	$tmp = "cat $sessiondir/*.hmm.tmp.result > $table_out";
	system($tmp) ;

	# parse the results to table hash
	my %table ;
	my @lines = genomics::read_lines($table_out);
	foreach my $line (@lines) {
		if ($line !~ m/^\#/g) {
			my @items = split " ", $line ;
			my $key = $items[0] ;
			if (defined($table{$key}{pfamname})) {
				$table{$key}{pfamname} .= ";$items[2] [$items[4]]" ;
				if ($items[3] ne '-') {
					$table{$key}{pfam} .= " $items[3]" ;
				} else {	
					$table{$key}{pfam} .= " $items[2]" ;
				}	
				
			} else {
				$table{$key}{pfamname} = $items[2] ;
				$table{$key}{pfam} = $items[2] ;
				$table{$key}{pfam} = $items[3] if ($items[3] ne '-')  ;
			}
			genomics::write_log("$sessiondir/$conf{logfile}", "$key\t$table{$key}{pfamname}\t$table{$key}{pfam}",'false');
		}
	}
	unlink @files ;  # remove all tmp files
	return %table ;
}

sub search_public_PFAMS {
	my $proteins 	= shift ;
	my $pfam_table 	= shift ;
	my $logfile 	= "$sessiondir/hmm_search.log" ;
	my $table_out = "$sessiondir/AOI.PFAM.results";
	my @files ;
	unlink ($table_out) if (-e $table_out) ; # remove old result file because pfam_scan does not do this job
	#$tmp = "$conf{pfam_scan} -cpu $conf{pfamscan_processors} -fasta $proteins -dir $conf{pfam_db} -outfile $table_out" ;
	$tmp = "hmmsearch --tblout $table_out --cpu 12 $conf{pfam_db}/Pfam-A.hmm $proteins >$sessiondir/null";
	print $tmp."\n";
	system($tmp) ;

	# parse the results to table hash
	my %table ;
	my @lines = genomics::read_lines($table_out);
	foreach my $line (@lines) {
		if ($line !~ m/^\#/g) {
			my @items = split " ", $line ;
			if (scalar @items>3) {
				my $key = $items[0] ;
				if (defined($table{$items[0]}{name}) ) {
					$table{$key}{pfam} .= ";$items[3] [$items[4]]" ;
					$table{$key}{pfamname} .= " $items[2]" if ($items[2] ne '-') ;
				} else {
					$table{$key}{pfam} = "$items[3] [$items[4]]" ;
					$table{$key}{pfamname} = "";
					$table{$key}{pfamname} = $items[2] if ($items[2] ne '-')  ;
				}
				genomics::write_log("$sessiondir/$conf{logfile}", "$key\t$table{$key}{pfamname}\t$table{$key}{pfam}",'false');
			}	
		}
	}	
	# my %table ;
	# my @lines = genomics::read_lines($table_out);
	# foreach my $line (@lines) {
		# if ($line !~ m/^\#/g) {
			# my @items = split " ", $line ;
			# if (scalar @items>3) {
				# my $key = $items[0] ;
				# if (defined($table{$items[0]}{name}) ) {
					# $table{$key}{pfam} .= ";$items[5] [$items[12]]" ;
					# $table{$key}{pfamname} .= " $items[6]" if ($items[6] ne '-') ;
				# } else {
					# $table{$key}{pfam} = "$items[5] [$items[12]]" ;
					# $table{$key}{pfamname} = "";
					# $table{$key}{pfamname} = $items[6] if ($items[6] ne '-')  ;
				# }
				# genomics::write_log("$sessiondir/$conf{logfile}", "$key\t$table{$key}{pfamname}\t$table{$key}{pfam}",'false');
			# }	
		# }
	# }
	return %table ;
}

sub read_table_class_properties {
	my %hmm_primary_rules = genomics::read_table_to_hash("$conf{bagel3dir}/tables/AOI_identification_rules.txt") ;
	my %pfams ;
	# Get the unique pfams from the table and convert rule to regular expression
	my $regexAND = '.*)(?=.*';
	foreach my $key (sort keys %hmm_primary_rules) {
		$hmm_primary_rules{$key}{rule} =~ s/ and / AND /g ;
		my @items = split ( / AND /, $hmm_primary_rules{$key}{rule} ) ; 
		foreach my $item (@items) {
			$item =~ s/(\(|\))//g ; # remove brackets
			my @sub_items = split ( /\|/, $item ) ; 
			foreach my $pfam (@sub_items) { $pfams{$pfam} = $hmm_primary_rules{$key}{threshold} ; }
		}
		# convert rule to regular expression
		$hmm_primary_rules{$key}{rule} =~ s/ AND /$regexAND/g ;
		$hmm_primary_rules{$key}{rule} = "(?=.*$hmm_primary_rules{$key}{rule}).*" ;		
	}
	return (\%hmm_primary_rules, \%pfams) ;
}

sub search_primary_PFAM_rules {
	# Search for AOIs on the basis of PFAM rules defined in table AOI_identification_rules.txt

	my ($hmm_primary_rules, $pfams) = read_table_class_properties() ;
	my %hmm_primary_rules = %$hmm_primary_rules ;
	my %pfams = %$pfams ;
	genomics::print_hash_of_hashes(%hmm_primary_rules);

	# 1. HMMsearch with all motifs
	my %hmmresult ;
	genomics::write_log($logfile,"Search primary motifs from primary PFAM rules table:",'true');
	foreach my $pfam ( keys %pfams ) {
		genomics::write_log($logfile, "Search for primary pfam $pfam with Evalue=$pfams{$pfam}",'true');
		my $tblout = "$sessiondir/$pfam.hmmprimary.result" ;
		unlink ($tblout) if (-e $tblout) ; # remove the old file
		$tmp = "$conf{hmmsearch} --noali --cpu $conf{hmmsearch_processors} --tblout $tblout -E $pfams{$pfam} $conf{hmm_models_dir}$pfam.hmm $sessiondir/PROT_$conf{orfs} >>$sessiondir/hmm_primary_search.log" ;
		system($tmp)  ;
		
		# Store the positions of the pfam domains
		my @lines = genomics::read_lines($tblout) ;
		chomp @lines ;
		foreach my $line (@lines) {
			if ($line =~ m/.*start=(\d+).*/) {	
				push @{$hmmresult{$pfam}{position}}, $1 ;
			}
		}
	}
	
	# 2. Check the rules
	my %pfamstrings ;
	foreach my $pfam (sort keys %hmmresult) {
		foreach my $pos (@{$hmmresult{$pfam}{position}}) {
			# check if other pfam is whitin range of this position
			my $pfamlist = $pfam ;
			foreach my $pfam2 (sort keys %hmmresult) {
				foreach my $pos2 (@{$hmmresult{$pfam2}{position}}) {
					$pfamlist .= " $pfam2" if ( ($pfam2 ne $pfam) and ($pfamlist !~ m/.*$pfam2.*/) and (abs($pos-$pos2) < $conf{contextsize} )) ;
				}
			}
			$pfamstrings{$pfamlist}{position} = $pos ;
			$pfamstrings{$pfamlist}{name} = $pfam ;
		}	
	}	
		
	# 3. Check the AOI using the primary pfam rules
	genomics::write_log($logfile, "\nCheck the AOI using the primary pfam rules:", 'true') ;
	my %results ;
	foreach my $pfamlist (sort  keys %pfamstrings) {
		genomics::write_log($logfile, "PFAMS wihtin one AOI at pos $pfamstrings{$pfamlist}{position}: $pfamlist", 'true') ;
		foreach my $key_rules (sort { $hmm_primary_rules{$a}{order} <=> $hmm_primary_rules{$b}{order} } keys %hmm_primary_rules) {
			if ($pfamlist =~ m/$hmm_primary_rules{$key_rules}{rule}/ ) {
				genomics::write_log($logfile, "\t$pfamlist fits $key_rules based on rule $hmm_primary_rules{$key_rules}{rule}", 'true');
				my $key = $pfamlist ;
				$results{$key}{pfamlist} = $pfamlist ;
				$results{$key}{position} = $pfamstrings{$pfamlist}{position} ;
				$results{$key}{name} = $hmm_primary_rules{$key_rules}{name} ;
			}
		}	
	}

	# 4. show results
	genomics::write_log($logfile, "\nThe AOIs identified:", 'true') ;
	foreach my $key (sort  keys %results) {
		genomics::write_log($logfile, "AOI for $results{$key}{name} identified on the basis of [$results{$key}{pfamlist}] domains at position $results{$key}{position}", 'true');
	}	
	
	return %results ;
}



sub orf_calling_prodigal {
	my $tmp = $conf{prodigal}."< $sessiondir/AOI.scaffolds.fna > $sessiondir/AOI_prodigal";
	system($tmp) ;

	my @lines = genomics::read_lines("$sessiondir/AOI_prodigal");
	chomp @lines ;
	
	my %tmp ;
	my $key = 'AOI';
	my $count = 0 ;
	open(FAA, ">$sessiondir/AOI.$conf{orfcalling_prog}.faa" ) or die("could not find the file $sessiondir/AOI.$conf{orfcalling_prog}.faa\n") ;
	foreach my $line (@lines) {
		# lines are from prodigal result
		if ($line =~ m/(AOI_.+)\ start=(\d+)\ end=(\d+)\"/g) {
			$key = $1 ;
		} elsif ($line =~ m/CDS\s+(complement\(|)(\d+)\.\.(\d+)/g) {
			$count++;
			my $orf = "$key;ORF_$count" ; 
			$tmp{$key}{$orf}{orf} 		= $orf ; 
			$tmp{$key}{$orf}{strand} 	= "+" ;
			$tmp{$key}{$orf}{strand} 	= "-" if ($1 eq "complement(") ;
			$tmp{$key}{$orf}{start} 	= $2 ;
			$tmp{$key}{$orf}{end} 		= $3 ;			
			$tmp{$key}{$orf}{dna} = substr $AOI{$key}{dna}, ($tmp{$key}{$orf}{start}-1),  (1+$tmp{$key}{$orf}{end}-$tmp{$key}{$orf}{start}) ;
			$tmp{$key}{$orf}{dna} = genomics::inverse_complement($tmp{$key}{$orf}{dna}) if ($tmp{$key}{$orf}{strand} eq "-") ;
			$tmp{$key}{$orf}{dna_len} = length ($tmp{$key}{$orf}{dna}) ;
			$tmp{$key}{$orf}{protein} = genomics::translate($tmp{$key}{$orf}{dna}) ;
			$tmp{$key}{$orf}{protein} =~ s/\-//g ;
			$tmp{$key}{$orf}{protein_len} = length ($tmp{$key}{$orf}{protein}) ;
			print FAA ">$tmp{$key}{$orf}{orf} strand$tmp{$key}{$orf}{strand} start=$tmp{$key}{$orf}{start}\n$tmp{$key}{$orf}{protein}\n";	# print proteins
		}
	}	
	return %tmp ;
}

sub orf_calling_glimmer {
	# ORF prediction using GLIMMER
	my ($dnafile, $dna_out, $prot_out, $RBS) = @_ ;
	my $tmp;	
	# depending on the presence of the trainings model glimmer.icm build the model or else execute glimmer3 using the model generated by bagel3_wrapper
	if ($glimmer_train ) {  # glimmer.icm present
		$tmp = $conf{glimmer}."/glimmer3  -o50 -g60 -t30 --linear $dnafile $sessiondir/glimmer.icm $sessiondir/glimmer";
	} else { # glimmer.icm should be made
		$tmp = "$conf{bagel3dir}/glimmer3_script.sh $dnafile $sessiondir/glimmer $conf{glimmer}";
	}
	print "GLIMMER: $tmp\n" ;
	system($tmp) ;

	# make one string of each dna sequence in the scaffold file
	my @lines = genomics::read_lines($dnafile);
	chomp @lines ;
	my %dna ;
	my $AOI = 'AOI_1';
	foreach my $line (@lines) {
		if ($line =~ m/^\>(.*) start.*/g) {
			$AOI = $1 ;
		} else {	
			$dna{$AOI} .= $line ;
		}	
	}	


	# read the glimmer results file
	my $glimmerfile = "$sessiondir/glimmer.predict" ;
	@lines = genomics::read_lines("$glimmerfile");
	chomp @lines ;
	
	my %tmp ;	
	open(FFN, ">$sessiondir/$dna_out" )  or die("could not find the file $sessiondir/$dna_out\n") ;
	open(FAA, ">$sessiondir/$prot_out" ) or die("could not find the file $sessiondir/$prot_out\n") ;
	foreach my $line (@lines) {
		# lines are from glimmer.predict
		if ($line =~ m/^\>(.*)\ start=(\d+)\ end=(\d+)/g) {
		#if ($line =~ m/^\>(.*)(\ |\,|\.).*/g) {
			$AOI = $1 ;
		} else {	
			my @col = split /\s+/, $line ;
			my $key = $AOI.";".$col[0] ; 
			$key =~ s/orf00/orf/ ;
			$tmp{$key}{orf} = $key ; 
			$tmp{$key}{orfcall} = 'glimmer' ; 
			$tmp{$key}{strand} = $col[3] ; 
			$tmp{$key}{rbs} = '';
			if ($tmp{$key}{strand} > 0) {
				$tmp{$key}{start} 	= $col[1] ; 	
				$tmp{$key}{end} 	= $col[2] ;
				if ($RBS eq 'rbs=on') { 
					$tmp{$key}{rbs} = substr $dna{$AOI}, ($tmp{$key}{start}-16), 18 ;
					if ($tmp{$key}{rbs} =~ m/(..GG.{5,14}(ATG|GTG|TTG))/g) { $tmp{$key}{rbs} = $1 ; }  # get the RBS
				}
			} else {
				$tmp{$key}{start} 	= $col[2] ; 	
				$tmp{$key}{end} 	= $col[1] ;
				if ($RBS eq 'rbs=on') { 
					$tmp{$key}{rbs} = genomics::inverse_complement( substr $dna{$AOI},($tmp{$key}{end}-3),19 ) ;
					if ($tmp{$key}{rbs} =~ m/(..GG.{5,14}(ATG|GTG|TTG))/g) { $tmp{$key}{rbs} = $1 ; }  # get the RBS
				}
			}			
			$tmp{$key}{dna} = substr $dna{$AOI}, ($tmp{$key}{start}-1),  (1+$tmp{$key}{end}-$tmp{$key}{start}) ;
			$tmp{$key}{dna} = genomics::inverse_complement($tmp{$key}{dna}) if ($tmp{$key}{strand} < 0) ;
			$tmp{$key}{dna_len} = length ($tmp{$key}{dna} ) ; 
			$tmp{$key}{protein} = genomics::translate($tmp{$key}{dna}) ;
			$tmp{$key}{protein} =~ s/\-//g ;
			$tmp{$key}{protein_len} = length ($tmp{$key}{protein}) ;
			print FFN ">$tmp{$key}{orf} strand$tmp{$key}{strand} start=$tmp{$key}{start} end=$tmp{$key}{end} RBS=$tmp{$key}{rbs} \n$tmp{$key}{dna}\n";	# print genes
			print FAA ">$tmp{$key}{orf} strand$tmp{$key}{strand} start=$tmp{$key}{start} end=$tmp{$key}{end} RBS=$tmp{$key}{rbs} \n$tmp{$key}{protein}\n";	# print proteins
		}
	}
	close FFN ;
	close FAA ;
	return %tmp ;
}

sub get_class_from_blast {
	# take the headers of the database: e.g. >Astexin |Asticcacaulis excentricus |Lasso peptide |
	my ($dbname, $subject) = @_ ;
	my %db ;
	$subject =~ s/\s+//g ; 
	my @lines = genomics::read_lines($dbname) ;
	chomp @lines ;
	foreach my $line (@lines) {
		if ($line =~ m/\>(.*)/) {  # the fasta header 
			my @items = split "|", $1 ;
			my $key = $items[0] ;
			$key =~ s/\s+//g ;
			$db{$key} = $items[2] ;
		}	
	}
	my $result = '';
	$result = $db{$subject} if (defined($db{$subject})) ;
	return $result ;
}

sub Select_AOI {
	# This function determine the Area of Interest and returns the AOI as a hash table
	my %starts ;
	my %class ;
	my $AOI_identifier = '';	
	foreach my $key (keys %AOI_identified_primaryPFAMrules) {
		my $start = $AOI_identified_primaryPFAMrules{$key}{position};
		$starts{$start} = 'primary_pfams' ;
		$class{$start} = $AOI_identified_primaryPFAMrules{$key}{name} ;
		&genomics::write_log($logfile,"Primary PFAM $key Start=$start Classname=$class{$start}", 'true') ;
		$AOI_identifier = 'primary_pfams' ;
	}

	# use the blast results of bacteriocinI bacteriocinII and bacteriocinIII to add AOIs
	# foreach my $key (keys %blast_bacteriocinI) { 
		# push @starts, $all_orfs{$key}{start} ;
		# &genomics::write_log($logfile,"Blast hit to bacteriocinI db $key Start=$all_orfs{$key}{start}", 'true') ;
	# }	

	
	
	foreach my $key (keys %blast_bacteriocinI) { 
		my $to_close = 0 ;
		if ($AOI_identifier eq 'primary_pfams') { # do not add AOI if to close to primary pfam
			foreach my $start (sort {$a<=>$b} keys %starts) {
				$to_close = 1 if (abs($all_orfs{$key}{start}-$start) < $conf{contextsize}) ;
			}
		}
		if (!$to_close) {
			$starts{$all_orfs{$key}{start}} = 'blast_bacteriocinI'  ;
			$class{$all_orfs{$key}{start}} = get_class_from_blast("$conf{bagel3dir}/$conf{bacteriocinI_db}", $blast_bacteriocinI{$key}{Subject}) ;
			$class{$all_orfs{$key}{start}} = 'unmodified' if ($class{$all_orfs{$key}{start}} eq '') ;

			&genomics::write_log($logfile,"Blast hit to bacteriocinI db $key Start=$all_orfs{$key}{start} class=$class{$all_orfs{$key}{start}}", 'true') ;
			
		}	
	}
	
	foreach my $key (keys %blast_bacteriocinII) { 
		my $to_close = 0 ;
		if ($AOI_identifier eq 'primary_pfams') { # do not add AOI if to close to primary pfam
			foreach my $start (sort {$a<=>$b} keys %starts) {
				$to_close = 1 if (abs($all_orfs{$key}{start}-$start) < $conf{contextsize}) ;
			}
		}
		if (!$to_close) {
			$starts{$all_orfs{$key}{start}} = 'blast_bacteriocinII'  ;
			$class{$all_orfs{$key}{start}} = get_class_from_blast("$conf{bagel3dir}/$conf{bacteriocinII_db}", $blast_bacteriocinII{$key}{Subject}) ;
			$class{$all_orfs{$key}{start}} = 'unmodified' if ($class{$all_orfs{$key}{start}} eq '') ;

			&genomics::write_log($logfile,"Blast hit to bacteriocinII db $key Start=$all_orfs{$key}{start} class=$class{$all_orfs{$key}{start}}", 'true') ;
			
		}	
	}

	foreach my $key (keys %blast_bacteriocinIII) { 
		my $to_close = 0 ;
		if ($AOI_identifier eq 'primary_pfams') { # do not add AOI if to close to primary pfam
			foreach my $start (sort {$a<=>$b} keys %starts) {
				$to_close = 1 if (abs($all_orfs{$key}{start}-$start) < $conf{contextsize}) ;
			}
		}
		if (!$to_close) {
			$starts{$all_orfs{$key}{start}} = 'blast_bacteriocinIII'  ;
			$class{$all_orfs{$key}{start}} = get_class_from_blast("$conf{bagel3dir}/$conf{bacteriocinIII_db}", $blast_bacteriocinIII{$key}{Subject}) ;
			$class{$all_orfs{$key}{start}} = 'bacteriocin >10kd' if ($class{$all_orfs{$key}{start}} eq '') ;

			&genomics::write_log($logfile,"Blast hit to bacteriocinIII db $key Start=$all_orfs{$key}{start}", 'true') ;
		}	
	}	
	
	my $count = 0 ;
	if (keys %starts > 0 ) {	# al least one start point is found
		my $key  ;
		my $current_start ;  
		my $current_end ;
		my $last_start ='';
		my $first_key = 1 ;
		foreach my $start (sort {$a<=>$b} keys %starts) { 
			if ($first_key) {
				$first_key = 0 ;
				# calculate the AOI start and end positions
				$current_start 	= $start - $conf{contextsize};  
				$current_start 	= 1 if ($current_start < 1 ) ;
				$current_end 	= $start + $conf{contextsize} ;
				$current_end 	= length($query_DNA) if ($current_end > length($query_DNA)) ;				
			} elsif ( $start-$current_end > $conf{contextsize} ) {  # in this case the new start is in a new AOI
				# new AOI found, time to save the previous start
				$count++ ;
				$key = "AOI_$count" ;
				$AOI{$key}{gene_start}	= $last_start ;
				$AOI{$key}{start}		= $current_start ;
				$AOI{$key}{end} 		= $current_end ;
				$AOI{$key}{dna} 		= substr $query_DNA, $current_start, $current_end-$current_start ;
				$AOI{$key}{aoi_select} 	= $starts{$last_start} ;
				$AOI{$key}{class} 		= $class{$last_start} ;
				# calculate the new position of the AOI
				$current_start 	= $start - $conf{contextsize};  
				$current_start 	= 0 if ($current_start < 0 ) ;
				$current_end 	= $start + $conf{contextsize};
				$current_end 	= length($query_DNA) if ($current_end > length($query_DNA)) ;
			} 
			$last_start = $start ;
		}
		# the first or last AOI:
		$count++ ;
		$key = "AOI_$count" ;
		$AOI{$key}{gene_start}	= $last_start ;
		$AOI{$key}{start}		= $current_start ;
		$AOI{$key}{end} 		= $current_end ;
		$AOI{$key}{dna} 		= substr $query_DNA, $current_start, $current_end-$current_start;
		$AOI{$key}{aoi_select} 	= $starts{$last_start} ;
		$AOI{$key}{class} 		= $class{$last_start} ;
		
		# write all AOIs to a fasta file
		open (FILE, "> $sessiondir/AOI.scaffolds.fna") or die ("Could not write $sessiondir/AOI.scaffolds.fna");
		foreach my $key (sort keys %AOI) {
			&genomics::write_log($logfile,"Area Of Interrest [$key] of class '$AOI{$key}{class}' identified on the basis of $AOI{$key}{aoi_select}: ".$key." start=".$AOI{$key}{start}." end=".$AOI{$key}{end}." based on gene start=".$AOI{$key}{gene_start},'true');
			print FILE ">$key start=".$AOI{$key}{start}." end=".$AOI{$key}{end}."\n";
			print FILE $AOI{$key}{dna}."\n" ;
		}	
		close FILE ;
	}
	return $count ;
}

sub parseparam {
    # reads the command line parameters
	my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        $sessiondir 	= shift(@arg) if($var eq '-s') ;
		$queryfile		= shift(@arg) if($var eq '-i') ;
		$reference_genome	= shift(@arg) if($var eq '-r') ;
		$glimmer_train	= shift(@arg) if($var eq '-glimmer') ;
		$make_graphics	= shift(@arg) if($var eq '-graphics') ;
        $outputfile 	= shift(@arg) if($var eq '-o') ;
    }
    die "No filename found\n$usage" if (!$queryfile) ;
	#remove the last / from sessiondir to make it universal
	$sessiondir =~ s/\/$// ;
}

sub read_conf {
	# returns the configuration file as a hash table
	my $filename = shift ;
	my %tmp ;
	open(FILE,"$filename") or die('could not find configuration file $filename\n') ;
	my(@lines) = <FILE>;
	chomp @lines ;
	foreach my $line (@lines) {
		$line =~ s/(\s|\t)//g ; # remove spaces and tabs
		if ($line =~ m/(.*)\=(.*)/g ) {	$tmp{$1}=$2; }	
	}
	return %tmp ;
}

sub read_DNA {
	# read plain DNA	
	my $filename = shift ;
	open (FILE,'<', "$filename") or die ("$filename does not exists\n") ;
	my(@lines) = <FILE>;					# Get the tablefile into a array
	chomp @lines;
	close FILE ;
	my $header ;
	my $DNA ;
	foreach my $line (@lines)  { 
		$line =~ s/(\n|\r|\x0d)//g;  # remove DOS/windows line breaks etc
		if ($line =~m/^\>/) {
			$header = $line ;
		} else {
			$DNA .= uc($line) ;
		}
	}
	return $DNA ;
}
