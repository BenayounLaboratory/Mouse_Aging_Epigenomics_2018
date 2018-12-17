#! /usr/bin/perl

use warnings;
use strict;

# 2017-12-12
# parse HOMER maps into one file

unless (@ARGV >= 2) {
	print "parse_HOMER_mappings.pl <outfile name> <kallistofile1> ... <kallistofileN>";

}

my $outname = shift @ARGV;

# to get gene length on the first file
my $start = 1;

my %outresults = ();

# prepare header
my $header = "Name\tlength\t";
$header .= join("\t",@ARGV);
$header .= "\n";

# loop on files
foreach my $file (@ARGV) {
	
	open(FILE,$file) or die "Could not open $file: $!\n";
	
	# skip first line
	my $headerline = <FILE>;
	
	if ($start == 1) {
		
		while (my $line = <FILE>) {
			#Name	Length	Counts
			#NM_026810	43380	29
			
			my @linedata = get_line_data($line);
			#print "@linedata\n:";

			push(@{$outresults{$linedata[0]}}, ($linedata[1],$linedata[2]) );
		}
		
		$start = 0;
	
		
	} else {
	
		while (my $line = <FILE>) {
			my @linedata = get_line_data($line);
			push(@{$outresults{$linedata[0]}}, $linedata[2] );
		}
		
	}
	
	close FILE;
}



#### output results
open(OUT,'>',$outname) or die "Could not open $outname: $!\n";

print OUT $header;

foreach my $gene (sort keys %outresults) {
	
	#print $gene."\n";
	
	my $outline = $gene."\t".join("\t",@{$outresults{$gene}})."\n";
	print OUT $outline;
}


close OUT;

exit;

###########################################################
# SUBROUTINES
###########################################################

###########################################################
# a subroutine that separates fields from a data line and
# returns them in an array

sub get_line_data {

    my $line = $_[0];
    
    chomp $line;  

    my @linedata = split(/\t/, $line);
        
    return @linedata;
}
