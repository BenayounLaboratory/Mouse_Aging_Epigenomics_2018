#! /usr/bin/perl

use warnings;
use strict;
use Math::Round;

# a script to create a unified count matrix for TE and genes

unless (scalar @ARGV >= 2) {
	die "usage: Parse_homer_counts.pl Repeat_file Gene_file \n";
}


# get file names
my $TEfile = shift @ARGV;
my $Genefile = shift @ARGV;

# matched pattern containing the sample
$TEfile =~ m/(.+)_RNA_repeats\.txt/;
my $sample = $1;

my $outfile = $sample."_gene_and_repeat_counts.txt";

# create output file
open(OUT, '>', $outfile) or die "Unable to create $outfile: $!\n";

my $header = "Name\tLength\tCounts\n";
print OUT $header;


## 1. parse genes
open(GENES, $Genefile) or die "Unable to open $Genefile: $!\n";

while (my $line = <GENES>) {
	next if ($. == 1);
	
	my @linedata = get_line_data($line);
	
	my $outline = $linedata[0]."\t".round($linedata[5])."\t".round($linedata[8])."\n";
	print OUT $outline;

}

close GENES;


## 2. parse repeats
open(TES, $TEfile) or die "Unable to open $TEfile: $!\n";

while (my $line = <TES>) {
	next if ($. == 1);
	
	my @linedata = get_line_data($line);
	
	my $outline = $linedata[0]."\t".round($linedata[5])."\t".round($linedata[8])."\n";
	print OUT $outline;

}

close TES;

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
