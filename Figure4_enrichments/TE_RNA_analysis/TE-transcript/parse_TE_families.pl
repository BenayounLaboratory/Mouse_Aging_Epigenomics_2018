#! /usr/bin/perl

use warnings;
use strict;

#2017-09-12
# only keep TE with known phylogeny in TEtrascript (":")


unless (scalar @ARGV == 1) {
	die "parse_TE_families.pl TE_names_list.txt\n";
}

my $TEs = shift @ARGV;

open (TENAMES,$TEs) or die "$TEs fiel doesn't exist: $!\n";


print "Rowname\tName\tFamily\tClass\n";


while (my $line = <TENAMES>) {
	
	# read and skip unless has :
	next unless ($line =~ m/:/);
	#print $line;
	
	chomp $line;
    
    my @linedata = split(/:/, $line);
    
    my $num = scalar @linedata;
	#print $num; # is 3
	# name/family/class
	
	print $line."\t".join("\t",@linedata)."\n";
	
}


close TENAMES;

exit;