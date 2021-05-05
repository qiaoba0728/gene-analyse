#!/usr/bin/perl

$inputfile = $ARGV[0];
$delta_cutoff = $ARGV[1];
$outputfile = $ARGV[2];
open OUT,'>',"$outputfile";
open IN ,'<',"$inputfile" or die "cant open inputfile:$! $inputfile";
while (<IN>) {
	chomp;
	$delta = (split)[5];
	if ($delta > $delta_cutoff) {
		print OUT "$_\n";
	}
}

