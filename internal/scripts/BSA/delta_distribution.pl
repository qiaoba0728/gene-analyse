#!/usr/bin/perl

use strict;
use Getopt::Std;
################################################################################################
our ($opt_i,$opt_o,$opt_s,$opt_r);
getopts('i:o:s:r:');

my $input = $opt_i;
my $output = $opt_o;
my $step = (defined $opt_s)?$opt_s:1000000;
my $range = (defined $opt_r)?$opt_r:3000000;

my $Usage = "\n$0 <INPUT> <option>
		-i <FILE> input file format : chromosome  snp_locus  delta
		-o <FILE> output file
		-s <INT>  the window step
		-r <INT>  the window range
		\n";

die $Usage unless ($opt_i && $opt_o);
################################################################################################
open IN,'<',"$input" or die "Cant open input file:$!";
open OUT,'>',"$output" or die "Cant create output file:$!";
my ($chromosome,$prechr,$locus,$delta,$counter,@array,@array2);
while (<IN>) {
	chomp;
	if (/\Amut/) { next };
	my ($chr11,$pos11,$s,$f,$scov,$fcov) = (split)[0,1,3,4,6,11];
	my $strue = $s * $scov + $f * $fcov;
	my $ftrue = (1-$s) * $scov + (1-$f) * $fcov;
	if ( ($s > 0.9 && $f > 0.9) || $strue < 4 || $ftrue < 4) { next };
	#if ($s > 0.8 && $f > 0.8) { next };
	#print "$chr11\t$pos11\t$s\t$f\t$scov\t$fcov\t$ftrue\t$strue\n";
	($chromosome,$locus,$delta) = (split)[0,1,5];
	if ($counter == 0) {
		push @array,"$chromosome\t$locus\t$delta\t$s\t$f";
		$counter = 1;	
		$prechr = $chromosome;
		next;
	}
	if ($chromosome eq $prechr) {
		push @array,"$chromosome\t$locus\t$delta\t$s\t$f";
		$prechr = $chromosome;
	}else{
		@array2 = &step(@array);
		foreach (@array2) {
			print OUT "$_\n";
			if ($_ eq @array2[-1]) {
				print OUT "\n";
			}
		}
		@array = ();
		@array2 = ();
		push @array,"$chromosome\t$locus\t$delta";
		$prechr = $chromosome;
	}
}

@array2 = &step(@array);
foreach (@array2) {
	print OUT "$_\n";
	if ($_ eq @array2[-1]) {
		print OUT "\n";
	}
}
close OUT;
################################################
open FORMAT,'<',"$output" or die;
open ROUT,'>',"R_plot.input";
while (<FORMAT>) {
	chomp;
	my ($chr,$pos,$delta) = (split)[0,1,2];
	if ($chr =~ /.*?(\d+)\z/) {
		my $num = $1;
		if ($num == 00 || $num == 0) {
			next;
		}
		print ROUT "$num\t$pos\t$delta\n";
	}
}
close FORMAT;
close ROUT;
#!system "/results/others/bin/R_BSA_plot.pl R_plot.input" or die "Something wrong with R:$!";
#####################--- sub progrem ---######################
sub step {
	my ($chr,$loc,$del,$start,$end,$sum,$count,$avr,$pos,@return,$outrange,$last,$s,$f);
	@return = ();
	$start = 0;
	$end = $range;
	$last = (split /\t/,@_[-1])[1];
	while (1) {	
		foreach (@_) {
			($chr,$loc,$del,$s,$f) = (split)[0,1,2,3,4];
			if ( $last < $end ) { 
				$outrange = 1;				
				last;
			}
			if ($loc > $start && $loc < $end) {
				if ($s > 0.9 && $f > 0.9) {		###filer snp
					next;
				}
				$sum += $del;
				$count ++;
			}elsif ($loc > $end) {
				if ($count == 0) {
					$avr = 0;
					$pos = ($start + $end) / 2000000;
					push @return,"$chr\t$pos\t$avr";
					$sum = ();
					$count = ();
					$start += $step;
					$end += $step;
					last;
				}else{
					$avr = $sum / $count;
					$pos = ($start + $end) / 2000000;
					push @return,"$chr\t$pos\t$avr";
					$sum = ();
					$count = ();
					$start += $step;
					$end += $step;
					last;
				}
			}
		}
		if ($outrange == 1) {
			$outrange = 0;
			last;
		}
	}
	$sum = ();
	$count = ();
	return @return;
}
