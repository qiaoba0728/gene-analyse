#!/usr/bin/perl

use strict;
my ($file,$depth) = @ARGV[0,1];
my $Usage = "\n\t$0 <mpileup file> <Depth>
\n";
die $Usage unless (@ARGV == 2);
#############################################
my $prefix = (split /\.part/,$file)[0];
my $out3 = "$prefix.txt";

open my $fh_out3,'>',"$out3";
&snpMapper($file);
close $fh_out3;

sub snpMapper {
	my $medfile_path = @_[0];
	open my $fh_pileup,'<',"$medfile_path" or die "Cannot open medfile:$! $medfile_path";
	
	while (<$fh_pileup>) {
		chomp;

	        next if /^\[/;
       	 	next if /^</;

		my @Contents = (split /\s+/,$_);

		# Omit any pos where only one sample mapped
		next if @Contents < 9;

		my ( $Chr, $Pos, $refBase ) = @Contents[ 0, 1, 2 ];

		# For convenience, Sample1 for mutant pool, Sample2 for wild pool
		my ( $mutcov0, $mutbases, $wtcov0, $wtbases ) = @Contents[ 3, 4, 6, 7 ];

		# Omit any intronic position?
		# next if ($Mut_Bases =~ />|</ or $Wt_Bases =~ />|</);

		# Omit low-coverage position

		# Calculate reads coverage for the position
		my @mutCounts = base_counter( $mutbases, $refBase );
		my @wtCounts  = base_counter( $wtbases,  $refBase );

		my %mut = (
			'A' => $mutCounts[0],
			'C' => $mutCounts[1],
			'G' => $mutCounts[2],
			'T' => $mutCounts[3]
		);
		$mut{$refBase} += $mutCounts[4];
		my %wt = (
			'A' => $wtCounts[0],
			'C' => $wtCounts[1],
			'G' => $wtCounts[2],
			'T' => $wtCounts[3]
		);
		$wt{$refBase} += $wtCounts[4];

		my $mutcov1 = $mut{'A'} + $mut{'C'} + $mut{'G'} + $mut{'T'};
		my $wtcov1  = $wt{'A'} + $wt{'C'} + $wt{'G'} + $wt{'T'};
		if ( $mutcov1 < $depth || $wtcov1 < $depth ) {
			next;
		}
		my @alleles =  sort { $mut{$b} <=> $mut{$a} or $wt{$b} <=> $wt{$a} } keys %mut;
		my $mut_freq = sprintf "%.4f", $mut{ $alleles[0] } / $mutcov1;    # 0.7
		my $wt_freq  = sprintf "%.4f", $wt{ $alleles[0] } / $wtcov1;      # 0.3
		my $deltaSNP = sprintf "%.4f", abs( $mut_freq - $wt_freq );
		
		if ($mut_freq > 0.95 && $wt_freq > 0.95) {
			next;
		}
	# delta SNP-index great than user-defined value will be considered as putative SNPs
		print $fh_out3 "$Chr\t$Pos\t$refBase\t$mut_freq\t$wt_freq\t$deltaSNP\t",
	"$mutcov1\t$mut{'A'}\t$mut{'C'}\t$mut{'G'}\t$mut{'T'}\t$wtcov1\t$wt{'A'}\t$wt{'C'}\t$wt{'G'}\t$wt{'T'}\n";
	}
	close $fh_pileup;
}

# base_counter: calculate base counts for each base in (A,C,G,T) order
sub base_counter {
	my ( $sample_bases, $refbase ) = @_;

	# Convert all dot and comma symbol to ref base
	#$sample_bases =~ s/\.|,/$refbase/gi;
	
	my $baseE = ($sample_bases =~ tr/.|,/./);
	# Remove patterns that represents INDELs
	while ( $sample_bases =~ /(.*?)[+-](\d+)[ATCG.,]+/ig ) {
		$sample_bases =~ s/(.*?)[+-](\d+)[ATCGNatcgn]{$2}(.*)/$1$3/i;
	}

	my $baseA = ($sample_bases =~ tr/A|a/A/);
	my $baseC = ($sample_bases =~ tr/C|c/C/);
	my $baseG = ($sample_bases =~ tr/G|g/G/);
	my $baseT = ($sample_bases =~ tr/T|t/T/);
	return ( $baseA, $baseC, $baseG, $baseT ,$baseE);
}

sub sum{
	my $sum = 0;
	for (@_){
		$sum += $_; 
	}
	return $sum;
}
