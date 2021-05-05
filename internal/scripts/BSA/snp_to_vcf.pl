#!/usr/bin/perl

##---2017-7-1 fix some refence bugs --##
use strict;
use Getopt::Std;
our($opt_i,$opt_o);
getopts('i:o:');

my $input = $opt_i;
my $out = $opt_o;

my $Usage = "\n$0 -i <snp_file> -o <out>
	-i <FILE> input snp file
	-o <FILE> output
	\n";
die $Usage unless ($opt_i && $opt_o);
######################################################################
open OUT,'>',"$out";
open IN,'<',"$input" or die "Cant open input file!:$!";
my ($chr,$pos1,$refbase,$cov1,$a1,$c1,$g1,$t1,$cov2,$a2,$c2,$g2,$t2);
while (<IN>) {
	chomp;
	if (/\A#/) { next };
	($chr,$pos1,$refbase,$cov1,$a1,$c1,$g1,$t1,$cov2,$a2,$c2,$g2,$t2) = (split)[0,1,2,6,7,8,9,10,11,12,13,14,15];
	my $aa = $a1+$a2;
	my $cc = $c1+$c2;
	my $gg = $g1+$g2;
	my $tt = $t1+$t2;
	my @SR1 = ("$aa\tA","$cc\tC","$gg\tG","$tt\tT");
	my @SR2 = sort {(split /\t/,$b)[0] <=> (split /\t/,$a)[0]} @SR1;
	my ($base1,$base2,$base3) = ( (split /\t/,$SR2[0])[1], (split /\t/,$SR2[1])[1], (split /\t/,$SR2[2])[1]);
	my ($base1_cov,$base2_cov,$base3_cov) = ( (split /\t/,$SR2[0])[0], (split /\t/,$SR2[1])[0], (split /\t/,$SR2[2])[0] );
	if ($base3_cov > 5) {
              	next;
        }else{
		if ($base1 ne $refbase && $base2 eq $refbase) {
			my $basex = $base1;
			my $basex_cov = $base1_cov;
			$base1 = $base2;
			$base1_cov = $base2_cov;
			$base2 = $basex;
			$base2_cov = $basex_cov;
		}
                print OUT "$chr\t$pos1\t.\t$base1\t$base2\t50\tDP4=$base1_cov,0,$base2_cov,0\n";
        }
}

