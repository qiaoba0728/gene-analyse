#!/usr/bin/perl
##--- zwy 2016.9.27 for web --##

use strict;
use Getopt::Std;
our($opt_f,$opt_b,$opt_a,$opt_1,$opt_2,$opt_r,$opt_s,$opt_3,$opt_4,$opt_5,$opt_6);
getopts('f:b:a:1:2:r:s:3:4:5:6:');

my $reference    = $opt_f;
my $bowtie_index = $opt_b;
my $lable1       = (defined $opt_1)?$opt_1:"S1";
my $lable2       = (defined $opt_2)?$opt_2:"S2";
my $threads      = (defined $opt_a)?$opt_a:6;
my $step         = (defined $opt_s)?$opt_s:1000000;
my $range        = (defined $opt_r)?$opt_r:3000000;
my $fq1_1	 = $opt_3;
my $fq1_2	 = $opt_4;
my $fq2_1	 = $opt_5;
my $fq2_2	 = $opt_6;

my $Usage = "\n$0 -f <Reference genome> -b <Index> -3 -4 -5 -6 <fastq files>
Input (required):
	-f <FILE>  reference genome file
	-b <FILE>  bowtie2 index [same as reference name]
	-3 <FILE>  pool1 left
	-4 <FILE>  pool1_right
	-5 <FILE>  pool2_left
	-6 <FILE>  pool2_right
	
Advanced options:	
	-a <INT>   number of threads [6]
	-s <INT>   window step [1000000]
	-r <INT>   window range [3000000]
	-1 <STR>   pool1's name [S1]
        -2 <STR>   pool2's name [S2]
	
	\n";
die $Usage unless ($opt_f && $opt_b && $opt_3 && $opt_4 && $opt_5 && $opt_6);
################################################################
my $time1 = time();
my $vcf = "$lable1\_$lable2";
my $fai  = "$reference.fai";

# ----------- index the genome using samtools ---------- #
unless (-e $fai){
	print "  indexing the genome with samtools\n";
	!system "samtools faidx $reference" or die "Something wrong with samtools faidx:$!";
}else{
	print "$fai -- found! indexing skipped!\n";
	!system "samtools faidx $reference" or die "Something wrong with samtools faidx:$!";
}

## check chromosome numbers ##
open CHRNUM,'<',"$reference" or die "Cannot open reference";
my $chrseq_numbers = 0;
while (<CHRNUM>) {
	chomp;
	if (/>/) {
		$chrseq_numbers ++;
	}
	if ($chrseq_numbers == 100) {
		last;
	}
}
close CHRNUM;


# -------- running snpMapper --------------#
!system "snpMapper-1.07_forDNABSA.pl -f $reference -b $bowtie_index -3 $fq1_1 -4 $fq1_2 -5 $fq2_1 -6 $fq2_2 -a $threads -d 0 -o /data/output/smsnpMapper_out > /data/log/snpMapper.log" or die "Something wrong with snpMapper-107_forDNABSA.pl:$!";
!system "delta.pl /data/output/smsnpMapper_out/smcandidatesnps_D15d0.txt 0.4 /data/output/smsnpMapper_out/smcandidatesnps_D15d0.4.txt" or die "Something wrong wich delta.pl!";
!system "snp_to_vcf.pl -i /data/output/smsnpMapper_out/smcandidatesnps_D15d0.txt -o /data/output/$vcf.vcf" or die "Somerthing wrong with snp2vcf:$! $vcf";
!system "Enzyme_ParaFly.pl -v /data/output/$vcf.vcf -s $reference -p $threads -a > /data/log/Enzyme.log" or die "Something wrong with Enzyme_parafly.pl!";
mkdir "Final_Results"; 

if ($chrseq_numbers <= 30) {
	!system "delta_distribution.pl -i /data/output/smsnpMapper_out/smcandidatesnps_D15d0.txt -o /data/output/delta_distribution.txt -r $range -s $step > /data/log/delta_plot.log" or die "Cannot run delta_distribution.pl:$!";
	print "\nstart plot:\n";
	!system "Rscript BSA_permutation_parallel.R /data/output/smsnpMapper_out/smcandidatesnps_D15d0.txt $range $step $threads" or die "ERROR with R_BSA_permutation.R:$!";
}

#system "cp smsnpMapper_out/smcandidatesnps_D15d0.txt Final_Results/";
#system "cp smsnpMapper_out/smcandidatesnps_D15d0.4.txt Final_Results/";
#system "mv Better_enzyme_site.txt /data/output/";
#system "mv Enzyme_site.txt /data/output/";

my $time2 = time();
my $time = sprintf "%.1f",($time2-$time1)/3600;
$time .= 'h';
print "Pipeline finished: $time\n";


