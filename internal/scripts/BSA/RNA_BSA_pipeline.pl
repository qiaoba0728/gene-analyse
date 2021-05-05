#!/usr/bin/perl
##--- zwy 2016.9.27 for web --##

use strict;
use Getopt::Std;
our($opt_f,$opt_b,$opt_a,$opt_1,$opt_2,$opt_r,$opt_s,$opt_3,$opt_4,$opt_5,$opt_6,$opt_g,$opt_d);
getopts('f:b:a:1:2:r:s:3:4:5:6:g:d:');

my $reference    = $opt_f;
my $bowtie_index = $opt_b;
my $lable1       = (defined $opt_1)?$opt_1:"S1";
my $lable2       = (defined $opt_2)?$opt_2:"S2";
my $threads      = (defined $opt_a)?$opt_a:6;
my $step         = (defined $opt_s)?$opt_s:1000000;
my $range        = (defined $opt_r)?$opt_r:3000000;
my $fq1_1        = $opt_3;
my $fq1_2        = $opt_4;
my $fq2_1        = $opt_5;
my $fq2_2        = $opt_6;
my $gff_file     = (defined $opt_g)?$opt_g:0;
my $delta_cut    = (defined $opt_d)?$opt_d:0.4;

my $Usage = "\n$0 -f <Reference genome> -b <Index> -3 <pool1 left reads> -4 <pool1 right reads> -5 <pool2 left reads> -6 <pool2 right reads>
Input (required):
        -f <FILE>  reference genome file
        -b <FILE>  hisat2 index [same as reference name]
        -3 <FILE>  pool1 left
        -4 <FILE>  pool1_right
        -5 <FILE>  pool2_left
        -6 <FILE>  pool2_right

Advanced options:
        -a <INT>   number of threads [6]
        -s <INT>   window step [1000000]
        -r <INT>   window range [3000000]
        -d <FLOAT> delta value cutoff [0.4]
        -1 <STR>   pool1's name [S1]
        -2 <STR>   pool2's name [S2]
        -g <FILE>  GTF file
        \n";
die $Usage unless ($opt_f && $opt_b && $opt_3 && $opt_4 && $opt_5 && $opt_6);
################################################################
my $time1 = time();
my $vcf = "$lable1\_$lable2";
my $fai  = "$reference.fai";

# ----------- index the genome using samtools ---------- #
unless (-e $fai){
        print "indexing the genome with samtools</br>";
        !system "samtools faidx $reference" or die "Something wrong with fasta file!";
}else{
        print "$fai -- found! indexing skipped!</br>";
        !system "samtools faidx $reference" or die "Something wrong with fasta file!";
}
############################################################## check chr numbers #
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
##
# -------- running snpMapper --------------#
if ($chrseq_numbers < 100) {
        !system "snpMapper-1.08.pl -f $reference -b $bowtie_index -3 $fq1_1 -4 $fq1_2 -5 $fq2_1 -6 $fq2_2 -a $threads -d 0 -o /data/output/smsnpMapper_out > /data/log/snpMapper.out" or die "Something wrong with snpMapper!";
}else{
        !system "snpMapper-1.07.pl -f $reference -b $bowtie_index -3 $fq1_1 -4 $fq1_2 -5 $fq2_1 -6 $fq2_2 -a $threads -d 0 -o /data/output/smsnpMapper_out > /data/log/snpMapper.out" or die "Something wrong with snpMapper!";
}
if ($delta_cut == 0.4) {
        !system "delta.pl /data/output/smsnpMapper_out/smcandidatesnps_D10d0.txt 0.4 /data/output/smsnpMapper_out/smcandidatesnps_D10d0.4.txt" or die "Something wrong wich delta.pl!";
}else{
        !system "delta.pl /data/output/smsnpMapper_out/smcandidatesnps_D10d0.txt $delta_cut /data/output/smsnpMapper_out/smcandidatesnps_D10d$delta_cut.txt" or die "Something wrong wich delta.pl!";
        !system "delta.pl /data/output/smsnpMapper_out/smcandidatesnps_D10d0.txt 0.4 /data/output/smsnpMapper_out/smcandidatesnps_D10d0.4.txt" or die "Something wrong with delta.pl!";
}

!system "snp_to_vcf.pl -i /data/output/smsnpMapper_out/smcandidatesnps_D10d0.txt -o /data/output/$vcf.vcf" or die "Somerthing wrong with snp2vcf:$!";
!system "Enzyme_ParaFly.pl -v /data/output/$vcf.vcf -s $reference -p $threads -a > /data/log/enzyme.log" or die "Something wrong with Enzyme_parafly.pl!";
#mkdir "Final_Results";

if ($chrseq_numbers < 100) {
        !system "delta_distribution.pl -i /data/output/smsnpMapper_out/smcandidatesnps_D10d0.txt -o /data/output/delta_distribution.txt -r $range -s $step > /data/log/delta.log" or die "Cannot run delta_distribution.pl:$!";
        print "\nstart plot:\n";
        !system "Rscript BSA_permutation_parallel.R /data/output/smsnpMapper_out/smcandidatesnps_D10d0.txt $range $step $threads" or die "ERROR with R_BSA_permutation.R:$!";
}

# expression #
if ($gff_file == 0) {
        !system "Stringtie.pl $threads $lable1 $lable2 $reference" or warn "ERROR with Stringtie:$!";
        $gff_file = "/data/output/Stringtie_out/merge.gtf";
}
my $bam1 = "/data/output/smsnpMapper_out/smhisat2_out1/acc.sorted.bam";
my $bam2 = "/data/output/smsnpMapper_out/smhisat2_out2/acc.sorted.bam";
!system "gfold.pl $bam1 $bam2 $gff_file" or warn "ERROR with GFOLD, please check GTF format:$!";

# move some results to Final_Results
#system "cp smsnpMapper_out/smcandidatesnps_D10d0.txt Final_Results/";
#system "cp smsnpMapper_out/smcandidatesnps_D10d0.4.txt Final_Results/";
##--- get Enzyme site of which delta value greater than 0.4---##
open BETTER,'<',"Better_enzyme_site.txt" or die;
my (%better_seq,%better_name,$id);
while (<BETTER>) {
        chomp;
        if (/>/) {
                $id = (split />/,(split)[0])[1];
                $better_seq{$id} = '';
                my ($chr,$pos) = (split /:/,$id)[0,2];
                my $newname = "$chr\t$pos";
                $better_name{$newname} = $id;
        }else{
                $better_seq{$id} .= $_;
        }
}
close BETTER;

open ENZYME,'<',"Enzyme_site.txt" or die;
my (%enzyme_seq,%enzyme_name,$id);
while (<ENZYME>) {
        chomp;
        if (/>/) {
                $id = (split />/,(split)[0])[1];
                $enzyme_seq{$id} = '';
                my ($chr,$pos) = (split /:/,$id)[0,2];
                my $newname = "$chr\t$pos";
                $enzyme_name{$newname} = $id;
        }else{
                $enzyme_seq{$id} .= $_;
        }
}
close ENZYME;

open DELTA04,'<',"/data/output/smsnpMapper_out/smcandidatesnps_D10d$delta_cut.txt" or die "Cannot open delta file:$!";
open BETTEROUT,'>',"Better_enzyme_site_$delta_cut.txt";
open ENZYMEOUT,'>',"Enzyme_site_$delta_cut.txt";
while (<DELTA04>) {
        chomp;
        my ($chr0,$pos0,$pf1,$pf2,$pd) = (split)[0,1,3,4,5];
        my $newname = "$chr0\t$pos0";
        if (exists $better_name{$newname}) {
                print BETTEROUT ">$better_name{$newname}\t$pf1\t$pf2\t$pd\n$better_seq{$better_name{$newname}}\n";
        }elsif (exists $enzyme_name{$newname}) {
                print ENZYMEOUT ">$enzyme_name{$newname}\t$pf1\t$pf2\t$pd\n$enzyme_seq{$enzyme_name{$newname}}\n";
        }
}
close DELTA04;
close BETTEROUT;
close ENZYMEOUT;

##---                          END                         ---##
my $time2 = time();
my $time = sprintf "%.1f",(($time2-$time1)/3600);
$time = "$time\h";
print "Pipeline finished! duration time: $time\n"