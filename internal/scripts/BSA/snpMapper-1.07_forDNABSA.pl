#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Path qw(make_path remove_tree);

our ($opt_z,$opt_o,$opt_D,$opt_d,$opt_c,$opt_f,$opt_3,$opt_4,$opt_5,$opt_6,$opt_a,$opt_b,$opt_q,$opt_Q,$opt_s,$opt_n);
getopt ("z:o:D:d:c:f:3:4:5:6:a:b:q:Q:s:n:");
########### snpMapper.pl Qun Hu and Weiyi Zhang
# Function: Get base counts from samtools' mpileup output
# the mpileup file should be produced with mutant.bam as first inputfile, wild_type.bam as second input file
# samtools command should be run as:
#     samtools mpile -ABf ref.fa mut.bam wt.bam >mut-wt_raw.mpileup
#
# 21:47 2014/3/24 v0.0.1
# 19:47 2014/3/26 v1.0.1 modify the program structure
# 19:42 2014/3/30 v1.0.3 Incorporate positions that have 'skipped reference' meaning
# 16:05 2014/4/9  v1.04
# 22:41 2014/5/21 v1.05
# 22:18 2016/9/27 v1.07-BOWTIE2 modify the tophat2 to bowtie2 and this version is for web
#
# Future version should add algorithm to consider contamination in each pool -- but future will not come
#
#############################################################################
# this program takes mpileup file as input, and outoput several files regarding the SNP
# It can be modified to subject tetraploid species RAN-seq analysis

my $time1 = time();
my $version = '1.07-bowtie2';
my $Usage = "Usage:\n  $0 -OPTIONS VALUES

options:
--input options
        -f REFGENOME    reference genome for mapping
        -b bowtie2index  bowtie2 index

--output options
        -o OUTDIR   output direcotry containg all output files [snpMappr_out]
        -3 FILE     bool1 left reads path
        -4 FILE     pool1 right reads path
        -5 FILE     pool2 left reads path
        -6 FILE     pool2 right reads path
        -p FILE     mpileup format file produced by samtools mpileup command, if this option is not provided, then -c,-m,-f have to be provided
        -s STR      prefix string [sm]

--criteria options
        -D INT      minimum read coverage to be considered for valid SNP [15]
        -n FLOAT    pos whose wildtype allele frequency gt FLOAT will be ignored [0.95]
        -d FLOAT    minimum delta snp index, required [0.2]
        -q INT      skip alignments with mapQ smaller than INT [20]
        -Q INT      skip bases with baseQ/BAQ smaller than INT [20]

--perfromanat options
        -a INT      cpus cores used for the analysis [6]
        \n";

die $Usage unless ($opt_f && $opt_b && $opt_3 && $opt_4 && $opt_5 && $opt_6);

# -------------- Start snpMapper --------------------
print "\n\nStarting $0 ...\n";
my $refgenome           = $opt_f; # reference genome
my $bowtie2_index       = $opt_b || 'bowtie2-index'; # intermediate outfile for bowtie2 index
my $fq1_left            = $opt_3;
my $fq1_right           = $opt_4;
my $fq2_left            = $opt_5;
my $fq2_right           = $opt_6;
my $user_defined_depth  = (defined $opt_D)?$opt_D:15;  # read coverage for a given position, default=15;
my $delta            = (defined $opt_d)?$opt_d:0; # delta of SNP-index of locus between mutant and wild type
my $nopoly           = (defined $opt_n)?$opt_n:0.8;
my $cpu              = (defined $opt_a)?$opt_a:6;
my $mapq             = (defined $opt_q)?$opt_q:40;
my $baseQ            = (defined $opt_Q)?$opt_Q:30;
my $prefix           = (defined $opt_s)?$opt_s:'sm';
my $outdir           = (defined $opt_o)?$opt_o:"${prefix}snpMapper_out"; # outfile
my $SNP_Locus_File   = "$outdir/${prefix}snpresults.snps";
my $Putative_SNP     = "$outdir/${prefix}candidatesnps_D${user_defined_depth}d$delta.txt";
my $merge_reads      = (defined $opt_z)?$opt_z:'F';
my $snpdepth         = 5;

###################### snpMapper.pl Main program starts from here ##################################

print "  Checking fastq files ...</br>";

make_path($outdir,{verbose=>1,mode=>0777}) unless (-e $outdir);
my $tophat_outdir1 = $prefix.'bowtie2_out1';
my $tophat_outdir2 = $prefix.'bowtie2_out2';
make_path($tophat_outdir1, {-verbose=>1,mode=>0777}) unless -e $tophat_outdir1 && -d $tophat_outdir1;
make_path($tophat_outdir2, {-verbose=>1,mode=>0777}) unless -e $tophat_outdir2 && -d $tophat_outdir2;


my @bowtie_index_files = ("$bowtie2_index.1.bt2","$bowtie2_index.2.bt2",
                                                  "$bowtie2_index.3.bt2","$bowtie2_index.4.bt2",
                                                  "$bowtie2_index.rev.1.bt2",
                                                  "$bowtie2_index.rev.2.bt2");


# checking whether genome has be indexed
my $flag_bowtie2 = 1;
foreach my $file (@bowtie_index_files){
        unless (-e $file){
        $flag_bowtie2 = 0;
        last;
        }
}
if ($flag_bowtie2 == 0){
        print "$bowtie2_index  -- Not Found! Creating...<.br>";
        &command_mkbowtie_index($refgenome,$bowtie2_index);
        print "done\n\n";
}else{
        print "$bowtie2_index -- Found! bowtie2 index building skipped!</br>";
}

my ($mutbamfile, $wtbamfile) = ("NA","NA");
        print "  mutbamfile -- not found\n  Runing bowtie2 for sample mutant ...";
        &command_tophat($bowtie2_index,$fq1_left,$fq1_right,$tophat_outdir1,$cpu);
        print "finished!\n\n";
        $mutbamfile = "$tophat_outdir1/acc.sorted.bam";
#=cut 1

        print "  wtbamfile -- not found\n  Runing bowtie2 for sample wildtype ...";
        &command_tophat($bowtie2_index,$fq2_left,$fq2_right,$tophat_outdir2,$cpu);
        print "finished!</br></br>";
        $wtbamfile = "$tophat_outdir2/acc.sorted.bam";
#=cut 2
# ----  check whether all required files are present ----- #
check_file($mutbamfile,$wtbamfile,$refgenome);
#exit(0);

print "  -- Starting converting pileup format to an elite SNP file for simplicity ... </br>";
open  my $fh_pileup, "samtools mpileup -q $mapq -Q $baseQ -ABf $refgenome $mutbamfile $wtbamfile|"  or die "Can not read input from samtools!\n";

my $header = "#chr      pos     refbase mut.cov mut.A   mut.C   mut.G   mut.T   wt.cov  wt.A    wt.C    wt.G    wt.T\n";
open my $fh_out1, ">$SNP_Locus_File" or die "Can not open $SNP_Locus_File\n";    # .snps
open my $fh_out3, ">$Putative_SNP"   or die "Can not open $Putative_SNP\n";    # .delta
print $fh_out1 $header;

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
                if ( $mutcov1 < $snpdepth || $wtcov1 < $snpdepth ) {
                        next;
                }
                my @alleles =  sort { $mut{$b} <=> $mut{$a} or $wt{$b} <=> $wt{$a} } keys %mut;
                my $mut_freq = sprintf "%.4f", $mut{ $alleles[0] } / $mutcov1;    # 0.7
                my $wt_freq  = sprintf "%.4f", $wt{ $alleles[0] } / $wtcov1;      # 0.3
                my $deltaSNP = sprintf "%.4f", abs( $mut_freq - $wt_freq );

                if ($mut_freq > $nopoly && $wt_freq > $nopoly) {
                        next;
                }

        # delta SNP-index great than user-defined value will be considered as putative SNPs
                print $fh_out3 "$Chr\t$Pos\t$refBase\t$mut_freq\t$wt_freq\t$deltaSNP\t",
        "$mutcov1\t$mut{'A'}\t$mut{'C'}\t$mut{'G'}\t$mut{'T'}\t$wtcov1\t$wt{'A'}\t$wt{'C'}\t$wt{'G'}\t$wt{'T'}\n";
}


close $fh_pileup;
close $fh_out1;
close $fh_out3;

my $time2 = time();
my $time  = ( $time2 - $time1 ) / 60;
print "$0 finished! Total time elapsed: $time</br>";
print "At this points, please find your candidate SNPs in $outdir! Wish you a good Luck!</br>";
exit(0);

############################ subroutine definition #######################
sub check_file{
  my @files = @_;

        my $flag = 1;
        foreach my $file (@files){
                next if $file =~ /^#/;
                next if $file =~ /^\s*$/;
                if (-e $file){
                        print "  $file -- found</br>";
                }else{
                        print "  $file -- not found</br>";
                        $flag = 0;
                }
        }
        if ($flag == 0){
                print ("some files does not exist,see above!</br>");
                exit(1);
        }
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

sub command_mkbowtie_index{
        my ($refgenome,$bowtie2_index) = @_;
        system "bowtie2-build $refgenome $bowtie2_index 2> bowtie_build.log";
}

sub command_tophat{
        my ($bowtie_index,$s_left,$s_right,$outdir,$cpu)  = @_;
        !system "bowtie2 -k 1 --quiet --no-unal --omit-sec-seq -p $cpu -x $bowtie_index -1 $s_left -2 $s_right -S $outdir/acc.sam 2> $outdir/bowtie2_aln.out" or die "Error in bowtie2:$!";
        !system "samtools view -bhS $outdir/acc.sam > $outdir/acc.bam" or die "Error in samtools view (sam to bam)";
        system "samtools sort -@ $cpu -o $outdir/acc.sorted.bam $outdir/acc.bam";
        system "samtools index $outdir/acc.sorted.bam";
        unlink "$outdir/acc.sam"; unlink "$outdir/acc.bam";
#       ~/local/bin/tophat2 -p 10  -o Tophat_HQ6 -i 50 -I 30000 /disk1/genome/lactuca-bowtie2 Sample_HQ6/HQ6_GGCTAC_L007_R1.fq Sample_HQ6/HQ6_GGCTAC_L007_R2.fq
}