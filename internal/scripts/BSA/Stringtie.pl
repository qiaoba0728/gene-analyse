#!/usr/bin/perl

use strict;
my ($threads,$p1,$p2,$genome) = @ARGV[0,1,2,3];
my $Usage = "\n\t$0 <threads> <pname1> <pname2> <genome>
\n";
die $Usage unless (@ARGV == 4);
##############################################################################################
#------ From here: main scripts ---------#
print "Start stringtie\n";
print "$threads,$p1,$p2,$genome\n";
if (-e "/data/output/Stringtie_out") {
        system "rm -rf Stringtie_out";
}
mkdir "/data/output/Stringtie_out";
print "Maked directory Stringtie_out\n";
##############################################################################################
my $bam1 = "/data/output/smsnpMapper_out/smhisat2_out1/acc.sorted.bam";
my $bam2 = "/data/output/smsnpMapper_out/smhisat2_out2/acc.sorted.bam";
print "$bam1\n$bam2\n";
!system "stringtie $bam1 -p $threads -o /data/output/Stringtie_out/$p1 > /data/log/Stringtie1.log" or die "Wrong with stringtie:$!";
!system "stringtie $bam2 -p $threads -o /data/output/Stringtie_out/$p2 > /data/log/Stringtie2.log" or die "Wrong with stringtie:$!";
!system "stringtie --merge /data/output/Stringtie_out/$p1 /data/output/Stringtie_out/$p2 -o /data/output/Stringtie_out/merge.gtf > /data/log/Stringtie_merge.log" or die "Wrong with stringtie merge:$!";
my $GTF = "/data/output/Stringtie_out/merge.gtf";

