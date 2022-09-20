#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    $0 genome.fasta > telomere_info.txt

    --split-length <int>    default: 100000
    --overlap-length <int>    default: 10000
   
    --repeat-unit <string>    default: TTAGGG
    please find target repeat unit in（http://telomerase.asu.edu/sequences_telomere.html
    vertebrate sp.      TTAGGG
    plants sp.          TTTAGGG
    Pezizomycotina      TTAGGG

    --min-repeat-num <int>    default: 4

USAGE
if (@ARGV==0){die $usage}

my ($splitLength, $overlapLength, $repeatunit, $minRepeatNum);
GetOptions(
    "split-length:i" => \$splitLength,
    "overlap-length:i" => \$overlapLength,
    "repeat-unit:s" => \$repeatunit,
    "min-repeat-num:s" => \$minRepeatNum,
);
$splitLength ||= 100000;
$overlapLength ||= 10000;
$repeatunit ||= "CCCTAAA"; #Here we tested all possible repeat units CCCTAAA/TTTAGGG, CCCCTAA/TTAGGGG, CCCTAA/TTTAGG, CCCTAA/TTAGGG and CCTAA/TTAGG
$repeatunit = uc($repeatunit);
my $repeatunit_reverse = reverse $repeatunit;
$repeatunit_reverse =~ tr/ATCG/TAGC/;
$minRepeatNum ||= 4;

# read genome sequences
open IN, $ARGV[0] or die "Can not open file $ARGV[0], $!\n";
my (%seq, $seq_id);
while (<IN>) {
    chomp;
    if (m/^>(\S+)/) {
        $seq_id = $1;
    }
    else {
        $_ = uc($_);
        $seq{$seq_id} .= $_;
    }
}
close IN;

# break genome sequences to several parts
my (%seq_split, %seq_length);
foreach my $id (keys %seq) {
    my $seq = $seq{$id};
    my $length = length($seq);
    $seq_length{$id} = $length;
    my $pos = 0;
    while ($pos < $length) {
        $seq_split{$id}{$pos} = substr($seq, $pos, $splitLength + $overlapLength);
        $pos += $splitLength;
    }
}

print "SeqID\tSeqLength\tStart\tEnd\tLength\tType\n";
foreach my $id (sort keys %seq_split) {
    foreach my $pos (sort {$a <=> $b} keys %{$seq_split{$id}}) {
        my $seq = $seq_split{$id}{$pos};
        while ($seq =~ m/(($repeatunit){$minRepeatNum,})/g) {
            my $length = length($1);
            my $end = pos($seq);
            $end = $end + $pos;
            my $start = $end - $length + 1;
            print "$id\t$seq_length{$id}\t$start\t$end\t$length\t$repeatunit\n";
        }
        while ($seq =~ m/(($repeatunit_reverse){$minRepeatNum,})/g) {
            my $length = length($1);
            my $end = pos($seq);
            $end = $end + $pos;
            my $start = $end - $length + 1;
            print "$id\t$seq_length{$id}\t$start\t$end\t$length\t$repeatunit_reverse\n";
        }
    }
}
