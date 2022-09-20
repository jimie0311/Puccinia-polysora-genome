#!/opt/perl/bin/perl

use warnings;
use strict;
use Bio::SeqIO;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Bio::Index::Fasta;
use File::Find;
use Data::Dumper;

my $seqID=shift or die "sequence ID?";
my $fastafile=shift or die "fasta file?";

my $s_wrote1;
my $s_read1;
my $print_it1;
$s_read1 = $s_wrote1 = $print_it1 = 0; 

open(ID,"<$seqID"); 

#open(FASTA,">$seqID.fasta");
my %ids1=();
while (<ID>) { 
  s/\r?\n//; /^>?(\S+)/; $ids1{$1}++; 
} 
my $num_ids1 = keys %ids1; 

open(F, "<$fastafile");
while (<F>) { 
  if (/^>(\S+)/) { 
    $s_read1++; 
    if ($ids1{$1}) { 
      $s_wrote1++; 
      $print_it1 = 1; 
      $_ = join("\n", $_);
      #$_=join("\n", $new, $_);
      #$_=join($new, $_);
      delete $ids1{$1} 
    } else { 
      $print_it1 = 0  
    }
  }
  if ($print_it1) { 
     print $_;
#    print FASTA $_; 
  }
} 
warn "Searched $s_read1 FASTA records.\nFound $s_wrote1 IDs out of $num_ids1 in the ID list.\n";
