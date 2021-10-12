#!/usr/bin/env perl
#a hack, but for some reason much faster than the EMBOSS equivalent on large alignments
use strict;
my $cons_name = shift;

my @seqs;
my $aln_len;
$/="\n>";
while (<>) {
    chomp;
    my ($header, $seq) = /([^\n]*)\n(.*)/s;
    $seq =~ s/\n//g;
    #we don't want to treat lack of sequence info on ends as we do internal gaps
    my ($starting_dashes, $ending_dashes) = ($seq =~ /^(-*).*?(-*)$/); 
    my $starting_ns = 'N'x length($starting_dashes);
    my $ending_ns = 'N'x length($ending_dashes);
    $seq =~ s/^-*/$starting_ns/;
    $seq =~ s/-*$/$ending_ns/;
    my @seq = split //, $seq;
    $aln_len = length($seq);
    push @seqs, \@seq;
}
my $cons;
for (my $i=0; $i < $aln_len; $i++) {
    my %counts;
    foreach my $seq (@seqs) {
        $counts{$seq->[$i]}++;
    }
    my (@cons_bp) = sort {$counts{$b} <=> $counts{$a}} keys %counts;
    my $cons_bp = ($cons_bp[0] eq "N" ? $cons_bp[1] : $cons_bp[0]);
    $cons .= uc($cons_bp) unless $cons_bp eq "-";
}
print ">$cons_name\n$cons\n";
