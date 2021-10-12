#!/usr/bin/env perl
use strict;
use Getopt::Long;
my $ploidy=2;
my $min_coverage=10;
my $min_hap_coverage=5;
my $max_hap_dist = 10;
my $min_minor_count = 5;
my $min_maf = 0.2;
GetOptions(
    "ploidy=i" => \$ploidy,
    "min_coverage=i" => \$min_coverage,
    "min_hap_coverage=i" => \$min_hap_coverage,
    "max_hap_dist=i" => \$max_hap_dist,
    "min_minor_count=i" => \$min_minor_count,
    "min_maf=n" => \$min_maf,
);
my $mpileup = shift;
my $reads_fasta = shift;
my $base_name = shift;

my %readname2seq;
open(F, $reads_fasta) || die $!;
{
    local $/="\n>";
    while (<F>) {
        chomp;
        s/^>//;
        my ($readname, $seq) = /^(\S+)[^\n]*\n(.*)/s;
        $readname2seq{$readname} = $seq;
    }
}
close F;

#keep track in case we don't end up building subsets
my %all_reads;
my %reads2haps;
my $var_count = 0;
open(M, $mpileup) || die $!;
while (<M>) {
    chomp;
    my @data = split /\t/;
    my @reads = split /,/, $data[$#data];
    next if scalar(@reads) < $min_coverage;
    map {$all_reads{$_} = 1} @reads;
    my @bp = split //, $data[4];
    #for now, let's only worry about SNPs
    if (scalar(@bp) != scalar(@reads)) {
        print STDERR "skipping non-SNP\n";
        next;
    }
    my %allele_counts;
    map {$allele_counts{$_}++;} @bp;
    my @alleles = sort {$allele_counts{$b} <=> $allele_counts{$a}} keys %allele_counts;
    #in theory, we should be able to handle these but their tendency to be frequent errors unlinked to the haplotype is causing some trouble in the current implementation, so skipping them as well
    if ($alleles[0] eq "*" || $alleles[1] eq "*") {
        print STDERR "skipping 1bp del\n";
        next;
    }
    if (!($allele_counts{$alleles[1]}/scalar(@bp) >= $min_maf || $allele_counts{$alleles[1]} >= $min_minor_count)) {
        print STDERR "skipping SNP below MAF and minor count\n";
        next;
    }
    for (my $i = 0; $i < @reads; $i++) {
        if (!defined $reads2haps{$reads[$i]}) {
            $reads2haps{$reads[$i]} = [];
        }
        if ($bp[$i] eq $alleles[0]) {
            $reads2haps{$reads[$i]}->[$var_count] = "A";
        }
        elsif ($bp[$i] eq $alleles[1]) {
            $reads2haps{$reads[$i]}->[$var_count] = "B";
        }
        else {
            $reads2haps{$reads[$i]}->[$var_count] = "N";
        }
    }
    $var_count++;
}
close M;
print STDERR "phaseable variant count: $var_count\n";

my %hap2count;
my %hap2reads;
foreach my $read (keys %reads2haps) {
    for (my $i=0; $i < $var_count; $i++) {
        if (!length($reads2haps{$read}->[$i])) {
            $reads2haps{$read}->[$i] = "N";
        }
    }
    my $hap = join("", @{$reads2haps{$read}});;
    if ($hap !~ /N/) {
        $hap2count{$hap}++;
    }
    $hap2reads{$hap}->{$read} = 1;
}
my @haps = sort {$hap2count{$b} <=> $hap2count{$a}} keys %hap2count;
my %haps_to_use;
#in this case, assume a single haplotype, shared by all reads
if (scalar(@haps) == 0) {
    $haps_to_use{hom} = "hom";
    $hap2reads{hom} = \%all_reads;
}
else {
    if (scalar(@haps) == 1 || $hap2count{$haps[1]} < $min_hap_coverage) {
        $haps_to_use{hom} = "hom";
        $haps_to_use{hom} = $haps[0];
    }
    else {
        for (my $i = 1; $i <= $ploidy && $i <= scalar(@haps); $i++) {
            #we need to take into account the fact that we may have fewer valid haplotypes than ploidy
            $haps_to_use{"het".$i} = $haps[$i-1];
        }
    }
    #for better quantification, assign reads not exactly matching one of the determined haps to their closest match; would it be better to do haplotype correction on a site-by-site basis? should we check for recombinant haplotypes (e.g. by looking for haps that were not sufficiently close to either form)
    my %haps_to_use_lookup = map {$_ => 1;} values %haps_to_use;
    my @haps_to_use = keys %haps_to_use_lookup;
    foreach my $read (keys %reads2haps) {
        my $hap = join("",@{$reads2haps{$read}});
        #if it isn't already associated with a well-supported hap, find the closest to assign it to
        if (! $haps_to_use_lookup{$hap}) {
            my @matched_haps;
            my $min_dist = ~0;
            for (my $i=0; $i < scalar @haps_to_use; $i++) {
                my $dist = &haplotype_distance($hap, $haps_to_use[$i]); 
                if ($dist < $min_dist && $dist <= $max_hap_dist) {
                    $min_dist = $dist;
                    @matched_haps = ($haps_to_use[$i]);
                }
                elsif ($dist == $min_dist && $dist <= $max_hap_dist) {
                    push(@matched_haps,$haps_to_use[$i]);
                }
            }
            if (scalar(@matched_haps) == 1) {
                $hap2reads{$matched_haps[0]}->{$read} = 1
            }
            else {
                print STDERR "$read looks like a recombinant!?\n";
            }
        }
    }
}
foreach my $hap_key (keys %haps_to_use) {
    my $hap_val = $haps_to_use{$hap_key};
    my $hap_count = scalar(keys %{$hap2reads{$hap_val}}); 
    next unless $hap_count >= $min_hap_coverage;
    open(H, ">$base_name.$hap_key.read_count-$hap_count.fna") || die $!;
    foreach my $readname (keys %{$hap2reads{$hap_val}}) {
        print H ">$readname\n$readname2seq{$readname}\n";
    }
    close H;
    #TODO: maybe this should be outside the script to enable different downstream applications of the basic premise of the script (ie to partition into readsets by haplotype)
    system("mafft $base_name.$hap_key.read_count-$hap_count.fna > $base_name.$hap_key.read_count-$hap_count.aln");
    system("cons.pl $base_name.$hap_key.read_count-$hap_count $base_name.$hap_key.read_count-$hap_count.aln | make_faidxable.pl > $base_name.$hap_key.read_count-$hap_count.cons.fna");
    system("getorf -find 1 $base_name.$hap_key.read_count-$hap_count.cons.fna stdout | getorf_longest.pl > $base_name.$hap_key.read_count-$hap_count.cons.faa");
}

#stolen with modification from report_recombination_events.pl in mpr_hacks
sub convert_to_binary() {
        my ($genotypes) = @_;
        my @genotypes = split //, $genotypes;
        my $retval = '';
        for (my $i = 0; $i < @genotypes; $i++) {
            if ($genotypes[$i] eq "A") {
                vec($retval,$i,2) = 0b00;
            }
            elsif ($genotypes[$i] eq "B") {
                vec($retval,$i,2) = 0b11;
            }
            else {
                vec($retval,$i,2) = 0b01;
            }
        }
        return $retval;
}

sub convert_from_binary() {
    my ($genotypes) = @_;
    my @genotype_codes = ("A", "N", "N", "B");
    my @bits = split(//, unpack("b*", $genotypes));
    my $retval = "";
    for (my $i = 0; $i < @bits; $i += 4) {
        $retval .= $genotype_codes[oct("0b".join("", $bits[$i], $bits[$i+1]))];
    }
    return $retval;
}

sub haplotype_distance() {
    my ($g1, $g2) = @_;
    #my $bits = unpack("b*", ($g1^$g2));
    #my $distance =()= $bits =~ /1/g;
    my $distance = unpack("%32b*", ($g1^$g2));
    return $distance;
}
