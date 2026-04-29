#!/usr/bin/env perl
use strict;
use warnings;

my $outdir = shift @ARGV or die "usage: make_synthetic.pl OUTDIR\n";
my @samples = map { "S$_" } 1..8;

sub header {
    my ($fh) = @_;
    print $fh "##fileformat=VCFv4.3\n";
    print $fh "##contig=<ID=chr22,length=50818468>\n";
    print $fh "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    print $fh "##FORMAT=<ID=AB,Number=1,Type=Float,Description=\"Allele balance\">\n";
    print $fh "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths\">\n";
    print $fh "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n";
    print $fh "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n";
    print $fh "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n";
    print $fh "##FORMAT=<ID=PGT,Number=1,Type=String,Description=\"Physical phase GT\">\n";
    print $fh "##FORMAT=<ID=PID,Number=1,Type=String,Description=\"Physical phase ID\">\n";
    print $fh "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype quality\">\n";
    print $fh "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum depth\">\n";
    print $fh "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Strand bias\">\n";
    print $fh "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihoods\">\n";
    print $fh "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Sample filter\">\n";
    print $fh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", join("\t", @samples), "\n";
}

sub open_vcf {
    my ($name) = @_;
    open my $fh, ">", "$outdir/$name.vcf" or die "$outdir/$name.vcf: $!\n";
    header($fh);
    return $fh;
}

sub genotype {
    my ($i, $s, $n_alt) = @_;
    return "./." if (($i + $s) % 29) == 0;
    return "0/0" if (($i + $s) % 5) == 0;
    return "1|1" if $n_alt == 1 && (($i + $s) % 11) == 0;
    return $n_alt > 1 && (($i + $s) % 7) == 0 ? "1/2" : "0/1";
}

sub ad {
    my ($i, $s, $n_allele) = @_;
    return "." if (($i + $s) % 37) == 0;
    return join(",", map { (($i * 3 + $s * 5 + $_ * 7) % 40) } 0..($n_allele - 1));
}

sub pl {
    my ($i, $s, $n_allele) = @_;
    return "." if (($i + $s) % 41) == 0;
    my $n = $n_allele * ($n_allele + 1) / 2;
    return join(",", map { (($i + $s + $_) * 13) % 500 } 0..($n - 1));
}

sub gl {
    my ($i, $s, $n_allele) = @_;
    return "." if (($i + $s) % 31) == 0;
    my $n = $n_allele * ($n_allele + 1) / 2;
    return join(",", map { sprintf("%.2f", -1 * ((($i + $s + $_) % 20) / 3.0)) } 0..($n - 1));
}

my $fh = open_vcf("synthetic_ccdg_likelihood");
for my $i (1..2000) {
    my $pos = 20000000 + $i;
    my $phase = $i % 2 == 0;
    my $fmt = $phase ? "GT:AB:AD:DP:GQ:PGT:PID:PL" : "GT:AB:AD:DP:GQ:PL";
    my @vals;
    for my $s (0..$#samples) {
        my $gt = genotype($i, $s, 1);
        my $ab = $gt eq "0/1" ? sprintf("%.2f", (($i + $s) % 90) / 100) : ".";
        my $base = join(":", $gt, $ab, ad($i, $s, 2), (($i+$s)%80), (($i+$s)%99));
        if ($phase) {
            push @vals, join(":", $base, ($gt =~ /\|/ ? $gt : "0|1"), "${pos}_A_T", pl($i, $s, 2));
        } else {
            push @vals, join(":", $base, pl($i, $s, 2));
        }
    }
    print $fh join("\t", "chr22", $pos, ".", "A", "T", 50, "PASS", ".", $fmt, @vals), "\n";
}
close $fh;

$fh = open_vcf("synthetic_reordered_likelihood");
for my $i (1..2000) {
    my $pos = 20100000 + $i;
    my @vals;
    for my $s (0..$#samples) {
        push @vals, join(":", (($i+$s)%80), (($i+$s)%99), genotype($i, $s, 1), ad($i, $s, 2), pl($i, $s, 2));
    }
    print $fh join("\t", "chr22", $pos, ".", "G", "C", 50, "PASS", ".", "DP:GQ:GT:AD:PL", @vals), "\n";
}
close $fh;

$fh = open_vcf("synthetic_fixed_numeric");
for my $i (1..2000) {
    my $pos = 20200000 + $i;
    my @vals;
    for my $s (0..$#samples) {
        my $hq = (($i+$s)%150) . "," . (($i+$s+9)%150);
        my $sb = join(",", map { ($i + $s + $_) % 30 } 0..3);
        push @vals, join(":", genotype($i, $s, 1), $hq, (($i+$s)%60), $sb);
    }
    print $fh join("\t", "chr22", $pos, ".", "C", "A", 50, "PASS", ".", "GT:HQ:MIN_DP:SB", @vals), "\n";
}
close $fh;

$fh = open_vcf("synthetic_float_string");
for my $i (1..2000) {
    my $pos = 20300000 + $i;
    my @vals;
    for my $s (0..$#samples) {
        my $ft = (($i+$s)%13) == 0 ? "LowQual" : "PASS";
        push @vals, join(":", genotype($i, $s, 2), gl($i, $s, 3), $ft, (($i+$s)%80), (($i+$s)%99));
    }
    print $fh join("\t", "chr22", $pos, ".", "A", "C,G", 50, "PASS", ".", "GT:GL:FT:DP:GQ", @vals), "\n";
}
close $fh;

$fh = open_vcf("synthetic_multiallelic_likelihood");
for my $i (1..1200) {
    my $pos = 20400000 + $i;
    my $n_alt = 1 + ($i % 3);
    my @alts = qw(C G T);
    my $alt = join(",", @alts[0..($n_alt - 1)]);
    my @vals;
    for my $s (0..$#samples) {
        push @vals, join(":", genotype($i, $s, $n_alt), ad($i, $s, $n_alt + 1), (($i+$s)%90), (($i+$s)%99), pl($i, $s, $n_alt + 1));
    }
    print $fh join("\t", "chr22", $pos, ".", "A", $alt, 50, "PASS", ".", "GT:AD:DP:GQ:PL", @vals), "\n";
}
close $fh;
