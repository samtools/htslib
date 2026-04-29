#!/usr/bin/env perl
use strict;
use warnings;

my $outdir = shift @ARGV or die "usage: make_large_synthetic.pl OUTDIR [NSAMPLES]\n";
my $nsamples = shift @ARGV || 2048;
my $scale = shift @ARGV || 1;
my @samples = map { "S$_" } 1..$nsamples;

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
    return "./." if (($i + $s) % 97) == 0;
    return "0/0" if (($i + $s) % 5) == 0;
    return "1|1" if $n_alt == 1 && (($i + $s) % 23) == 0;
    return $n_alt > 1 && (($i + $s) % 7) == 0 ? "1/2" : "0/1";
}

sub ad {
    my ($i, $s, $n_allele) = @_;
    return "." if (($i + $s) % 131) == 0;
    return join(",", map { (($i * 3 + $s * 5 + $_ * 7) % 120) } 0..($n_allele - 1));
}

sub pl {
    my ($i, $s, $n_allele) = @_;
    return "." if (($i + $s) % 137) == 0;
    my $n = $n_allele * ($n_allele + 1) / 2;
    return join(",", map { (($i + $s + $_) * 13) % 700 } 0..($n - 1));
}

sub gl {
    my ($i, $s, $n_allele) = @_;
    return "." if (($i + $s) % 127) == 0;
    my $n = $n_allele * ($n_allele + 1) / 2;
    return join(",", map { sprintf("%.2f", -1 * ((($i + $s + $_) % 30) / 4.0)) } 0..($n - 1));
}

sub write_ccdg_like {
    my ($name, $records) = @_;
    my $fh = open_vcf($name);
    for my $i (1..$records) {
        my $pos = 21000000 + $i;
        my $phase = $i % 2 == 0;
        my $fmt = $phase ? "GT:AB:AD:DP:GQ:PGT:PID:PL" : "GT:AB:AD:DP:GQ:PL";
        my @vals;
        for my $s (0..$#samples) {
            my $gt = genotype($i, $s, 1);
            my $ab = $gt eq "0/1" ? sprintf("%.2f", (($i + $s) % 90) / 100) : ".";
            my $base = join(":", $gt, $ab, ad($i, $s, 2), (($i+$s)%160), (($i+$s)%99));
            if ($phase) {
                push @vals, join(":", $base, ($gt =~ /\|/ ? $gt : "0|1"), "${pos}_A_T", pl($i, $s, 2));
            } else {
                push @vals, join(":", $base, pl($i, $s, 2));
            }
        }
        print $fh join("\t", "chr22", $pos, ".", "A", "T", 50, "PASS", ".", $fmt, @vals), "\n";
    }
    close $fh;
}

sub write_reordered {
    my ($name, $records) = @_;
    my $fh = open_vcf($name);
    for my $i (1..$records) {
        my $pos = 22000000 + $i;
        my @vals;
        for my $s (0..$#samples) {
            push @vals, join(":", (($i+$s)%160), (($i+$s)%99), genotype($i, $s, 1), ad($i, $s, 2), pl($i, $s, 2));
        }
        print $fh join("\t", "chr22", $pos, ".", "G", "C", 50, "PASS", ".", "DP:GQ:GT:AD:PL", @vals), "\n";
    }
    close $fh;
}

sub write_multiallelic {
    my ($name, $records) = @_;
    my $fh = open_vcf($name);
    for my $i (1..$records) {
        my $pos = 23000000 + $i;
        my $n_alt = 1 + ($i % 3);
        my @alts = qw(C G T);
        my $alt = join(",", @alts[0..($n_alt - 1)]);
        my @vals;
        for my $s (0..$#samples) {
            push @vals, join(":", genotype($i, $s, $n_alt), ad($i, $s, $n_alt + 1), (($i+$s)%160), (($i+$s)%99), pl($i, $s, $n_alt + 1));
        }
        print $fh join("\t", "chr22", $pos, ".", "A", $alt, 50, "PASS", ".", "GT:AD:DP:GQ:PL", @vals), "\n";
    }
    close $fh;
}

sub write_float_string {
    my ($name, $records) = @_;
    my $fh = open_vcf($name);
    for my $i (1..$records) {
        my $pos = 24000000 + $i;
        my @vals;
        for my $s (0..$#samples) {
            my $ft = (($i+$s)%17) == 0 ? "LowQual" : "PASS";
            push @vals, join(":", genotype($i, $s, 2), gl($i, $s, 3), $ft, (($i+$s)%160), (($i+$s)%99));
        }
        print $fh join("\t", "chr22", $pos, ".", "A", "C,G", 50, "PASS", ".", "GT:GL:FT:DP:GQ", @vals), "\n";
    }
    close $fh;
}

sub write_phase_width_variation {
    my ($name, $records) = @_;
    my $fh = open_vcf($name);
    for my $i (1..$records) {
        my $pos = 25000000 + $i;
        my @vals;
        for my $s (0..$#samples) {
            my $gt = genotype($i, $s, 1);
            my $pgt = (($i + $s) % 29) == 0 ? "." : ($gt =~ /\|/ ? $gt : "0|1");
            my $pid;
            if (($i + $s) % 31 == 0) {
                $pid = ".";
            } elsif (($i + $s) % 7 == 0) {
                $pid = "${pos}_${s}_A_T_LONG_PHASE_SET";
            } elsif (($i + $s) % 5 == 0) {
                $pid = "${pos}_A_T";
            } else {
                $pid = "P" . (($i + $s) % 97);
            }
            push @vals, join(":", $gt, ad($i, $s, 2), (($i+$s)%160),
                             (($i+$s)%99), $pgt, $pid, pl($i, $s, 2));
        }
        print $fh join("\t", "chr22", $pos, ".", "A", "T", 50, "PASS", ".",
                       "GT:AD:DP:GQ:PGT:PID:PL", @vals), "\n";
    }
    close $fh;
}

sub write_mixed_likelihood {
    my ($name, $records) = @_;
    my $fh = open_vcf($name);
    for my $i (1..$records) {
        my $pos = 26000000 + $i;
        my $n_alt = ($i % 17) == 0 ? 8 : (($i % 11) == 0 ? 2 : 1);
        my @alts = qw(C G T AA AC AG AT GA);
        my $alt = join(",", @alts[0..($n_alt - 1)]);
        my $n_allele = $n_alt + 1;
        my @vals;
        for my $s (0..$#samples) {
            my $ad = ad($i, $s, $n_allele);
            my $pl = pl($i, $s, $n_allele);
            if (($i % 19) == 0 && $ad ne ".") {
                my @ad = split /,/, $ad;
                pop @ad;
                $ad = join(",", @ad);
            }
            if (($i % 23) == 0 && $pl ne ".") {
                my @pl = split /,/, $pl;
                pop @pl;
                $pl = join(",", @pl);
            }
            push @vals, join(":", genotype($i, $s, $n_alt), $ad,
                             (($i+$s)%160), (($i+$s)%99), $pl);
        }
        print $fh join("\t", "chr22", $pos, ".", "A", $alt, 50, "PASS", ".",
                       "GT:AD:DP:GQ:PL", @vals), "\n";
    }
    close $fh;
}

sub write_gt_first_reordered {
    my ($name, $records) = @_;
    my $fh = open_vcf($name);
    for my $i (1..$records) {
        my $pos = 27000000 + $i;
        my @vals;
        for my $s (0..$#samples) {
            push @vals, join(":", genotype($i, $s, 1), (($i+$s)%160),
                             ad($i, $s, 2), (($i+$s)%99), pl($i, $s, 2));
        }
        print $fh join("\t", "chr22", $pos, ".", "G", "C", 50, "PASS", ".",
                       "GT:DP:AD:GQ:PL", @vals), "\n";
    }
    close $fh;
}

sub write_two_string_float {
    my ($name, $records) = @_;
    my $fh = open_vcf($name);
    for my $i (1..$records) {
        my $pos = 28000000 + $i;
        my @vals;
        for my $s (0..$#samples) {
            my $ft = (($i+$s)%17) == 0 ? "LowQual" : "PASS";
            my $pid = (($i+$s)%13) == 0 ? "." : "PS" . (($i * 11 + $s) % 100000);
            push @vals, join(":", genotype($i, $s, 2), $ft, $pid,
                             gl($i, $s, 3), (($i+$s)%160));
        }
        print $fh join("\t", "chr22", $pos, ".", "A", "C,G", 50, "PASS", ".",
                       "GT:FT:PID:GL:DP", @vals), "\n";
    }
    close $fh;
}

unless ($ENV{SYNTHETIC_ONLY_NEW}) {
    write_ccdg_like("large_ccdg_likelihood_${nsamples}s", 20000 * $scale);
    write_reordered("large_reordered_likelihood_${nsamples}s", 20000 * $scale);
    write_multiallelic("large_multiallelic_likelihood_${nsamples}s", 16000 * $scale);
    write_float_string("large_float_string_${nsamples}s", 16000 * $scale);
}
write_phase_width_variation("large_phase_width_variation_${nsamples}s", 12000 * $scale);
write_mixed_likelihood("large_mixed_likelihood_${nsamples}s", 12000 * $scale);
write_gt_first_reordered("large_gt_first_reordered_${nsamples}s", 12000 * $scale);
write_two_string_float("large_two_string_float_${nsamples}s", 12000 * $scale);
