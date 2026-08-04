#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# --------------- Parse command-line arguments ---------------
my $bam;
GetOptions('bam=s' => \$bam) or die "Usage: $0 --bam INPUT.bam > output.bed\n";
die "Usage: $0 --bam INPUT.bam > output.bed\n" unless $bam;

# --------------- Open BAM stream ---------------
open my $fh, '-|', 'samtools', 'view', $bam
  or die "Cannot open BAM '$bam': $!\n";

my $count = 0;
while (<$fh>) {
    $count++;
    print STDERR "- processed $count reads...\n" if $count % 1_000_000 == 0;

    chomp;
    my @col = split "\t";
    my $flag = $col[1];

    # --------------- Filter out unmapped reads only ---------------
    next if $flag & (1 << 2);

    # --------------- Extract BC tag (or default to NA) ---------------
    my ($bc) = $_ =~ /BC:Z:(\S+)/;
    $bc //= 'NA';

    # --------------- Parse CIGAR and sum only M/D/N operations ---------------
    my $cigar_ref = parseCIGAR($col[5]);
    my %cig      = %{ $cigar_ref };
    my $consumed = 0;
    while (my ($act, $len) = each %cig) {
        my (undef, $op) = split /_/, $act;
        $consumed += $len if $op eq 'M' or $op eq 'D' or $op eq 'N';
    }

    # --------------- Convert to 0-based coordinates & apply Tn5 shift (+4/-5) ---------------
    my $start0 = $col[3] - 1;
    my $end0   = $start0 + $consumed;
    my ($cut, $strand);
    if ($flag & (1 << 4)) {
        $strand = '-';
        $cut    = $end0 - 5;
    } else {
        $strand = '+';
        $cut    = $start0 + 4;
    }
    next if $cut < 0;

    # --------------- Output 5 columns: chr, start, end, BC, strand ---------------
    my $chr = $col[2];
    my $beg = $cut;
    my $fin = $cut + 1;
    print join("\t", $chr, $beg, $fin, $bc, $strand), "\n";
}
close $fh;

# ----------------- Subroutine: parse CIGAR -----------------
sub parseCIGAR {
    my ($cigar) = @_;
    my @parts = split(/(?<=\d)(?=\D)|(?<=\D)(?=\d)/, $cigar);
    my %map;
    my $cnt = 0;
    for (my $i = 0; $i < @parts; $i += 2) {
        my $len = $parts[$i];
        my $op  = $parts[$i+1];
        $cnt++;
        $map{"${cnt}_$op"} = $len;
    }
    return \%map;
}

