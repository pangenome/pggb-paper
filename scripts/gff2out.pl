#!/usr/bin/perl

use strict;
use warnings;

while (my $line = <>) {
    chomp $line;
    if ($line =~ /^#/) {
        next;
    }
    my @fields = split "\t", $line;
    if (@fields != 9) {
        next;
    }
    my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = @fields;
    my ($id) = $attributes =~ /ID=([^;]+)/;
    my ($name) = $attributes =~ /Name=([^;]+)/;
    my ($classification) = $attributes =~ /Target=([^;]+)/;
    my ($identity) = $attributes =~ /Identity=([^;]+)/;

    my $length = $end - $start + 1;
    my $class_family = $name || $type;

    print "10000\t0.001\t0.001\t0.001\t$seqid\t$start\t$end\t($length)\t$strand\t$length\t$classification\t1\t50\t(99)\t+\n";
}
