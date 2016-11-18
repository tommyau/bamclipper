#!/usr/bin/env perl
# inject separator after all alignment of a read pair
use strict;
use warnings;
my $readname;
while (<>) {
    my $original_line = $_;
    if (substr($original_line, 0, 1) ne "@") {
        my @fields = split("\t", $original_line);
	print "__\n" if defined $readname && $readname ne $fields[0];
        $readname = $fields[0];
    }
    print $original_line;
}