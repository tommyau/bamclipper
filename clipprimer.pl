#!/usr/bin/env perl
# clipprimer.pl - the workhorse of BAMClipper
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use IO::File;

##
# internal logic below
my $phredoffset = 33;
my $bedpe;
my $debug = 0;
my $window_upstream = 1;
my $window_downstream = 5;
GetOptions ("in=s" => \$bedpe,
	    "phredoffset=i" => \$phredoffset,
	    "upstream=i" => \$window_upstream,
	    "downstream=i" => \$window_downstream,
	    "debug" => \$debug,);
my %position2amplicon_positive;
my %position2amplicon_negative;
my $bedpe_fh = IO::File->new($bedpe) || die "ERROR: failed to open .bedpe panel description file";
BEDPE:while (<$bedpe_fh>) {
    chomp;
    next BEDPE if length($_) == 0 || substr($_,0,1) eq "#";
    my @fields = split ("\t", $_);
    die "ERROR: unrecognized format of the .bedpe file provided" if scalar @fields < 6;
    my $name = scalar @fields >= 7 ? $fields[6] : join("-",@fields[0..5]);
    # lookup up table from position to amplicon
    # 1-based pos of left and right primers
    my $primers = [$fields[0], $fields[1]+1, $fields[2], $fields[3], $fields[4]+1, $fields[5], $name];
    OFFSET:foreach my $offset ($window_upstream * -1 .. $window_downstream) {
        $position2amplicon_positive{$fields[0]}{$fields[1]+1+$offset} = $primers;
    }
    # 1-based pos of right primer
    OFFSET:foreach my $offset ($window_downstream * -1 .. $window_upstream) {
        $position2amplicon_negative{$fields[3]}{$fields[5]+$offset} = $primers;
    }
}
my $current_readname = "";
my @line_buffer;
LINE:while (<>) {
    my $original_line = $_;
    my $line = $original_line;
    chomp($line);
    if (length($line) == 0 || substr($line,0,1) eq "@") {
	# SAM header lines
	print $original_line;
	next LINE;
    }
    my @output;
    my @fields = split("\t", $line);
    if ($fields[0] ne $current_readname && $current_readname ne ""){
	process_line_buffer();
    }
    $current_readname = $fields[0];
    if ($fields[1] & 0x4 || $fields[1] & 0x0100) {
	# unmapped, or not primary alignment
	push @line_buffer, [\@fields,[0]];
	next LINE;
    }
    my ($cigar_total_length, @cigar) = &parse_cigar($fields[5]);
    my $rawrefseq = "N"x$cigar_total_length;
    my $readseq = $fields[9];
    my $readqual = $fields[10];
    my $direction = $fields[1] & 0x0010 ? "-" : "+";
    my @reconstruct_alignmentoutput = reconstruct_alignment($cigar_total_length, \@cigar, uc($readseq), $readqual, uc($rawrefseq), $fields[3], $fields[2], $direction);
    my @new_fields = @fields;
    if ($reconstruct_alignmentoutput[6] == 0) {
	# no M in CIGAR at all
	print "CIGAR $reconstruct_alignmentoutput[5] becomes *\n" if $debug;
	$new_fields[5] = "*";
	push @line_buffer, [\@new_fields,[1]];
    } elsif ($reconstruct_alignmentoutput[7] == 1) {
	# not trimmed, do nothing at the moment
	push @line_buffer, [\@fields,[0]];
    } else {
	# alignment updated
	$new_fields[3] = $reconstruct_alignmentoutput[3];
	$new_fields[5] = $reconstruct_alignmentoutput[5];
	if ($debug) {
	    print "to remove NM and MD from updated alignments: ".join(" ", grep {$_ =~ /^(?:NM|MD):/} @new_fields)."\n";
	}
	@new_fields = map {$new_fields[$_]} grep {$_ <= 10 || (substr($new_fields[$_],0,5) ne "NM:i:" && substr($new_fields[$_],0,5) ne "MD:Z:")} (0..$#new_fields);
	push @line_buffer, [\@new_fields,[2]];
    }
    if ($debug) {
	print "\@fields\n";
	print join("\t", @fields). "\n";
	print "\@new_fields\n";
	print join("\t", @new_fields). "\n";
	print "\n";
	
    }
}
if (scalar @line_buffer >= 1){
    process_line_buffer();
}

sub process_line_buffer {
    # flush buffer
    # [*]->[1]->[0] flag
    # 0: print original line
    # 1: unmapped after clipping
    # 2: clipped
    foreach my $i (0..$#line_buffer) {
	# unmapped after clipping
	if ($line_buffer[$i]->[1]->[0] == 1 && $line_buffer[$i]->[0]->[5] eq "*") {
	    # flag
	    if ($debug) {
		print "unmapped after clipping for read: ".$line_buffer[$i]->[0]->[0]."\n";
		print "original flag: ".$line_buffer[$i]->[0]->[1]."\n";
	    }
	    $line_buffer[$i]->[0]->[1] = $line_buffer[$i]->[0]->[1] | 4; # turn on unmapped bit 4 (0x4)
	    if ($debug) {
	        print "modified flag: ".$line_buffer[$i]->[0]->[1]."\n";
	    }
	}
    }
    # if >= 1 line is clipped
	# identify read paired, mapped (not unmapped) && primary alignment (not not primary alignment)
	#    save pos for first in pair as X
	#    save pos for second in pair as Y
	# assign PNEXT of every first in pair alignment as Y
	# assign PNEXT of every second in pair alignment as X
    my @id_of_clipped_lines = grep {$line_buffer[$_]->[1]->[0] == 2} (0..$#line_buffer);
    if (scalar @id_of_clipped_lines >= 1) {
	my @id_of_primary_alignments = grep {($line_buffer[$_]->[0]->[1] & 0x1) && !($line_buffer[$_]->[0]->[1] & 0x4) && !($line_buffer[$_]->[0]->[1] & 0x100)} (0..$#line_buffer);
	if (scalar @id_of_primary_alignments == 2){
	    my $pos_of_first_in_pair;
	    my $pos_of_second_in_pair;
	    if (($line_buffer[$id_of_primary_alignments[0]]->[0]->[1] & 0x40) && ($line_buffer[$id_of_primary_alignments[1]]->[0]->[1] & 0x80)) {
		# first in pair, then second in pair
		$pos_of_first_in_pair = $line_buffer[$id_of_primary_alignments[0]]->[0]->[3];
		$pos_of_second_in_pair = $line_buffer[$id_of_primary_alignments[1]]->[0]->[3];
	    } elsif (($line_buffer[$id_of_primary_alignments[1]]->[0]->[1] & 0x40) && ($line_buffer[$id_of_primary_alignments[0]]->[0]->[1] & 0x80)) {
		# second in pair, then first in pair
		$pos_of_first_in_pair = $line_buffer[$id_of_primary_alignments[1]]->[0]->[3];
		$pos_of_second_in_pair = $line_buffer[$id_of_primary_alignments[0]]->[0]->[3];
	    }
	    
	    my @id_of_first_in_pair_alignments = grep {$line_buffer[$_]->[0]->[1] & 0x40} (0..$#line_buffer);
	    my @id_of_second_in_pair_alignments = grep {$line_buffer[$_]->[0]->[1] & 0x80} (0..$#line_buffer);
	    
	    if (defined $pos_of_second_in_pair) {
		# assign PNEXT of every first in pair alignments as $pos_of_second_in_pair
		foreach my $i (@id_of_first_in_pair_alignments) {
		    print "original PNEXT is ".$line_buffer[$i]->[0]->[7]."\n" if $debug;
		    $line_buffer[$i]->[0]->[7] = $pos_of_second_in_pair;
		    print "modified PNEXT is ".$line_buffer[$i]->[0]->[7]."\n" if $debug;
		}
	    }
	    if (defined $pos_of_first_in_pair) {
		# assign PNEXT of every second in pair alignments as $pos_of_first_in_pair
		foreach my $i (@id_of_second_in_pair_alignments) {
		    print "original PNEXT is ".$line_buffer[$i]->[0]->[7]."\n" if $debug;
		    $line_buffer[$i]->[0]->[7] = $pos_of_first_in_pair;
		    print "modified PNEXT is ".$line_buffer[$i]->[0]->[7]."\n" if $debug;
		}
	    }
	}
    }
    print "\@line_buffer\n" if $debug;
    map { print join("\t",@{$_->[0]})."\n" } @line_buffer;
	
    @line_buffer = ();
}

sub reconstruct_alignment {
    my ($cigar_total_length, $cigar_arrayref, $readseq, $readqual, $refseq, $refpos_start, $chrom, $alignmentdirection) = @_;
    my @output;
    my @readseq = split ("", $readseq); # 1-based, first element is not used
    unshift @readseq, undef;
    
    my @cigar_op; # 1-based, first element is not used
    my $cigar_pos_offset = 1;
    
    my @cigar_refpos;
    
    my $readbase_pos_offset = 0;
    my $refpos_pos_offset = $refpos_start;
    foreach my $cigarop (@$cigar_arrayref) {
	my $len = $cigarop->[0];
	my $op = $cigarop->[1];
	map {$cigar_op[$_] = $op} ($cigar_pos_offset .. $cigar_pos_offset + $len - 1);
	map {$cigar_refpos[$_] = ($op eq "I" || $op eq "S") ? "*" : $refpos_pos_offset + $_ - $cigar_pos_offset} ($cigar_pos_offset .. $cigar_pos_offset + $len - 1);
	$readbase_pos_offset -= $len if $op eq "D";
	$refpos_pos_offset += $len if ($op ne "I" && $op ne "S");
	$cigar_pos_offset += $len;
    }
    if ($debug) {
	print join("\t", map {defined $_ ? $_ : "_"} @cigar_refpos[1..$cigar_total_length])."\n";
	print join("\t", map {defined $_ ? $_ : "_"} @cigar_op[1..$cigar_total_length])."\n";
    }
    
    my $print_original_line = 0;
    my ($left_alignpos, $right_alignpos) = min_max(grep {defined $_ && $_ ne "*"} @cigar_refpos[1..$cigar_total_length]);
    my $primers = undef;
    if ($alignmentdirection eq "+") {
	if (defined $position2amplicon_positive{$chrom}{$left_alignpos}) {
	    $primers = $position2amplicon_positive{$chrom}{$left_alignpos};
	}
    } elsif ($alignmentdirection eq "-") {
	if (defined $position2amplicon_negative{$chrom}{$right_alignpos}) {
	    $primers = $position2amplicon_negative{$chrom}{$right_alignpos};
	}
    }
    if (!defined $primers) {
	$print_original_line = 1;
    } else {
	my $trim_5prime_up_to;
	my $trim_3prime_up_to;
	# 5' trim up to end position of left primer
	POS: foreach my $i (1..$cigar_total_length) {
	    if ($cigar_refpos[$i] ne "*") {
		if ($cigar_refpos[$i] <= $primers->[2]) {
		    # keep last base of primer if there is deletion right after primer
		    $trim_5prime_up_to = $i unless ($i != $cigar_total_length && $cigar_op[$i] eq "M" && $cigar_op[$i+1] eq "D");
		} else {
		    last POS;
		}
	    }
	}
	if (defined $trim_5prime_up_to) {
	    print "trim_5prime_up_to: $trim_5prime_up_to\n" if $debug;
	    map {$cigar_op[$_] = $cigar_op[$_] eq "D" ? "*" : "S"; $cigar_refpos[$_] = "*"} (1..$trim_5prime_up_to);
	}
	# 3' trim up to start position of right primer
	POS: foreach my $i (reverse(1..$cigar_total_length)) {
	    if ($cigar_refpos[$i] ne "*") {
		if ($cigar_refpos[$i] >= $primers->[4]) {
		    # keep first base of primer if there is deletion right before primer
		    $trim_3prime_up_to = $i unless ($i != 1 && $cigar_op[$i-1] eq "D" && $cigar_op[$i] eq "M");
		} else {
		    last POS;
		}
	    }
	}
	if (defined $trim_3prime_up_to) {
	    print "trim_3prime_up_to: $trim_3prime_up_to\n" if $debug;
	    map {$cigar_op[$_] = $cigar_op[$_] eq "D" ? "*" : "S"; $cigar_refpos[$_] = "*"} ($trim_3prime_up_to..$cigar_total_length);
	}
	if ($debug) {
	    print join("\t", map {defined $_ ? $_ : "_"} @cigar_refpos[1..$cigar_total_length])."\n";
	    print join("\t", map {defined $_ ? $_ : "_"} @cigar_op[1..$cigar_total_length])."\n";
	}
    }
    my ($collapsed_cigar, $containM) = collapse_cigar(grep {$_ ne "*"} @cigar_op[1..$cigar_total_length]);
    my ($new_left_alignpos, $new_right_alignpos) = min_max(grep {defined $_ && $_ ne "*"} @cigar_refpos[1..$cigar_total_length]);
    return (defined $primers ? $primers->[6] : "undef",
	    $left_alignpos,
	    $right_alignpos,
	    $new_left_alignpos,
	    $new_right_alignpos,
	    $collapsed_cigar,
	    $containM,
	    $print_original_line);
}

sub min_max {
    my @array = sort {$a <=> $b} @_;
    return ($array[0], $array[$#array]);
}

sub parse_cigar {
    my ($cigar_string) = @_;
    my @cigar;
    my $i = 0;
    while ($cigar_string =~ m/([0-9]+)([MIDNSHPX=])/g) {
	# skip hard-clipping
	die "ERROR: Unexpected CIGAR operation $2" if $2 eq "N" || $2 eq "P";
	if ($2 ne "H") {
	    push @cigar, [$1, $2];
	    $i += $1;
	}
    }
    return ($i, @cigar);
}

sub collapse_cigar {
    my $output = "";
    my $current_op = "";
    my $current_len = 0;
    my $containM = 0;
    OP:foreach my $op (@_) {
	if ($op ne $current_op && $current_op ne "") {
	    $output .= sprintf("%d%s",$current_len,$current_op);
	    $current_len = 0;
	    $containM = 1 if $current_op eq "M" || $current_op eq "=" || $current_op eq "X";
	}
	$current_op = $op;
	$current_len++;
    }
    if ($current_op ne "") {
	$output .= sprintf("%d%s",$current_len,$current_op);
	$containM = 1 if $current_op eq "M" || $current_op eq "=" || $current_op eq "X";
    }
    return ($output, $containM);
}
