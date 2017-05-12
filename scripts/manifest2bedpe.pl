#!/usr/bin/env perl
# Convert primer design file from illumina manifest format to BEDPE format
# Usage: manifest2bedpe.pl < TruSight-Myeloid-Manifest.txt > trusight_myeloid.bedpe
use strict;
use warnings;
use Tie::IxHash;

my $debug = 0;
my $current_block;
my @column_headers;
my %required_column_headers = ("Probes" => ["Target ID", "ULSO Sequence", "DLSO Sequence"],
			       "Targets" => ["Target Number", "TargetA", "TargetB", "Start Position", "End Position", "Probe Strand"],
			       );
my (%probes, %targets);
my $ixhash_probes = tie(%probes, 'Tie::IxHash');
my $ixhash_targets= tie(%targets, 'Tie::IxHash');

# parsing
LINE:while (<>) {
    print "# $_" if $debug;
    chomp;
    s/\r$//;
    my @fields = split(/\t/);

    next LINE if scalar @fields == 0 || length($fields[0]) == 0;
    
    # define block
    if (substr($fields[0],0,1) eq "[") {
	if ($fields[0] eq "[Probes]") {
	    $current_block = "Probes";
	    @column_headers = ();
	} elsif ($fields[0] eq "[Targets]") {
	    $current_block = "Targets";
	    @column_headers = ();
	} else {
	    $current_block = undef;
	    @column_headers = ();
	}
	next LINE;
    }

    # define column header (only if block is defined)
    if (defined $current_block && scalar @column_headers == 0 && grep {$_ eq "Chromosome"} @fields) {
	foreach my $required_column (@{$required_column_headers{$current_block}}) {
		die "ERROR: Unexpected $current_block line that lacks $required_column: ".join("\t", @fields) unless grep {$_ eq $required_column} @fields;
	}
	print "# defined column for $current_block: ". join(",", @fields). "\n"  if $debug;
	@column_headers = @fields;
	next LINE;
    }
    
    #print join("\t", $current_block, @fields)."\n";
    
    if (defined $current_block && scalar @column_headers >= 1) {
	if ($current_block eq "Probes") {
	    my %probe = (map {$column_headers[$_] => $fields[$_]} 0..$#column_headers);
	    die "ERROR: Duplicated Target ID found: ".$probe{"Target ID"} if exists $probes{$probe{"Target ID"}};
	    $probes{$probe{"Target ID"}} = \%probe;
	} elsif ($current_block eq "Targets") {
	    my %target = (map {$column_headers[$_] => $fields[$_]} 0..$#column_headers);
	    next LINE if $target{"Target Number"} != 1 || $target{"TargetA"} ne $target{"TargetB"};
	    die "ERROR: This Target ID is defined in Targets but not in Probes: ".$target{"TargetA"} if !exists $probes{$target{"TargetA"}};
	    $targets{$target{"TargetA"}} = \%target;
	}
    }
}

# parsing done, convert to BEDPE
foreach my $target_key (keys %targets) {
    my $target = $targets{$target_key};
    print "# processing target: $target_key\n" if $debug;
    print join("\t"
	       , $target->{"Chromosome"}
	       , $target->{"Start Position"} - 1
	       , $target->{"Start Position"} + length($target->{"Probe Strand"} eq "+" ? $probes{$target_key}->{"DLSO Sequence"} : $probes{$target_key}->{"ULSO Sequence"}) - 1
	       , $target->{"Chromosome"}
	       , $target->{"End Position"} - length($target->{"Probe Strand"} eq "+" ? $probes{$target_key}->{"ULSO Sequence"} : $probes{$target_key}->{"DLSO Sequence"})
	       , $target->{"End Position"}
	       , $target_key
	       , $target->{"End Position"} - $target->{"Start Position"} + 1
	       , "+"
	       , "-"
	       )."\n";
}