#!/usr/bin/perl -w 

use strict;

my $feature = "piRNA";
my $strand = 0; ## Strand of read with respect to annotation as defined by featureCounts. Values can be 0, 1, 2.
my $pdir = "/home/tvavouri/Projects/mouse-piRNA-variation/out";

my $sampleCountsTableF = "$pdir/featureCounts/$feature/norepeats/$feature.norepeats.featureCounts.s$strand.tsv";

my %ft_sample_count = (); ## Hash of hashes with feature (eg piRNA cluster id) as key and value a hash  where sample is the key and count the value
my @sampleIDs = (); ## An array of all sample names found here
my @fts = (); ## An array of features, i.e. piRNA clusters



## Exit if output counts file exists.
if (-e $sampleCountsTableF){
    print "File $sampleCountsTableF already exists\n";
    exit;
}

my $ftdir = $pdir."/featureCounts/".$feature."/norepeats";

print $ftdir . "\n";

foreach my $countsf (glob("$ftdir/sample??.s".$strand.".counts")) {
   	print $countsf . "\n";
    $countsf =~ /(sample\d\d)/;
    my $sampleID = $1;
    push(@sampleIDs,$sampleID);
    @fts = (); ## An array of features, e.g. piRNA clusters
    open my $fh, "<", $countsf or die "can't read open $countsf\n";
    while (<$fh>) {
	my $line = $_;
	next if $line =~ /Program/;
	next if $line =~ /Geneid/;	      
	chomp($line);
	my @info = split(/\t/,$line);
	$ft_sample_count{$info[0]}{$sampleID} = $info[-1]; 
	push(@fts,$info[0]);
    }
    close $fh or die "can't read close $countsf\n";
}

print $sampleCountsTableF . "\n";
my $fh;
open($fh,">", $sampleCountsTableF) or die "Could not open file $sampleCountsTableF for writing: $!";
foreach my $sampleID(@sampleIDs){
    print $fh "\t$sampleID";
}
print $fh "\n";
foreach my $ft (@fts){
    print $fh $ft;
    foreach my $sampleID(@sampleIDs){
	print $fh "\t$ft_sample_count{$ft}{$sampleID}";
    }
    print $fh "\n";
}
close $fh or warn "Close failed: $!";
