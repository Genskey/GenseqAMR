#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opt = ();
GetOptions(\%opt,"sub_coverage:s");
$opt{sub_coverage} ||= 0.6;

@ARGV || die "Usage: perl $0 <m6.filter.anno> > arg.anno\n";
open IN,$ARGV[0];
my (%hash,%contig_region,%arg_name,%model_type,%arg_len,%cover_len,%snp,%insert,%delete);
while(<IN>){
    chomp;
    my @ll=split /\t/;
    my $hit_mark = ($ll[4] == $ll[6] && $ll[5] == 100) ? "exact" : "no_exact";
    my $arg=(split /;/,$ll[19])[-2];
    $arg_name{$ll[3]} = $arg;
    $model_type{$ll[3]} = $ll[20];
    $arg_len{$ll[3]} = $ll[4];
    my $hit_region = "$ll[0]\t$ll[1]\t$ll[2]";
    ##
    $cover_len{$hit_region}{$ll[3]} += $ll[6];
    $hash{$hit_region}{$ll[3]} .= "$hit_mark--$arg;";
    my $strand = $ll[12] < $ll[13] ? "+" : "-";
    $contig_region{$hit_region}{$ll[3]} .= "$ll[10]-$ll[11]|$strand,";
    ##var
    $ll[16] ne "0" && ($snp{$hit_region}{$ll[3]} .= "$ll[16],");
    $ll[17] ne "0" && ($insert{$hit_region}{$ll[3]} .= "$ll[17],");
    $ll[18] ne "0" && ($delete{$hit_region}{$ll[3]} .= "$ll[18],");
}

for my $k(sort keys %hash){
    my @ref_ids = sort keys %{$hash{$k}};
    $hash{$k}{$ref_ids[0]} =~ s/;$//;
    my $cover = $cover_len{$k}{$ref_ids[0]}/$arg_len{$ref_ids[0]};
    $cover >= $opt{sub_coverage} || next;
    my $mark = $hash{$k}{$ref_ids[0]} =~ /;|no_exact/ ? "no_exact" : "exact";
    $contig_region{$k}{$ref_ids[0]} =~ s/,$//;
    $snp{$k}{$ref_ids[0]} ||= 0;
    $snp{$k}{$ref_ids[0]} =~ s/,$//;
    $insert{$k}{$ref_ids[0]} ||= 0;
    $insert{$k}{$ref_ids[0]} =~ s/,$//;
    $delete{$k}{$ref_ids[0]} ||= 0;
    $delete{$k}{$ref_ids[0]} =~ s/,$//;
    if(@ref_ids == 1){
	print "$k:$contig_region{$k}{$ref_ids[0]}\t$mark\t$arg_name{$ref_ids[0]}\t$model_type{$ref_ids[0]}\t$snp{$k}{$ref_ids[0]}\t$insert{$k}{$ref_ids[0]}\t$delete{$k}{$ref_ids[0]}\t-\n";
    }else{
	my $other_arg;
	for my $i(1..$#ref_ids){
	    $other_arg .= "$arg_name{$ref_ids[$i]};";
	}
	$other_arg =~ s/;$//;
	print "$k:$contig_region{$k}{$ref_ids[0]}\t$mark\t$arg_name{$ref_ids[0]}\t$model_type{$ref_ids[0]}\t$snp{$k}{$ref_ids[0]}\t$insert{$k}{$ref_ids[0]}\t$delete{$k}{$ref_ids[0]}\t$other_arg\n";
    }
}

