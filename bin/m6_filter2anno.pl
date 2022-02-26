#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use Getopt::Long;
use List::Util qw(min max);
use POSIX;

my %opt = ();
GetOptions(\%opt,"add_arg_anno:s","long_query_seq","identity:s","que_coverage:s","sub_coverage:s","mis:s","gap:s","evalue:s");
@ARGV==2 || die "Usage:perl $0 <origin.m6> <filter.m6> [--options]
Description: script to filter m6 file and add annotation infomation
Parameters:
    <origin.m6>			input m6 file
    <filter.m6>			output filtered m6 file with or without annotation
    --long_query_seq		if set ,it means query seq is long, such as ONT reads or assembled genome contig(usuallly length > 300bp)
    --identity	<float>		set identity threshold, by which filter the hits
    --evalue       <float>      set e-value threshold, by which filter the hits, default is 1e-5
    --mis	<int>		set max mismatch to filter hits
    --gap	<int>		set max gap to filter hits 
    --que_coverage <float>	set query coverage threshold, by which filter the hits
    --sub_coverage <float>	set sub coverage threshold, by which filter the hits
    --add_arg_anno <file>	set to add arg annotation info

Example: 
    perl $0 <origin.m6> <filter.m6>\n\n";

#======== set default parameters =============
$opt{add_arg_anno} ||= "$Bin/../DB/GenseqResDB/card.info.xls";
$opt{evalue} ||= 1e-5;
$opt{identity} ||= $opt{long_query_seq} ? 90 : 90;
$opt{que_coverage} ||= 0.5; #
$opt{sub_coverage} ||= 0.4; #

#======== ARG annotation infomation ==========
my (%ST2lingage,%ID_anno,%ID_ST);
open IN1,$opt{add_arg_anno} || die $!;  #L6 rank info
while(<IN1>){
    chomp;
    next if(/^seq_nucl_id/);
    my @ll = split /\t/;
    $ID_ST{$ll[0]} = $ll[10];
    $ID_anno{$ll[0]} = join("\t",@ll[2,3,11..$#ll]);
    my $lingage = join ";",@ll[5..10];
    $ST2lingage{$ll[10]} = $lingage;
}
close IN1;

#=========== main function ======================
my($input,$output)=@ARGV;
my (%uniq,%uniq2,%retain_firsthit,%filter_hit,%region,%ref_region,$num,$num_tmp,%score,%line);
open IN,$input || die $!;
while(<IN>){     #filter
    chomp;
    my @ll = /\t/ ? (split /\t/) : (split /\s+/);
    my ($que_len,$sub_len,$align_iden,$align_len,$nident,$mis,$gap,$evalue,$score) = @ll[1,3,4,5,6,7,8,13,14];
    my $que_cover = $align_len >= $que_len ? 1 : sprintf("%.3f",$align_len/$que_len);
    my $sub_cover = $align_len >= $sub_len ? 1 : sprintf("%.3f",$align_len/$sub_len);
    my $que_end_yes = ($ll[9] == 1 || $ll[9] == $ll[1] || $ll[10] == $ll[1] || $ll[10] == 1) ? 1 : 0;
    my $sub_end_yes = ($ll[11] == 1 || $ll[11] == $ll[3] || $ll[12] == 1 || $ll[12] == $ll[3]) ? 1 : 0;
    $evalue > $opt{evalue} && next;
    $align_iden < $opt{identity} && next;
    $opt{mis} && $mis > $opt{mis} && next;
    $opt{gap} && $gap > $opt{gap} && next;
    $uniq2{"$ll[0]\t$ll[2]"}++;
    if($opt{long_query_seq}){
	if($uniq2{"$ll[0]\t$ll[2]"} == 1){
            if($sub_end_yes){
	        my $align_len_cutoff = $que_end_yes ? 0.2*$sub_len : $opt{sub_coverage}*$sub_len;
	        $align_len < $align_len_cutoff && next;
		$retain_firsthit{"$ll[0]\t$ll[2]"} = 1;
            }else{
	        $sub_cover < $opt{sub_coverage} && next; #
		$retain_firsthit{"$ll[0]\t$ll[2]"} = 1;
            }
	}else{
	    $retain_firsthit{"$ll[0]\t$ll[2]"} || next;
	}
    }else{
        if($sub_end_yes){
	    ;
	}else{
	    ($opt{que_coverage} && $que_cover < $opt{que_coverage}) && next; #
	}
    }

    $uniq{$ll[0]}++;
    if($que_len <= 300){ #
	$uniq2{"$ll[0]\t$ll[2]"} > 1 && next; 
	$score{"$ll[0]_R1"}{$ll[2]} = $ll[14];
	$line{"$ll[0]_R1"}{$ll[2]} = join("\t",@ll[0,1],"R1",@ll[2..$#ll])."\n";
	next;
    }
    ###
    if($uniq{$ll[0]} == 1){
	$num = 1;
	##
	my ($start,$end) = $ll[9] < $ll[10] ? @ll[9,10] : @ll[10,9];
	for my $i (${start}..$end){ $region{$ll[0]}{$num}{$i} = 1; }
	##
	my ($ref_start,$ref_end) = $ll[11] < $ll[12] ? @ll[11,12] : @ll[12,11];
        for my $j(${ref_start}..$ref_end ){ $ref_region{$ll[0]}{$num}{$ll[2]}{$j} = 1; }
	$score{"$ll[0]_R$num"}{$ll[2]} = (($ll[6]/$ll[3])**3)*$ll[14];
	$line{"$ll[0]_R$num"}{$ll[2]} = join("\t",@ll[0,1],"R1",@ll[2..$#ll])."\n";
	next;
    }
    if($uniq{$ll[0]} >= 2){ #
        my $mark_new = 1;
	my @nums = sort {$b <=> $a} keys %{$region{$ll[0]}};
	my ($common_len,$com_percent);
	my ($start,$end) = $ll[9] < $ll[10] ? @ll[9,10] : @ll[10,9];
	for my $n(@nums){
	    ($common_len,$com_percent) = (0,0);
	    for my $j (${start}..$end){	$region{$ll[0]}{$n}{$j} && ($common_len++); }
	    $com_percent = sprintf "%.3f",$common_len/($end-$start+1);
	    if($com_percent > 0.05){ #
		$mark_new = 0;
		$num_tmp = $n;
		last;
	    }
	}
	my ($ref_start,$ref_end) = $ll[11] < $ll[12] ? @ll[11,12] : @ll[12,11];
	if($mark_new){   # 
	    my @gid = keys %{$ref_region{$ll[0]}{$num}}; #
	    my $exist = 0;
	    for (@gid){$_ eq $ll[2] && ($exist = 1);}
	    if($exist){
		my ($c_len,$c_percent) = (0,0);
		my ($ref_start,$ref_end) = $ll[11] < $ll[12] ? @ll[11,12] : @ll[12,11];
		for my $j (${ref_start}..$ref_end){$ref_region{$ll[0]}{$num}{$ll[2]}{$j} && ($c_len++);}
		$c_percent = sprintf "%.3f",$c_len/($ref_end-$ref_start);
		if($c_percent > 0.1){  #
		    $sub_cover > $opt{sub_coverage} || next; #
		    $num++;
		}
	    }else{
		$num++;
	    }
	    for my $j (${start}..$end){ $region{$ll[0]}{$num}{$j} = 1; }
            for my $j(${ref_start}..$ref_end ){ $ref_region{$ll[0]}{$num}{$ll[2]}{$j} = 1; }
            $score{"$ll[0]_R$num"}{$ll[2]} = $score{"$ll[0]_R$num"}{$ll[2]} ? ($score{"$ll[0]_R$num"}{$ll[2]} + (($ll[6]/$ll[3])**3)*$ll[14]) : (($ll[6]/$ll[3])**3)*$ll[14];
            $line{"$ll[0]_R$num"}{$ll[2]} = $line{"$ll[0]_R$num"}{$ll[2]} ? ($line{"$ll[0]_R$num"}{$ll[2]}.join("\t",@ll[0,1],"R$num",@ll[2..$#ll])."\n") : join("\t",@ll[0,1],"R$num",@ll[2..$#ll])."\n";
	}else{
	    for my $j (${start}..$end){ $region{$ll[0]}{$num_tmp}{$j} = 1; }
	    for my $j(${ref_start}..$ref_end ){ $ref_region{$ll[0]}{$num_tmp}{$ll[2]}{$j} = 1; }
            $score{"$ll[0]_R$num_tmp"}{$ll[2]} = $score{"$ll[0]_R$num_tmp"}{$ll[2]} ? ($score{"$ll[0]_R$num_tmp"}{$ll[2]} + (($ll[6]/$ll[3])**3)*$ll[14]) : (($ll[6]/$ll[3])**3)*$ll[14];
            $line{"$ll[0]_R$num_tmp"}{$ll[2]} = $line{"$ll[0]_R$num_tmp"}{$ll[2]} ? ($line{"$ll[0]_R$num_tmp"}{$ll[2]}.join("\t",@ll[0,1],"R$num_tmp",@ll[2..$#ll])."\n") : join("\t",@ll[0,1],"R$num_tmp",@ll[2..$#ll])."\n";
	}
    }
}
close IN;

### output
open OUT,">$output" || die $!;
foreach my $query_region (sort keys %score){
    my %score_hit;
    foreach my $sub_id (sort keys %{$score{$query_region}}){
	$score_hit{$score{$query_region}{$sub_id}} .= "$query_region,$sub_id\t";
    }
    my @scoces = sort {$b<=>$a} keys %score_hit;
    $score_hit{$scoces[0]} =~ s/\t$//;
    my @hits = split /\t/,$score_hit{$scoces[0]};
    for my $hit(@hits){
	my ($query_id,$sub_id) = split /,/,$hit;
	my @lines = split /\n/,$line{$query_id}{$sub_id};
	for (@lines){
	    print OUT "$_\t$ST2lingage{$ID_ST{$sub_id}}\t$ID_anno{$sub_id}\n";
	}
    }
}
close OUT;

