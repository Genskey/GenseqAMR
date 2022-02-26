#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd qw(abs_path);

###################### get options ################################
my ($prefix,$outdir);
GetOptions("prefix:s"=>\$prefix,"outdir:s"=>\$outdir);
$outdir ||= ".";
-d $outdir || `mkdir -p $outdir`;
$outdir = abs_path($outdir);

###################### get path ####################################
my $glen_file = "$Bin/../DB/target_pathogens_Refgenome/all.sort.info";
my $speciseid = "$Bin/../DB/target_pathogens_Refgenome/DB.tax.format";

###################### help info ####################################
@ARGV == 1 || die"Usage:perl cover_depth.pl <sam> [option]	
    --outdir <str>  output dir
    --prefix <str>  prefix of output file
Example:
    perl cover_depth.pl sample.sam --outdir result  --prefix sample\n\n"; 

########################## prepare genome info #################################
my (%contig_len,%genome_len,%gbk2assID,%assID);
open IN,$glen_file;
while(<IN>){
    chomp;
    my @tmp=/\t/ ? split /\t/ : split /\s+/;
    $contig_len{$tmp[0]}= @tmp==7 ? $tmp[-3] : $tmp[-2]; #gbkID contigLen
    $genome_len{$tmp[2]} += @tmp==7 ? $tmp[-3] :$tmp[-2]; #genome_assID genomeLen
    $gbk2assID{$tmp[0]} = $tmp[2]; #gbkID genome_assID
    push @{$assID{$tmp[2]}},$tmp[0]; #genome_assID contigID
}
close IN;

######### get speciese info #######################################################
my %spe;
open IN,$speciseid;
while (<IN>){
    /^#/ && next;
    chomp;
    my @tmp=/\t/ ? split /\t/ : split /\s+/;
    $spe{$tmp[0]}=$tmp[5]; #assID speciese_name
}
close IN;

####################### read bam ###############################
my $sam=shift @ARGV;
open IN,$sam || die $!;
my (%uniq,%cvg);
while(<IN>){
    /^\@/ && next;
    chomp;
    my @arr=/\t/ ? split /\t/ : split /\s+/;
    $arr[2] eq "*" && next;
    $uniq{$gbk2assID{$arr[2]}}{$arr[0]}++;
    $uniq{$gbk2assID{$arr[2]}}{$arr[0]} > 1 && next;
    my $maplen = 0;
    for my $n($arr[5] =~ /(\d+)[MID]/g){
	$maplen += $n;
    }
    my $start = $arr[3];
    my $end = $arr[3] + $maplen - 1;
    for my $pos($start..$end){
	$cvg{$gbk2assID{$arr[2]}}{$arr[2]}{$pos}++;
    }
}
close IN;

open OUT,">","$prefix.cvgstat" or die $!;
foreach my $ass_id (sort keys %cvg){
    my $reads_num = keys %{$uniq{$ass_id}};
    my ($cover_tlen,$cover_bases_tn);
    for my $cid(sort keys %{$cvg{$ass_id}}){
	my @POSs = sort {$a<=>$b} keys %{$cvg{$ass_id}{$cid}};
	my $cover_clen = @POSs;
	$cover_tlen += $cover_clen;
	for my $p(@POSs){
	    $cover_bases_tn += $cvg{$ass_id}{$cid}{$p};
	}
    }
    my $cvstr = "$cover_tlen/$genome_len{$ass_id}";
    my $cvrate = sprintf "%.3f",($cover_tlen/$genome_len{$ass_id})*100;
    my $dep = sprintf "%.1f",($cover_bases_tn/$cover_tlen);
    print OUT "$spe{$ass_id}\t$ass_id\t$reads_num\t$cvstr\t$cvrate\t$dep\n";
}
close OUT;

