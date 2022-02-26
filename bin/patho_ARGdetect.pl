#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use Getopt::Long;

my %opt=();
GetOptions(\%opt,"feature_weight:s","list_focus_species:s","patho_argCN:s","detect_patho:s","arg2taxon:s","pathogen:s");
@ARGV == 2 || die "Usage: perl $0 <arg.detect> <coverage.stat> [--option] >drug_detect.xls
Version: 1.0, 2022025
Option:
  <arg.detect>		  	set the file of detected ARGs
  <coverage.stat>		set the file of detected genome coverage and cover depth for targetd pathogen
  --feature_weight	[file]	set the feature weight file
  --patho_argCN		[file]  set the ARG copy number for target pathogen
  --list_focus_species	[file]  set the list of focus pathogens
  --detect_patho	[file]  set the profile file of detected pahtogens
  --arg2taxon		[file]  set the profile file of ARG's specied source obtained from read's annotation
  --pathogen		[str]	set the pathogen's name to output

Example:
  perl $0 <arg.detect> <coverage.stat> >drug_detect.xls
  perl $0 <arg.detect> <coverage.stat> --detect_patho <patho.detect> --arg2taxon <predict_arg2taxon> >drug_detect.xls\n\n";

### config
$opt{feature_weight} ||= "$Bin/../DB/feature_weight.xls";
$opt{patho_argCN} ||= "$Bin/../DB/patho_ARG.CopyNumber.stat";
$opt{list_focus_species} ||= "$Bin/../DB/focus_patho.list";

### Key ARG and ARG family with weight coefficients
open IN,$opt{feature_weight} || die $!;
my %key_gf;
while(<IN>){
    /^#/ && next;
    my @ll = split /\t/;
    $key_gf{$ll[4]} = 1;
}
close IN;

### the species source of ARGs obtained from read's annotation 
my %arg_taxon;
if($opt{arg2taxon}){
    open IN,$opt{arg2taxon} || die $!;
    while(<IN>){
        chomp;
        /^#/ && next;
        my @ll = split /\t/;
	$arg_taxon{$ll[0]} = "$ll[2]\t$ll[4]";
    }
    close IN;
}

### ARG copy number for target pathogen
my %spe_max_argCN;
if($opt{patho_argCN}){
    open IN,$opt{patho_argCN} || die $!;
    while(<IN>){
	chomp;
	my @ll = split /\t/;
	$spe_max_argCN{$ll[0]}{$ll[7]} = (split /,/,$ll[3])[1];
    }
    close IN;
}

### focus pathogen
open IN,$opt{list_focus_species} || die $!;
my (%focus_spe,%spe_rn,%spe_cov,%spe_covxdepth);
while(<IN>){
    chomp;
    s/^\s+//;s/\s+$//;
    $focus_spe{$_} = 1;
    $spe_rn{$_} ||= 1;  #
    $spe_covxdepth{$_} ||= 0.000001;
}
close IN;

################################ main process #################################
my ($arg_detectf,$patho_detectf) = @ARGV;
open IN,$patho_detectf || die $!;
my (%spe2genus_cov_rn,%total_bac_rn);
while(<IN>){
    chomp;
    my @ll = split /\t/;
    $spe_rn{$ll[0]} = $ll[2];
    $spe_covxdepth{$ll[0]} = ($ll[5]*$ll[4])/100;
    $spe2genus_cov_rn{$ll[0]} = "$ll[4]\t$ll[5]\t$ll[2]\tNA";
    $total_bac_rn{"Bacteria"} = "NA";
}
close IN;

if($opt{detect_patho} && -s $opt{detect_patho}){
    open IN,$opt{detect_patho} || die $!;
    my %uniq = ();
    $total_bac_rn{"Bacteria"} = 0;
    my ($bac_rn_index,$genusName_index,$genus_rn_index,$speName_index,$spe_rn_index,$cov_index,$depth_index);
    my $h = <IN>;chomp $h;
    my @hh = split /\t/,$h;
    for my $i(0..$#hh){
        $hh[$i] eq "Kingdom" && ($bac_rn_index = $i);
        $hh[$i] eq "Genus拉丁名" && ($genusName_index = $i);
        $hh[$i] eq "属原始reads数" && ($genus_rn_index = $i);
        $hh[$i] eq "种拉丁名" && ($speName_index = $i);
        $hh[$i] eq "种原始reads数" && ($spe_rn_index = $i);
        $hh[$i] eq "Coverage(ratio)" && ($cov_index = $i);
        $hh[$i] eq "Depth" && ($depth_index = $i);
    }
    while(<IN>){
        chomp;
        my @ll = split /\t/;
        $ll[$cov_index] eq "-" && next;
        my $cov = (split /\(/,$ll[$cov_index])[0];$cov =~ s/\%//;
        my $depth = $ll[$depth_index];
        if($focus_spe{$ll[$speName_index]}){
	    $spe_rn{$ll[$speName_index]} = $ll[$spe_rn_index]; 
	    $spe_covxdepth{$ll[$speName_index]} = ($cov*$depth)/100; 
	    $spe2genus_cov_rn{$ll[$speName_index]} = "$cov\t$depth\t$ll[$spe_rn_index]\t$ll[$genus_rn_index]"; 
        }
        $ll[$genusName_index] eq "-" && ($ll[$genusName_index] = $ll[$speName_index]); 
        $uniq{$ll[$genusName_index]}++;
        $uniq{$ll[$genusName_index]} < 2 && ($total_bac_rn{$ll[$bac_rn_index]}+=$ll[$genus_rn_index]); 
    }
    close IN;
}

### 
open IN,$arg_detectf || die $!;
my (%Patho_args,%is_locate,%not_locate);
my @focus_patho_sort = sort {$spe_rn{$b} <=> $spe_rn{$a}} keys %focus_spe;
while(<IN>){
    chomp;
    /^#/ && next;
    my @ll = split /\t/;
    my ($arg_taxons,$arg_born) = @ll[20,22];
    my @predict_spe;
    $arg_taxon{$ll[4]} && (@predict_spe = split /\t/,$arg_taxon{$ll[4]});
    my ($locate_mark,$cycle) = (0,1);
    CY:{;}
    for my $spe(@focus_patho_sort){
	$spe_covxdepth{$spe} == 0.000001 && next;
	my $arg_copy_n = sprintf "%.2f",($ll[16]*$ll[17])/$spe_covxdepth{$spe};
	$arg_copy_n = $arg_copy_n > 1 ? $arg_copy_n/$cycle : $arg_copy_n*$cycle;
	$spe_max_argCN{$spe}{$ll[6]} ||= 3;
	$cycle == 1 && (print STDERR "$ll[6]\t$arg_copy_n\t$spe\t$spe_max_argCN{$spe}{$ll[6]}+3\tCycle$cycle\n");
	$cycle > 1 && (print STDERR "$ll[6]\t$arg_copy_n\t$spe\t$spe_max_argCN{$spe}{$ll[6]}+3\tCycle$cycle\n");
	$cycle >= 50 && ($arg_copy_n = 1);
	($arg_copy_n >= $spe_max_argCN{$spe}{$ll[6]}+3 || $arg_copy_n <= 0.5) && next; 
	my $focus_genus = (split /\s+/,$spe)[0]; 
	if($arg_taxons =~ /$spe/ig || $arg_taxons =~ /$focus_genus/ig ){
	    $Patho_args{$spe} .= "$spe\tDB_taxon_anno\t$ll[4]\t$ll[5]\t$ll[6]\t$ll[7]|$ll[10]\t$arg_copy_n\t$ll[19]\t$ll[18]\n";
	    $is_locate{"$ll[4]\t$ll[6]"} = 1;
	    $locate_mark = 1;
	}
    }
    if($locate_mark == 0 && $arg_born =~ /plasmid/){
        for my $spe(@focus_patho_sort){
	    $spe_covxdepth{$spe} == 0.000001 && next;
            my $arg_copy_n = sprintf "%.2f",($ll[16]*$ll[17])/$spe_covxdepth{$spe};
            $arg_copy_n = $arg_copy_n > 1 ? $arg_copy_n/$cycle : $arg_copy_n*$cycle;
            $spe_max_argCN{$spe}{$ll[6]} ||= 3;
	    $cycle >= 50 && ($arg_copy_n = 1);
	    ($arg_copy_n >= $spe_max_argCN{$spe}{$ll[6]}+3 || $arg_copy_n <= 0.4) && next; 
	    $Patho_args{$spe} .= "$spe\tplasmid\t$ll[4]\t$ll[5]\t$ll[6]\t$ll[7]|$ll[10]\t$arg_copy_n\t$ll[19]\t$ll[18]\n";
            $is_locate{"$ll[4]\t$ll[6]"} = 1;
	    $locate_mark = 1;
        }
    }
    if($locate_mark == 0 && $arg_taxon{$ll[4]}){
        for my $spe(@focus_patho_sort){
            $spe_covxdepth{$spe} == 0.000001 && next;
            my $arg_copy_n = sprintf "%.2f",($ll[16]*$ll[17])/$spe_covxdepth{$spe};
            $arg_copy_n = $arg_copy_n > 1 ? $arg_copy_n/$cycle : $arg_copy_n*$cycle;
            $spe_max_argCN{$spe}{$ll[6]} ||= 3;
            $cycle >= 50 && ($arg_copy_n = 1);
            ($arg_copy_n >= $spe_max_argCN{$spe}{$ll[6]}+3 || $arg_copy_n <= 0.5) && next; 
            my $focus_genus = (split /\s+/,$spe)[0]; 
	    if($predict_spe[1] =~ /$spe/ig || $predict_spe[0] =~ /$focus_genus/ig){
		$Patho_args{$spe} .= "$spe\tread_taxon_anno\t$ll[4]\t$ll[5]\t$ll[6]\t$ll[7]|$ll[10]\t$arg_copy_n\t$ll[19]\t$ll[18]\n";
	        $is_locate{"$ll[4]\t$ll[6]"} = 1;
	        $locate_mark = 1;
	    }
        }
    }

    my $mark = $key_gf{$ll[4]} ? "Key" : "Non_key";
    if(!$is_locate{"$ll[4]\t$ll[6]"} && $cycle < 50){
	print STDERR "#$ll[4]\t$ll[6]\t$mark\tnot_locate\n";
	if($key_gf{$ll[4]}){$cycle++;goto (CY);}
    }else{
	print STDERR "#$ll[4]\t$ll[6]\t$mark\tlocate\n";
    }
    print STDERR "\n";
}
close IN;

###output
print "#focus_patho\tmethod\tdetect_gf\tgf_rn\tdetect_gt\tgt_rn(specific|ARG-like)\tcopy_num\tmodel_type\tVaration\tgenome_coverage(%)\tAvgDepth\tspe_rn\tgenus_rn\ttotal_Bac_rn\n";
for my $spe(sort keys %Patho_args){
    $opt{pathogen} && $opt{pathogen} ne $spe && next;
    my @lines = split /\n/,$Patho_args{$spe};
    for my $line(@lines){
	print "$line\t$spe2genus_cov_rn{$spe}\t$total_bac_rn{Bacteria}\n";
    }
}

