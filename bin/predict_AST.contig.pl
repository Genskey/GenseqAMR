#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use Getopt::Long;

my %opt=();
GetOptions(\%opt,"pathogen:s");
@ARGV || die "Usage: perl $0 <patho_ARG.detect> --pathogen <name> > predict.xls
Version: 3.0,20220103
Contact: p.han\@genskey.com\n\n";

### config
my $report = "$Bin/../DB/focus_patho_drug.report_WGScontig";
my $feature_weight = "$Bin/../DB/feature_weight.xls";
my $ARG_ppv = "$Bin/../DB/feature_ppv.stat";
my $card_info = "$Bin/../DB/GenseqResDB/card.info.xls";
my $gt_similarity = "$Bin/../DB/GenseqResDB/pairGene.similarity.xls";

### Each ARG's ppv obtained from training data set
open IN,$ARG_ppv || die $!;
my (%ARG_ppv,%ARG_ppv_freq,%nonkey_gt);
while(<IN>){
    chomp;
    my @ll = split /\t/;
    $ll[5] eq "0/0" && next;
    my $freq = (split /\//,$ll[5])[1];
    $ll[4] ||= 0.001;
    $ARG_ppv{$ll[0]}{$ll[1]}{$ll[3]} = $ll[4];
    $ARG_ppv_freq{$ll[0]}{$ll[1]}{$ll[3]} = (split /\//,$ll[5])[1];
    $freq >= 10 && $ll[4] < 0.6 && ($nonkey_gt{$ll[0]}{$ll[1]}{$ll[2]}{$ll[3]} = 1);
}
close IN;

### All pair ARG's similarity calculated based on the alignment
open SI,$gt_similarity;
my (%gt_sim);
while(<SI>){
    /^ARG1_gf/ && next;
    chomp;
    my @ll = split /\t/;
    $gt_sim{$ll[1]}{$ll[3]} ||= $ll[5];
    $gt_sim{$ll[1]}{$ll[3]} =  $gt_sim{$ll[1]}{$ll[3]} >= $ll[5] ? $gt_sim{$ll[1]}{$ll[3]} : $ll[5];
    $gt_sim{$ll[3]}{$ll[1]} = $gt_sim{$ll[1]}{$ll[3]};
}
close SI;

### Key ARGs with weight coefficients
open IN,$feature_weight || die $!;
my (%arg_weight,%argfamily_weight,%argfamily_ppv,%gf_keygt);
while(<IN>){
    chomp;
    /^#/ && next;
    my @ll = split /\t/;
    $arg_weight{$ll[0]}{$ll[1]}{$ll[7]} = $ll[8];
    $ll[3] =~ /variant/ && next;
    $argfamily_weight{$ll[0]}{$ll[1]}{$ll[4]} = $ll[5];
    $argfamily_ppv{$ll[0]}{$ll[1]}{$ll[4]} = (split /\(/,$ll[6])[0];
    $gf_keygt{$ll[0]}{$ll[1]}{$ll[4]}{$ll[7]} = 1;
}
close IN;

### card_info
open IN,$card_info || die $!;
my %arg2family;
while(<IN>){
    chomp;
    my @ll = split /\t/;
    $arg2family{$ll[9]} = $ll[8];
}
close IN;

############################### main process ###############################
my $pathogen_args = shift @ARGV;
open IN,$pathogen_args || die $!;
my (@args,%model_type);
while(<IN>){
    chomp;
    /^#/ && next;
    my @ll = split /\t/;
    if($ll[2] =~ /variant/){
	if($ll[3]){
	    my @aa = split /,/,$ll[3];
	    for my $v(@aa){
		my $f = "$ll[0]:$v";
		push @args,$f;
		$model_type{$f} = "variant";
	    }
	}
    }else{
	push @args,$ll[0];
    }
}
close IN;

### Predict AST for targeted pathogens
open IN,$report || die $!;
my %uniq_ARG_family;
print "#Patho\tdrug\tdetect_ARGs\tscores\tcutoff\tpredict\taccuracy\n";
while(<IN>){
    chomp;
    my @ll = split /\t/;
    my $index = 2;
    $ll[$index] eq "-" && next;
    $ll[0] eq $opt{pathogen} || next;
    my ($report_args,$scores) = ("",0);
    for my $arg(@args){
	if($model_type{$arg} && $model_type{$arg} eq "variant"){ #varient model
            if($arg_weight{$ll[0]}{$ll[1]}{$arg} && $arg_weight{$ll[0]}{$ll[1]}{$arg} > 0){
                $report_args .= "$arg;";
                $scores += $arg_weight{$ll[0]}{$ll[1]}{$arg};
            }
	}else{
	    my $arg_family = $arg2family{$arg};
	    $arg_weight{$ll[0]}{$ll[1]}{$arg} ||= 0;
	    $argfamily_weight{$ll[0]}{$ll[1]}{$arg_family} ||= 0; 
	    $argfamily_weight{$ll[0]}{$ll[1]}{$arg_family} <= 0 && next;
	    if($arg_weight{$ll[0]}{$ll[1]}{$arg} > 0){
	        $uniq_ARG_family{$ll[0]}{$ll[1]}{$arg_family}++;
	        $report_args .= "$arg;";
	        $scores += $arg_weight{$ll[0]}{$ll[1]}{$arg};
		my @tmp_args = split /;/,$report_args;
		my ($add,$report_args_tmp) = (0,"");
		for my $aa(@tmp_args){
                    $aa eq $arg_family && ($add = 1);
                    $aa ne $arg_family && ($report_args_tmp .= "$aa;");
                }
		$report_args = $report_args_tmp;
		$add && ($scores -= $argfamily_weight{$ll[0]}{$ll[1]}{$arg_family});
	    }else{
		$ARG_ppv_freq{$ll[0]}{$ll[1]}{$arg} && $ARG_ppv_freq{$ll[0]}{$ll[1]}{$arg} >= 10 && $ARG_ppv{$ll[0]}{$ll[1]}{$arg} < 0.6 && next;

		my $threhold = 0.95;
	        my @keygt = sort keys %{$gf_keygt{$ll[0]}{$ll[1]}{$arg_family}};
		my ($high_sim_key,$arg_weight) = (0,0);
		for my $gt(@keygt){
                    if($gt_sim{$gt}{$arg} && $gt_sim{$gt}{$arg} > $threhold){
                        $high_sim_key ||= $gt_sim{$gt}{$arg}; #
                        $arg_weight ||= "$gt\t$arg_weight{$ll[0]}{$ll[1]}{$gt}"; #
                        $high_sim_key < $gt_sim{$gt}{$arg} && ($arg_weight = "$gt\t$arg_weight{$ll[0]}{$ll[1]}{$gt}");
                    }
                }

		my @nonkey_gt = sort keys %{$nonkey_gt{$ll[0]}{$ll[1]}{$arg_family}};
                my $high_sim_nonkey;
                for my $gt(@nonkey_gt){
                    if($gt_sim{$gt}{$arg} && $gt_sim{$gt}{$arg} > $threhold){
                        $high_sim_nonkey ||= $gt_sim{$gt}{$arg}; #
                        $high_sim_nonkey < $gt_sim{$gt}{$arg} && ($high_sim_nonkey = $gt_sim{$gt}{$arg});
                    }
                }

		if($high_sim_key){ #gt与筛选出来的重要gt具有较高的序列相似性，往往功能上也具有相似性
                    my ($keygt,$keygt_weight) = (split /\t/,$arg_weight)[0,1];
                    my @tmp_args = split /;/,$report_args;
                    my $add = 0;
                    for my $aa(@tmp_args){$aa eq $keygt && ($add = 1);}
                    if(!$add){
                        $uniq_ARG_family{$ll[0]}{$ll[1]}{$arg_family}++;
                        if($uniq_ARG_family{$ll[0]}{$ll[1]}{$arg_family} == 1){
                            $report_args .= "$arg_family;";
                            $scores += $argfamily_weight{$ll[0]}{$ll[1]}{$arg_family};
                        }
                    }
                }else{
                    $high_sim_nonkey && next;
                    $uniq_ARG_family{$ll[0]}{$ll[1]}{$arg_family}++;
                    if($uniq_ARG_family{$ll[0]}{$ll[1]}{$arg_family} == 1){
                        $report_args .= "$arg_family;";
                        $scores += $argfamily_weight{$ll[0]}{$ll[1]}{$arg_family};
                    }
		}
	    }
        }
    }
    $scores = sprintf "%.6f", $scores;
    $report_args ||= "ND";
    $report_args =~ s/;$//;

    my $cutoff_R = $ll[2];
    my $cutoff_S = $ll[2];
    $cutoff_R = sprintf "%.6f", $cutoff_R;
    $cutoff_S = $cutoff_S ne "-" ? sprintf "%.6f", $cutoff_S : "-";
    my $drug_ppv = $ll[3];
    my $drug_npv = $ll[4];

    my $predict_ast = $scores >= $cutoff_R ? "R" : (($cutoff_S ne "-" && $scores < $cutoff_S) ? "S" : "/");
    my $accuracy = $predict_ast eq "R" ? $drug_ppv : ($predict_ast eq "S" ? $drug_npv : "/");
    my $cutoff = $predict_ast eq "R" ? $cutoff_R : ($predict_ast eq "S" ? $cutoff_S : "/");
    print "$ll[0]\t$ll[1]\t$report_args\t$scores\t$cutoff\t$predict_ast\t$accuracy\n";
}
close IN;

