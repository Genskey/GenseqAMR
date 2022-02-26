#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use Getopt::Long;

my ($WGSread,$mNGSread,$min_rn,$ppv_cutoff,$out_gf_score);
GetOptions(
  "WGSread" => \$WGSread,
  "mNGSread" => \$mNGSread,
  "min_rn:i" => \$min_rn,
  "ppv_cutoff:i" => \$ppv_cutoff,
  "out_gf_score" => \$out_gf_score,
);
$ppv_cutoff ||= 0.6;
$min_rn ||= 2000;

@ARGV || die "Usage: perl $0 <patho_ARG.detect> [--option] > predict.xls
Version: 3.0,20220103
Contact: p.han\@genskey.com
Option:
  --WGSread		  if set, select report templete file for WGS read
  --mNGSread		  if set, select report templete file for mNGS read
  --min_rn	[int]	  set the minimum number of detected reads for target pathogen to predict AST result, defalt is 2000
  --ppv_cutoff	[float]   set ppv cutoff for non-key ARG, defalt is 0.6
  --out_gf_score	  if set, output the score calculated based on ARG family
\n\n";

### config
my $report = $mNGSread ? "$Bin/../DB/focus_patho_drug.report_mNGSread" : ($WGSread ? "$Bin/../DB/focus_patho_drug.report_WGSread" : "$Bin/../DB/focus_patho_drug.report_WGScontig");
my $feature_weight = "$Bin/../DB/feature_weight.xls";
my $ARG_ppv = "$Bin/../DB/feature_ppv.stat";
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
    $freq >= 10 && $ll[4] < $ppv_cutoff && ($nonkey_gt{$ll[0]}{$ll[1]}{$ll[2]}{$ll[3]} = 1);
}
close IN;

### All pair ARG's similarity calculated based on the alignment
open SI,$gt_similarity || die $!;
my %gt_sim;
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

############################### main process ###############################
my $pathogen_args = shift @ARGV;
open IN,$pathogen_args || die $!;
my (%patho_ARGdetect,%model_type,%patho_cov,%patho_rn,%arg2family);
while(<IN>){
    chomp;
    /^#/ && next;
    my @ll = split /\t/;
    $patho_cov{$ll[0]} = $ll[9];  #the coverage of targeted pathogen is as a inicator of Sequencing data size
    $patho_rn{$ll[0]} = $ll[11];
    $arg2family{$ll[4]} = $ll[2];
    if($ll[7]=~/variant/){
	if($ll[8]){
	    my @aa = split /,/,$ll[8];
	    for my $v(@aa){
                my $f = "$ll[4]:$v";
                push @{$patho_ARGdetect{$ll[0]}},$f;
                $model_type{$f} = "variant";
            }
	}
    }else{
	push @{$patho_ARGdetect{$ll[0]}},$ll[4];
    }
}
close IN;

### Predict AST for targeted pathogens
open IN,$report || die $!;
my (%uniq_ARG_family,%uniq_gf);
my $header = $out_gf_score ? "#Patho\tdrug\tdetect_ARGs\tscores\tcutoff\tpredict\taccuracy\tpatho_rn\tpatho_coverage\tgf_scores\n" : "#Patho\tdrug\tdetect_ARGs\tscores\tcutoff\tpredict\taccuracy\tpatho_rn\tpatho_coverage\n";
print $header;
while(<IN>){
    chomp;
    my @ll = split /\t/;
    my $index = $mNGSread ? 10 : ($WGSread ? 6 : 2);
    $ll[$index] eq "-" && next;  #
    if($patho_cov{$ll[0]}){
	my ($report_args,$scores,$scores_gf) = ("",0,0);
	for my $arg(@{$patho_ARGdetect{$ll[0]}}){
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
		$uniq_gf{$ll[0]}{$ll[1]}{$arg_family}++;
		if($uniq_gf{$ll[0]}{$ll[1]}{$arg_family}==1){$scores_gf += $argfamily_weight{$ll[0]}{$ll[1]}{$arg_family};}
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
		    if($patho_cov{$ll[0]} >= 95){ #The amount of equencing data is sufficient
			$ARG_ppv_freq{$ll[0]}{$ll[1]}{$arg} && $ARG_ppv_freq{$ll[0]}{$ll[1]}{$arg} >= 10 && $ARG_ppv{$ll[0]}{$ll[1]}{$arg} < $ppv_cutoff && next;
		    }
		    my $threhold = 0.95; #similarity threhold
		    my @keygt = sort keys %{$gf_keygt{$ll[0]}{$ll[1]}{$arg_family}};
                    my ($high_sim_key,$arg_weight) = (0,0);
                    for my $gt(@keygt){
                        if($gt_sim{$gt}{$arg} && $gt_sim{$gt}{$arg} > $threhold){
                            $high_sim_key ||= $gt_sim{$gt}{$arg}; #
                            $arg_weight ||= "$gt\t$arg_weight{$ll[0]}{$ll[1]}{$gt}"; #
                            $high_sim_key < $gt_sim{$gt}{$arg} && ($arg_weight = "$gt\t$arg_weight{$ll[0]}{$ll[1]}{$gt}");
                        }
                    }

		    ##
                    my @nonkey_gt = sort keys %{$nonkey_gt{$ll[0]}{$ll[1]}{$arg_family}};
                    my $high_sim_nonkey;
                    for my $gt(@nonkey_gt){
                        if($gt_sim{$gt}{$arg} && $gt_sim{$gt}{$arg} > $threhold){
                            $high_sim_nonkey ||= $gt_sim{$gt}{$arg}; 
                            $high_sim_nonkey < $gt_sim{$gt}{$arg} && ($high_sim_nonkey = $gt_sim{$gt}{$arg});
                        }
                    }

                    if($high_sim_key){ #
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
                        $high_sim_nonkey && next; #
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
        $scores_gf = sprintf "%.6f", $scores_gf;
	$report_args ||= "ND";
        $report_args =~ s/;$//;
	$patho_cov{$ll[0]} eq "-" && ($patho_cov{$ll[0]} = 1);

	my $cutoff_R = $mNGSread ? $ll[10] : ($WGSread ? $ll[6] : $ll[2]);
	$cutoff_R = sprintf "%.6f", $cutoff_R;
	my $drug_ppv = $mNGSread ? ($ll[11] ne "-" ? sprintf "%.3f",$ll[11] : "-") : ($WGSread ? $ll[7] : $ll[3]);
	my ($cutoff_S,$drug_npv);
	if($mNGSread){
	    if($ll[12] ne "-" && $patho_cov{$ll[0]} >= $ll[12]){
		$cutoff_S = $ll[13];
		$drug_npv = $ll[14]
	    }elsif($ll[15] ne "-" && $patho_cov{$ll[0]} >= $ll[15]){
		$cutoff_S = $ll[16];
		$drug_npv = $ll[17];
	    }else{
		$cutoff_S = "-";
	    }
	}elsif($WGSread){
	    $cutoff_S = $ll[6];
	    $drug_npv = $ll[8];
	}else{
	    $cutoff_S = $ll[2];
	    $drug_npv = $ll[4];
	}
	$cutoff_S = $cutoff_S ne "-" ? sprintf "%.6f", $cutoff_S : "-";
	my $predict_ast = $scores >= $cutoff_R ? "R" : (($cutoff_S ne "-" && $scores < $cutoff_S) ? "S" : "/");
	my $accuracy = $predict_ast eq "R" ? $drug_ppv : ($predict_ast eq "S" ? $drug_npv : "/");

	my $out = $out_gf_score ? "$ll[0]\t$ll[1]\t$report_args\t$scores\t$cutoff_R\t$predict_ast\t$accuracy\t$patho_rn{$ll[0]}\t$patho_cov{$ll[0]}\t$scores_gf\n" : "$ll[0]\t$ll[1]\t$report_args\t$scores\t$cutoff_R\t$predict_ast\t$accuracy\t$patho_rn{$ll[0]}\t$patho_cov{$ll[0]}\n";
	$patho_rn{$ll[0]} >= $min_rn && (print $out);
    }
}
close IN;

