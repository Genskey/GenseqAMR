#!/usr/bin/perl -w
use FindBin qw($Bin);

@ARGV || die "Usage: perl $0 <ARG.stat> > ARG.stat.filter\n";
my $gt_similarity = "$Bin/../DB/GenseqResDB/pairGene.similarity.xls";
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

###
open IN,$ARGV[0] || die $!;
my (%uniq,%uniq2,%group_top1,%gt_cover);
while(<IN>){
    chomp;
    if(/^#/){print $_,"\n";next;}
    my @ll = split /\t/;
    $uniq{$ll[0]}{$ll[2]}++;
    $uniq2{$ll[0]}{$ll[2]}{$ll[4]}++;
    $group_top1{$ll[2]} ||= $ll[6];
    my $mu = $ll[10]/$ll[1]; #
    $mu > 2 && next; #
    my @cov = split /\//,$ll[15];
    ($cov[0] == 0 && $ll[16]<=0.6) && next;
    if($ll[10] >= 3 && $ll[16] <= 0.3){ #
	($ll[14] >= 2 || $ll[17] <3) || next;
    }
    $ll[16] <= 0.6 && $ll[17] > 5 && next;
    $gt_cover{$ll[6]}= $ll[16]; 
    my @detect = keys %retain;
    if($uniq{$ll[0]}{$ll[2]} < 2){
	$ll[3] < 1 && next; #
    }else{
	if($retain{$group_top1{$ll[2]}}){  #
	    my $gt_sim_max=0;
	    my $gt_sim_max_arg="";
	    for my $d(@detect){
		if($gt_sim{$ll[6]}{$d} && $gt_sim_max <= $gt_sim{$ll[6]}{$d}){
		    $gt_sim_max = $gt_sim{$ll[6]}{$d};
		    $gt_sim_max_arg = $d;
		}
	    }
	    if($gt_sim_max > 0.95){ #
		if($cov[0] == $cov[2]){ #
		    ;
		}else{
		    if($cov[1] == $cov[2]){
			my $spe_rn_cutoff = sprintf "%.0f",$ll[10]*((100-$gt_sim_max)/100);
			$ll[7] > $spe_rn_cutoff || next;
		    }else{
			next;
		    }
		}
	    }else{
		if($cov[0] == $cov[2] || $cov[1] == $cov[2]){
		    ;
		}else{
		    if($uniq2{$ll[0]}{$ll[2]}{$ll[4]} < 2){
		        $ll[5] < 1 && next;
		    }else{
		        if($gt_sim_max_arg){
		            $ll[16] >= $gt_cover{$gt_sim_max_arg} || next; #
		        }else{ 
		            $ll[9] > 2 || next; #
		        }
		    }
		}
	    }
	}else{
	    $cov[0] == $cov[2] || next; 
	}
    }
    $retain{$ll[6]} = 1; #
    print join("\t",@ll),"\n";
}
close IN;

