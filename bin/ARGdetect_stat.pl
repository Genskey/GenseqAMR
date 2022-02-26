#!usr/bin/perl -w
use strict;
use Getopt::Long;
my %opt = ();
@ARGV || die "Usage: perl $0 <read.anno.xls> <ARG.stat.xls> <reads.anno.lca> <reads.anno.lca.exact> 2> log
Description: a script to calculate reads number, coverage and depth for detected ARGs\n\n";

my ($anno,$statf,$lcaf,$lcaf_exact) = @ARGV;
### LCA annotation for read and stat
open IN,$anno || die $!;
my (%arg_len,%arg_match_status,%specific_match_var,%arg_cover,%arg_cover_exact,%arg_name,%L1_L2,%L1_L3,%hit_n_all,%hit_n_all_exact,%lingage,%anno_info,$reads_lca_anno,$reads_lca_anno2);
my ($lastid,$lastfull,$lastanno,$lastid2,$lastfull2,$lastanno2) = ('','','','','','','','');
while(<IN>){
    chomp;
    my @ll = split /\t/;
    my $sub_end_yes = 0;
    if($ll[12] < $ll[13]){
	(($ll[12] == 1 && $ll[11] == $ll[1]) || ($ll[13] == $ll[4] && $ll[10] == 1)) && ($sub_end_yes = 1);
    }else{
	(($ll[12] == $ll[4] && $ll[10] == 1) || ($ll[13] == 1 && $ll[11] == $ll[1])) && ($sub_end_yes = 1);
    }
    my %specific_match = ();
    if($ll[16] ne 0){
	my @tmp_X = split /,/,$ll[16];
	for my $i(0..$#tmp_X){ my ($p,$snp) = ($1,$2) if($tmp_X[$i]=~/([\d]+)\((\S+)\)/); $specific_match{$ll[3]}{$p} = "X";$specific_match_var{$ll[3]}{$p} = "$snp"; }
    }
    if($ll[17] ne 0){
	my @tmp_I = split /,/,$ll[17];
	for my $i(0..$#tmp_I){ my ($p,$bases) = ($1,$2) if($tmp_I[$i]=~/([\d\.]+)\((\S+)\)/); $specific_match{$ll[3]}{$p} = "I";$specific_match_var{$ll[3]}{$p} = "$bases"; }
    }
    if($ll[18] ne 0){
	my @tmp_D = split /,/,$ll[18];
	for my $i(0..$#tmp_D){
	    my ($p,$bases) = ($1,$2) if($tmp_D[$i]=~/([\d]+)\((\S+)\)/);
	    $specific_match_var{$ll[3]}{$p} = "$bases";
	    for my $p0($p .. ($p + length($bases) - 1) ){ $specific_match{$ll[3]}{$p0} = "D"; }
	}
    }
    my ($start,$end) = $ll[12] < $ll[13] ? @ll[12,13] : @ll[13,12];
    if($ll[5] == 100){
        for my $pos($start .. $end){
	    $arg_match_status{$ll[3]}{$pos} = "M";
	    my $pos0 = $pos - 0.5;
	    (($arg_match_status{$ll[3]}{$pos0} && $arg_match_status{$ll[3]}{$pos0} eq "I") || !$arg_match_status{$ll[3]}{$pos0}) && ($arg_match_status{$ll[3]}{$pos0} = "M");
	    ($sub_end_yes || $ll[1] == $ll[6]) && ($arg_cover_exact{$ll[3]}{$pos}++);
	    $arg_cover{$ll[3]}{$pos}++;
	} 
    }else{
        for my $pos($start .. $end){
	    $arg_match_status{$ll[3]}{$pos} ||= $specific_match{$ll[3]}{$pos} ? $specific_match{$ll[3]}{$pos} : "M";
	    my $pos0 = $pos - 0.5;
	    $specific_match{$ll[3]}{$pos0} && ($arg_match_status{$ll[3]}{$pos0} ||= $specific_match{$ll[3]}{$pos0});
	    $arg_cover{$ll[3]}{$pos}++; #
	}
    }

    $arg_len{$ll[3]} = $ll[4]; #refer genome length
    my @ranks = split /;/,$ll[19];
    $arg_name{$ll[3]} = $ranks[-1]; 
    $L1_L2{$ranks[-1]} = $ranks[-2];
    $L1_L3{$ranks[-1]} = $ranks[-3];
    $anno_info{$ranks[-1]} = join "\t",@ll[20..$#ll];  #add annotation
    $hit_n_all{$ranks[-1]}++; #stat ARG-like reads number
    $lingage{$ranks[-5]}{$ranks[-4]}{$ranks[-3]}{$ranks[-2]}{$ranks[-1]}=1;

    ### LCA annotation
    my ($cid,$fann)= ("$ll[0]\t$ll[2]",$ll[19]);
    if ($lastid eq ''){
        ($lastid,$lastfull) = ($cid,$fann);
    }else{
        if($lastid ne $cid){
            my @full = split /;/,$lastfull;
            $reads_lca_anno .= join("\t",$lastid,$lastfull)."\n";
            ($lastid,$lastfull)= ($cid,$fann);
        }else{
            $lastfull = &min_same_str($lastfull,$fann);
        }
    }

    ### LCA annotation with exact match 
    $ll[5] == 100 && ($sub_end_yes || $ll[1] == $ll[6]) || next;
    $hit_n_all_exact{$ranks[-1]}++; #
    my ($cid2,$fann2)= ("$ll[0]\t$ll[2]",$ll[19]);
    if ($lastid2 eq ''){
        ($lastid2,$lastfull2) = ($cid2,$fann2);
    }else{
        if($lastid2 ne $cid2){
            my @full = split /;/,$lastfull2;
            $reads_lca_anno2 .= join("\t",$lastid2,$lastfull2)."\n";
            ($lastid2,$lastfull2)= ($cid2,$fann2);
        }else{
            $lastfull2 = &min_same_str($lastfull2,$fann2);
        }
    }

}
$reads_lca_anno .= join("\t",$lastid,$lastfull)."\n";
$reads_lca_anno2 .= join("\t",$lastid2,$lastfull2)."\n";
close IN;
open OUT,">$lcaf" || die $!;
print OUT $reads_lca_anno;
close OUT;
open OUT,">$lcaf_exact" || die $!;
print OUT "###exact blast LCA anno###\n$reads_lca_anno2";
close OUT;

### state reads number for each classification levels
my @reads_annos = split /\n/,$reads_lca_anno;
my %hit_n_lca;
for my $line(@reads_annos){
    my @aa = split /\t/,$line;
    my @ranks = split /;/,$aa[2];
    my $n = @ranks;
    for my $i(0..$#ranks){ my $j = $n-$i;$ranks[$i] eq "__" && next;$hit_n_lca{$j}{$ranks[$i]}++; }
}

my @reads_annos2 = split /\n/,$reads_lca_anno2;
my %hit_n_lca_exact;
for my $line2(@reads_annos2){
    $line2 =~ /^#/ && next;
    $line2 =~ /^\s+/ && next;
    my @aa2 = split /\t/,$line2;
    my @ranks2 = split /;/,$aa2[2];
    my $n = @ranks2;
    for my $i2(0..$#ranks2){ my $j2 = $n-$i2;$ranks2[$i2] eq "__" && next;$hit_n_lca_exact{$j2}{$ranks2[$i2]}++; }
}

### calculate coverage and depth ...
my (%hit_n_all_adjust,%arg_detect_info);
for my $arg_id(sort keys %arg_cover){
    my @POSs= sort {$a<=>$b} keys %{$arg_cover{$arg_id}};
    my $arg_coverlen_exact = keys %{$arg_cover_exact{$arg_id}};
    my ($arg_coverlen,$total_base_num,$varation);
    my %arg_cover0 = ();
    my %del = ();
    for my $i(0..$#POSs){
	if($arg_match_status{$arg_id}{$POSs[$i]} ne "D"){
	    $arg_coverlen++;
	    $arg_cover0{$POSs[$i]} = 1;
	    $total_base_num += $arg_cover{$arg_id}{$POSs[$i]};
	}
	$del{$POSs[$i]} && next;
	if($arg_match_status{$arg_id}{$POSs[$i]} eq "X"){$varation .= "$POSs[$i]($specific_match_var{$arg_id}{$POSs[$i]}),"; }
	my $pos = $POSs[$i] - 0.5;
	if($arg_match_status{$arg_id}{$pos} && $arg_match_status{$arg_id}{$pos} eq "I"){$varation .= "$pos($specific_match_var{$arg_id}{$pos}),";}
	if($arg_match_status{$arg_id}{$POSs[$i]} eq "D"){
	    $varation .= "$POSs[$i]($specific_match_var{$arg_id}{$POSs[$i]}),";
	    my $del_len = length($specific_match_var{$arg_id}{$POSs[$i]});
	    for my $p($POSs[$i] .. ($POSs[$i] + $del_len -1)){$del{$p}=1;}
	}
    }
    $varation ||= 0; $varation =~ s/,$//;
    my $discrete_blocknum = 1;
    my @POSs0 = sort {$a<=>$b} keys %arg_cover0; 
    for my $j(1..$#POSs0){$POSs[$j] > $POSs[$j-1]+1 &&  ($discrete_blocknum++);}
    my $avgDepth = sprintf "%.2f",$total_base_num/$arg_coverlen;
    my $arg_coverage = sprintf "%.3f",$arg_coverlen/$arg_len{$arg_id};
    $hit_n_all_exact{$arg_name{$arg_id}} ||= 0;
    $hit_n_lca_exact{1}{$arg_name{$arg_id}} ||= 0;
    $hit_n_all_adjust{$arg_name{$arg_id}} = sprintf "%.3f",$hit_n_all{$arg_name{$arg_id}}*($arg_coverage**3); 
    $arg_detect_info{$arg_name{$arg_id}} = "$hit_n_all{$arg_name{$arg_id}}\t$hit_n_all_adjust{$arg_name{$arg_id}}\t$hit_n_lca_exact{1}{$arg_name{$arg_id}}\t$hit_n_all_exact{$arg_name{$arg_id}}\t$discrete_blocknum\t$arg_coverlen_exact/$arg_coverlen/$arg_len{$arg_id}\t$arg_coverage\t$avgDepth\t$varation";
}

### output detected ARG profile
open OUT,">$statf" || die $!;
print OUT "#L5_Mechanism\tL5_rn_lca\tL4_family\tL4_rn_lca\tL3_subfamily\tL3_rn_lca\tL2_GeneType\tL2_rn_lca\tL1_GeneST\tL1_rn_lca\tHitRN_all\tHitRN_all_adjust\tHitRN_lca_exact\tHitRN_all_exact\tDiscrete_block_n\tCoverage\tCoverage_ratio\tCoverDepth\tVaration\tmodel_type\ttaxon_name\tEffective_enzyme_inhibitors\tgene_coded\taro_des\tAMR_Gene_Family\tAMR_Gene_Family_des\tAntibiotic\tDrug_Class\tAdjuvant\tEfflux_Component\tEfflux_Regulator\tResistance_Mechanism\tResistance_Mechanism_des\n";
my %focus_level_uniq;
for my $L5 (sort {$hit_n_lca{5}{$b} <=> $hit_n_lca{5}{$a}} keys %{$hit_n_lca{5}}){
    my @L4s = keys %{$lingage{$L5}};
    for (@L4s){$hit_n_lca{4}{$_} ||=0;}
    for my $L4 (sort {$hit_n_lca{4}{$b} <=> $hit_n_lca{4}{$a}} @L4s){  #
	my @L3s = keys %{$lingage{$L5}{$L4}};
	for (@L3s){$hit_n_lca{3}{$_} ||=0;}
        my @L1s = ();
 	for my $L3 (@L3s){
	    my @L2s = keys %{$lingage{$L5}{$L4}{$L3}};
	    for (@L2s){$hit_n_lca{2}{$_} ||=0;}
    	    for my $L2(@L2s){ 
    	        my @L1s_tmp = keys %{$lingage{$L5}{$L4}{$L3}{$L2}}; 
    	        for (@L1s_tmp){ $hit_n_all{$_} || next;push @L1s,$_; }
    	    }
        }
        for my $L1 (sort {$hit_n_all_adjust{$b} <=> $hit_n_all_adjust{$a}} @L1s){  #
            my ($L2,$L3) = ($L1_L2{$L1},$L1_L3{$L1});
    	    $hit_n_lca{1}{$L1} ||= 0;$hit_n_lca{2}{$L2} ||= 0;$hit_n_lca{3}{$L3} ||= 0;$hit_n_lca{4}{$L4} ||= 0;
            if($hit_n_lca_exact{1}{$L1}){
	        $focus_level_uniq{3}{$L3}++;$focus_level_uniq{4}{$L4}++;$focus_level_uniq{5}{$L5}++;
	        print OUT "$L5\t$hit_n_lca{5}{$L5}\t$L4\t$hit_n_lca{4}{$L4}\t$L3\t$hit_n_lca{3}{$L3}\t$L2\t$hit_n_lca{2}{$L2}\t$L1\t$hit_n_lca{1}{$L1}\t$arg_detect_info{$L1}\t$anno_info{$L1}\n";
	    }elsif($hit_n_lca{2}{$L2}){
	        $focus_level_uniq{3}{$L3}++;$focus_level_uniq{4}{$L4}++;$focus_level_uniq{5}{$L5}++;
	        $focus_level_uniq{3}{$L3} > 5 && next;
	        print OUT "$L5\t$hit_n_lca{5}{$L5}\t$L4\t$hit_n_lca{4}{$L4}\t$L3\t$hit_n_lca{3}{$L3}\t$L2\t$hit_n_lca{2}{$L2}\t$L1\t0\t$arg_detect_info{$L1}\t$anno_info{$L1}\n";
	    }elsif($hit_n_lca{3}{$L3}){
	        $focus_level_uniq{3}{$L3}++;$focus_level_uniq{4}{$L4}++;$focus_level_uniq{5}{$L5}++;
	        $focus_level_uniq{3}{$L3} > 5 && next;
	        print OUT "$L5\t$hit_n_lca{5}{$L5}\t$L4\t$hit_n_lca{4}{$L4}\t$L3\t$hit_n_lca{3}{$L3}\t$L2\t0\t$L1\t0\t$arg_detect_info{$L1}\t$anno_info{$L1}\n";
	    }elsif($hit_n_lca{4}{$L4}){
	        $focus_level_uniq{4}{$L4}++;$focus_level_uniq{5}{$L5}++;
	        $focus_level_uniq{4}{$L4} > 5 && next;
	        print OUT "$L5\t$hit_n_lca{5}{$L5}\t$L4\t$hit_n_lca{4}{$L4}\t$L3\t0\t$L2\t0\t$L1\t0\t$arg_detect_info{$L1}\t$anno_info{$L1}\n";
            }elsif($hit_n_lca{5}{$L5}){
	        $focus_level_uniq{5}{$L5}++;
	        $focus_level_uniq{5}{$L5} > 5 && next;
	        print OUT "$L5\t$hit_n_lca{5}{$L5}\t$L4\t0\t$L3\t0\t$L2\t0\t$L1\t0\t$arg_detect_info{$L1}\t$anno_info{$L1}\n";
            }
	}
    } 
}
close OUT;

####################sub function#########################
sub min_same_str {
    my ($s1,$s2) = @_;
    my @a1 = split /;/,$s1;
    my @a2 = split /;/,$s2;
    my $min = @a1 < @a2 ? @a1 : @a2;
    my @same;
    for my $i (0..$min-1) {
        if($a1[$i] ne $a2[$i]){
                push @same,"__";
        }else{
                push @same,$a1[$i];
        }
    }
    my $same_str = join(";",@same);
    $same_str ||= "";
    return $same_str;
}

