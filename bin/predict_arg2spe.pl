#!/usr/bin/perl -w
use strict;

@ARGV || die "Usage:perl $0 <reads_arg.lca.anno> <reads_tax.anno> >arg_spe.xls
Version: 1.0,20210316
Contact: p.han\@genskey.com\n\n";

my ($arg_anno_lca,$tax_anno)= @ARGV;
open IN,$arg_anno_lca || die $!;
my (%reads2arg,%arg_readsnum,%args);
while(<IN>){
    chomp;
    my @ll = split /\t/;
    $ll[0] =~ s/_1$//;
    $ll[0] =~ s/_2$//;
    my ($L3,$L2) = (split /;/,$ll[2])[-3,-2];
    $L3 eq "__" && $L2 eq "__" && next;
    $arg_readsnum{$L3}++; 
    $arg_readsnum{"$L3\t$L2"}++;  
    $args{"$L3\t$L2"} = 1;
    $reads2arg{$ll[0]} = "$L3\t$L2"; 
}
close IN;

### 
open IN,$tax_anno || die $!;
my (%taxon_rn,%taxon_rn2,%uniq,%spe_complex,%spe_genus,%complex_genus);
my ($genus_index,$complex_index,$spe_index);
my $header = <IN>;chomp $header;
my @head = split /\t/,$header;
for my $i(0..$#head){
    $head[$i] eq "Genus" && ($genus_index = $i);
    $head[$i] eq "Complex" && ($complex_index = $i);
    $head[$i] eq "Species" && ($spe_index = $i);
}
while(<IN>){
    /^#/ && next;
    chomp;
    my @ll = split /\t/;
    my ($genus,$complex,$spe) = @ll[3,4,5];
    if($reads2arg{$ll[0]}){
	my $L3 = (split /\t/,$reads2arg{$ll[0]})[0];
	$uniq{$ll[0]}{$spe}++;
	$uniq{$ll[0]}{$spe} > 1 && next;
	$taxon_rn2{$L3}{"spe"}{$spe}++;
        if($complex ne "-" ){
            $spe_complex{$spe} = $complex;	
	    $uniq{$ll[0]}{$complex}++;
	    $uniq{$ll[0]}{$complex} > 1 && next;
	    $taxon_rn2{$L3}{"complex"}{$complex}++;
        }
	if($genus ne "-"){
	    $spe_genus{$spe} = $genus;
	    $uniq{$ll[0]}{$genus}++;
            $uniq{$ll[0]}{$genus} > 1 && next;
	    $taxon_rn2{$L3}{"genus"}{$genus}++;
	}
    }
}
close IN;

### Output
my %genetype_detect;
print "#L3_GF\tL1_GeneType\tgenus\tspe_complex\tspe\n";
for my $arg(sort keys %args){
    my ($L3,$L1) = split /\t/,$arg;
    $L1 ne "__" && ($genetype_detect{$L3} = 1);
}
my $ratio_cutoff = 0.2;
for my $arg(sort keys %args){
    my ($L3,$L1) = split /\t/,$arg;
    $L1 eq "__" && $genetype_detect{$L3} && next; 
    my @spes = sort {$taxon_rn2{$L3}{"spe"}{$b} <=> $taxon_rn2{$L3}{"spe"}{$a}} keys %{$taxon_rn2{$L3}{"spe"}};
    my @complexs = sort {$taxon_rn2{$L3}{"complex"}{$b} <=> $taxon_rn2{$L3}{"complex"}{$a}} keys %{$taxon_rn2{$L3}{"complex"}};
    my @genuss = sort {$taxon_rn2{$L3}{"genus"}{$b} <=> $taxon_rn2{$L3}{"genus"}{$a}} keys %{$taxon_rn2{$L3}{"genus"}};
    my ($spe_ratio,$complex_ratio,$genus_ratio,$spe_r,$complex_r,$genus_r);
    if(@spes){
	for my $i(0..$#spes){
	    $i <= 10 || next;
	    $spe_r = $taxon_rn2{$L3}{spe}{$spes[$i]}/$arg_readsnum{$L3};
	    $spe_r >= $ratio_cutoff && ($spe_ratio .= "$spes[$i]($taxon_rn2{$L3}{spe}{$spes[$i]}/$arg_readsnum{$L3})|");
	    if($spe_complex{$spes[$i]}){
		$complex_r = $taxon_rn2{$L3}{complex}{$spe_complex{$spes[$i]}}/$arg_readsnum{$L3};
	        $complex_r >= $ratio_cutoff && ($complex_ratio .= "$spe_complex{$spes[$i]}($taxon_rn2{$L3}{complex}{$spe_complex{$spes[$i]}}/$arg_readsnum{$L3})|");
	    }else{
	        $complex_ratio .= "-|";
  	    }
	    if($spe_genus{$spes[$i]}){
		$genus_r = $taxon_rn2{$L3}{genus}{$spe_genus{$spes[$i]}}/$arg_readsnum{$L3};
	        $genus_r >= $ratio_cutoff && ($genus_ratio .= "$spe_genus{$spes[$i]}($taxon_rn2{$L3}{genus}{$spe_genus{$spes[$i]}}/$arg_readsnum{$L3})|");
	    }else{
	        $genus_ratio .= "-|";
	    }
	}
	$spe_ratio ||= "$spes[0]($taxon_rn2{$L3}{spe}{$spes[0]}/$arg_readsnum{$L3})";
	$complex_ratio ||= "-";
	$genus_ratio ||= "-";
    }elsif(@complexs){
	$spe_ratio = "-";
	for my $i(0..$#complexs){
	    $i <= 10 || next;
	    $complex_r = $taxon_rn2{$L3}{complex}{$complexs[$i]}/$arg_readsnum{$L3};
	    $complex_r >= $ratio_cutoff && ($complex_ratio .= "$complexs[$i]($taxon_rn2{$L3}{complex}{$complexs[$i]}/$arg_readsnum{$L3})|");
	    if($complex_genus{$complexs[$i]}){
		$genus_r = $taxon_rn2{$L3}{genus}{$complex_genus{$complexs[$i]}}/$arg_readsnum{$L3};
	        $genus_r >= $ratio_cutoff && ($genus_ratio .= "$complex_genus{$complexs[$i]}($taxon_rn2{$L3}{genus}{$complex_genus{$complexs[$i]}}/$arg_readsnum{$L3})|");
	    }else{
	        $genus_ratio .= "-|";
	    }
	}
	$complex_ratio ||= "-";
	$genus_ratio ||= "-";
    }elsif(@genuss){
	$spe_ratio = "-";
	$complex_ratio = "-";
	for my $i(0..$#complexs){
	    $i <= 10 || next;
	    $genus_r = $taxon_rn2{$L3}{genus}{$genuss[$i]}/$arg_readsnum{$L3};
	    $genus_r >= $ratio_cutoff && ($genus_ratio .= "$genuss[$i]($taxon_rn2{$L3}{genus}{$genuss[$i]}/$arg_readsnum{$L3})|");
	}
	$genus_ratio ||= "-";
    }else{
	$spe_ratio = "-";
	$complex_ratio = "-";
	$genus_ratio = "-";
    }
    my (%unique,$com_ratio,$gen_ratio);
    my @com = split /\|/,$complex_ratio;
    for (@com){$_ eq "-" && next;$unique{$_}++;$unique{$_} < 2 && ($com_ratio .= "$_|");}
    my @gen = split /\|/,$genus_ratio;
    for (@gen){$_ eq "-" && next;$unique{$_}++;$unique{$_} < 2 && ($gen_ratio .= "$_|");}
    $gen_ratio ||= "-";
    $com_ratio ||= "-";
    $com_ratio =~ s/\|$//;
    $gen_ratio =~ s/\|$//;
    $spe_ratio =~ s/\|$//;
    print "$arg\t$gen_ratio\t$com_ratio\t$spe_ratio\n";
}

