#!/usr/bin/perl -w
use strict;

@ARGV || die "perl $0 <detect_arg.anno> > candidate_feature.stat\n\n";
my $file = shift @ARGV;
open IN,$file || die $!;
my (%model_type,%copyn,%var);
while(<IN>){
    chomp;
    my @ll = split /\t/;
    $model_type{$ll[4]} = $ll[5];
    $copyn{$ll[4]}++;
    $ll[6] ne "0" && ($var{$ll[4]} .= "$ll[6],");
    $ll[7] ne "0" && ($var{$ll[4]} .= "$ll[7],");
    $ll[8] ne "0" && ($var{$ll[4]} .= "$ll[8],");
}
close IN;

my %uniq;
for my $arg(sort keys %copyn){
    if($var{$arg}){
        $var{$arg} =~ s/,$//;
	my @aa = split /,/,$var{$arg};
	my $vars;
	for my $v(@aa){
	    my $k = "$arg:$v";
	    $uniq{$k}++;
	    $uniq{$k} > 1 && next;
	    $vars .= "$v,";
	}
	$vars =~ s/,$//;
	print "$arg\t$copyn{$arg}\t$model_type{$arg}\t$vars\n";
    }else{
	print "$arg\t$copyn{$arg}\t$model_type{$arg}\t0\n";
    }
} 

