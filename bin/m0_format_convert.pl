#!/usr/binperl -w
use strict;

@ARGV == 2 || die "Usage: perl $0 <blastn.m0> <out:blast.m6>
Description: a script to covert m0 to m6 and obtain mutants info\n\n";

my ($inf,$outf) = @ARGV;
open IN,$inf || die $!;
open OUT,">$outf" || die $!;
$/= "Query= ";<IN>;
while(my $query0 = <IN>){
    $query0 =~ /\* No hits found \*/ && next;
    chomp $query0;
    my @blocks = split />/,$query0;
    my @overall = split /\n/,$blocks[0];
    my $query_id = (split /\s+/,$overall[0])[0];
    my $query_len = $1 if($blocks[0]=~/Length=([\d]+)/);

    for my $i(1..$#blocks){
	my @blocks2 = split /Score = /,$blocks[$i];
	my $sub_id = $1 if($blocks2[0]=~/^(\S+)/);
	my $sub_len = $1 if($blocks2[0] =~ /Length=([\d]+)/);
	for my $m(1..$#blocks2){
	    my ($sub_score,$sub_evalue) = ($1,$2) if($blocks2[$m] =~ /^([\d\.]+) bits\s+\([\d]+\),\s+Expect = (\S+)/);
	    my ($sub_ident_ratio,$sub_gap)=($1,$2) if($blocks2[$m] =~ /Identities =\s+(.+?)\s+\(\S+\),\s+Gaps =\s+([\d]+)\/.+?\)/);
	    my ($nident,$align_len) = split /\//,$sub_ident_ratio;
	    my $sub_ident = sprintf "%.3f",($nident/$align_len)*100;
	    my $strand = $1 if($blocks2[$m]=~/Strand=(\S+)/);

	    my @hits = split /\n/,$blocks2[$m];
	    my ($query_start,$query_end,$sub_start,$sub_end);
	    my $mis = "";
	    my %ins = ();
	    my %del = ();
	    for my $j(3..$#hits){
	        $hits[$j]=~ /^Query/ || next;
	        my @que = split /\s+/,$hits[$j];
	        my @sub = split /\s+/,$hits[$j+2];
	        my @aa = split //,$que[2];
	        my @bb = split //,$sub[2];
	        $query_start ||= $que[1];$query_end = $que[3];
	        $sub_start ||= $sub[1];$sub_end = $sub[3];
	        for my $k(0..$#bb){
		    $aa[$k] = uc($aa[$k]);
		    $bb[$k] = uc($bb[$k]);
		    if($aa[$k] ne $bb[$k]){
		        my $seq0 = join("",@bb[0..$k]);
		        $seq0 =~ s/-//g;
		        my $var_pos = $strand=~ /\/Minus/ ? ($sub[1] - length($seq0) +1) : ($sub[1] + length($seq0) - 1);
			$strand=~ /\/Minus/ && ($bb[$k] =~ tr/AGCTagct/TCGAtcga/ && $aa[$k] =~ tr/AGCTagct/TCGAtcga/);
		        if($aa[$k] ne "-"){
			    if($bb[$k] ne "-"){
			        $mis .= "$var_pos($bb[$k]\->$aa[$k]),";
			    }else{
				$var_pos = $strand=~ /\/Minus/ ? ($var_pos - 0.5) : ($var_pos + 0.5);
			        $ins{$var_pos} .= "$aa[$k]";
			    }
		        }else{
			    $del{$var_pos} = $bb[$k];
		        }
		    }
	        }
	    }
	    $mis =~ s/,$//;$mis ||= 0;
	    my $mis_n= $mis eq 0 ? 0 : (split /,/,$mis);
	    my $inss = "";
	    my @strands = split /\//,$strand;
	    for my $k(sort { $a<=>$b } keys %ins){
		if($strands[0] eq $strands[1]){
		    $inss .= "$k($ins{$k}),";
		}else{
		    $ins{$k} =~ tr/AGCT/TCGA/;$ins{$k} = reverse($ins{$k});
		    $inss .= "$k($ins{$k}),";
		}
	    }
	    $inss =~ s/,$//;
	    $inss ||= 0;

	    my @del_pos = sort { $a<=>$b } keys %del;
	    my $dels = "";
	    if(@del_pos > 1){
	        for my $i(1..$#del_pos){
		    my $j = $i - 1;
		    if($del_pos[$i] == ($del_pos[$j]+1)){
		        $del{$del_pos[$i]} .= $del{$del_pos[$j]};
		        delete $del{$del_pos[$j]};
		    }
	        }
	        for my $k(sort {$a<=>$b} keys %del){
		    my $p = $k - length($del{$k}) + 1;
		    my $bases = reverse $del{$k};
		    $dels .= "$p($bases),";
	        }
	        $dels =~ s/,$//;
	    }elsif(@del_pos == 1){
		$dels = "$del_pos[0]($del{$del_pos[0]})";
	    }
	    $dels ||= 0;

	    print OUT join("\t",$query_id,$query_len,$sub_id,$sub_len,$sub_ident,$align_len,$nident,$mis_n,$sub_gap,$query_start,$query_end,$sub_start,$sub_end,$sub_evalue,$sub_score,$mis,$inss,$dels),"\n";
        }
    }
}
$/= "\n";
close IN;
close OUT;

