#!usr/bin/perl -w
use strict;

(@ARGV==3)||die "perl $0 <fq.file> <fa.file>\n";

my ($infile,$outfile,$pe)=@ARGV;  
$infile =~ /gz$/ ? open IN,"gzip -dc $infile |" : open IN,$infile;
open OUT,">$outfile";
while ( my $id=<IN>){
    chomp $id;
    my $id2= (split /\s+/,$id)[0];
    $id2 =~ s/^\@//;
    chomp (my $seq=<IN>);
    my $len = length ($seq);
    <IN>;
    <IN>;
    print OUT ">$id2\_$pe\n$seq\n";
}

close IN;
close OUT;
