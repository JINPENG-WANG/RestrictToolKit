#!/usr/bin/perl -w
use strict;
my ($path,$ref,$front_enzyme,$later_enzyme,$result_dir);
while(@ARGV){
$_=shift @ARGV;
if($_=~/^-p$/){$path=shift @ARGV;}
elsif($_=~/^-ref$/){$ref=shift @ARGV;}
elsif($_=~/^-e1$/){$front_enzyme=shift @ARGV;}
elsif($_=~/^-e2$/){$later_enzyme=shift @ARGV;}
elsif($_=~/^-dir$/){$result_dir=shift @ARGV;}
}
my $filename="selected_blastn_alignments_of_${ref}_by_${front_enzyme}_and_${later_enzyme}.txt";
open FRAGMULTICOUNTS,"$path/run_results/$result_dir/$filename" or die "cannot open $filename:$!";
open STATISTICS,">>$path/run_results/$result_dir/characteristics_of_${ref}_by_${front_enzyme}_and_${later_enzyme}.txt"or die "cannot write to characteristics_of_${ref}_by_${front_enzyme}_and_${later_enzyme}.txt:$!";
my $line;
my @line;
my $frag_name;
my $multi_align_counts;
my $all_frag_nums=0;
my $frag_nums_unique=0;
my $frag_nums_more_2_5=0;
my $frag_nums_more_6_10=0;
my $frag_nums_more_11_20=0;
my $frag_nums_more_21_50=0;
my $frag_nums_more_51_100=0;
my $frag_nums_more_100=0;
my $rate_unique=0;
my $rate_more_2_5=0;
my $rate_more_6_10=0;
my $rate_more_11_20=0;
my $rate_more_21_50=0;
my $rate_more_51_100=0;
my $rate_more_100=0;
foreach(<FRAGMULTICOUNTS>){
$line=$_;
chomp $line;
@line=split /\t/,$line;
$frag_name=$line[0];
$multi_align_counts=$line[1];
$all_frag_nums++;
if($multi_align_counts ==1){
$frag_nums_unique++;
}
if(1<$multi_align_counts && $multi_align_counts< 6){
$frag_nums_more_2_5++;
}
if(5 <$multi_align_counts && $multi_align_counts< 11){
$frag_nums_more_6_10++;
}
if(10<$multi_align_counts && $multi_align_counts< 21){
$frag_nums_more_11_20++
}
if(20<$multi_align_counts && $multi_align_counts<51){
$frag_nums_more_21_50++;
}
if(50<$multi_align_counts && $multi_align_counts<101){
$frag_nums_more_51_100++;
}
if(100<$multi_align_counts){
$frag_nums_more_100++;
}
}
$rate_unique=$frag_nums_unique/$all_frag_nums;
$rate_more_2_5=$frag_nums_more_2_5/$all_frag_nums;
$rate_more_6_10=$frag_nums_more_6_10/$all_frag_nums;
$rate_more_11_20=$frag_nums_more_11_20/$all_frag_nums;
$rate_more_21_50=$frag_nums_more_21_50/$all_frag_nums;
$rate_more_51_100=$frag_nums_more_51_100/$all_frag_nums;
$rate_more_100=$frag_nums_more_100/$all_frag_nums;
print "There are $frag_nums_unique fragments that can align reference uniquely\n";
printf STATISTICS "GG\tunique alignment rate is\t%1.4f\n",$rate_unique;
