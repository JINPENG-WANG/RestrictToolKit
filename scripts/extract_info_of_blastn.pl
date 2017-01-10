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
my $filename="blastn_alignments_of_${ref}_by_${front_enzyme}_and_${later_enzyme}.txt";
open BLASTRLT,"$path/run_results/$result_dir/$filename" or die "cannot open $filename:$!";
#open EXTRACTRLT,">>$path/run_results/selected_$filename" or die "cannot write to selected_${filename}:$!"; 
open FIRSTLINE,"$path/run_results/$result_dir/$filename" or die "cannot open $filename:$!";
open MULTIALIGNCOUNTS,">>$path/run_results/$result_dir/selected_${filename}"or die "cannot write to selected_${filename}:$!";
my $frag_serial_num=1;
my $frag_name;
my $frag_length;
my ($querry,$subject,$length,$length_down);
foreach(<FIRSTLINE>){
my $first_line=$_;
chomp $first_line;
my @first_line=split /\t/, $first_line;
$frag_name=$first_line[0];
$frag_length=$first_line[3];
last;
}
my $multialigncounts=0;
foreach(<BLASTRLT>){
my $line=$_;
chomp $line;
my @line=split /\t/,$line;
$querry=$line[0];
$subject=$line[1];
$length=$line[3];
$length_down=$frag_length*0.8;
if ($querry eq $frag_name){
if ($length > $length_down){
#print EXTRACTRLT "$line\n"; 
$multialigncounts++;
}
}
unless ( $querry eq $frag_name){
#print"$frag_serial_num $frag_name is done\n";
if ($multialigncounts > 0){
print MULTIALIGNCOUNTS"$frag_name\t$multialigncounts\n";
}
if ($multialigncounts ==0){
print MULTIALIGNCOUNTS"$frag_name\t1\n";
}
$multialigncounts=0;
$frag_name=$querry;
$frag_length=$line[3];
#print EXTRACTRLT "$line\n";
$frag_serial_num++;
}
}


