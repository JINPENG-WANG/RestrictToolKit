#!/usr/bin/perl -w 
use strict;
#input the filename behind parameter -ref; the later enzyme behind -e2; the front enzyme behind -e1;
my @input=@ARGV;
my($filename,$front_enzyme,$later_enzyme,$path,$result_dir);
my $range1=90; my $range2=750;
while(@ARGV){
$_=shift @ARGV;
if($_=~/^-ref$/){$filename=shift @ARGV;}
elsif($_=~/^-e1$/){$front_enzyme=shift @ARGV;}
elsif($_=~/^-e2$/){$later_enzyme=shift @ARGV;}
elsif($_=~/^-p$/){$path=shift @ARGV;}
elsif($_=~/^-r1$/){$range1=shift @ARGV;}
elsif($_=~/^-r2$/){$range2=shift @ARGV;}
elsif($_=~/^-dir$/){$result_dir=shift @ARGV;}
}
#set the paths of input files and output files;
open SCAFD,"$path/reference/$filename" or die "cannot open $filename :$!";
open LOCSINRG,">>$path/run_results/$result_dir/restrict_locs_in_range_of_${filename}_by_${front_enzyme}_and_${later_enzyme}.txt" or die "cannot write to restrict_locs_in_range_of_${filename}_by_${front_enzyme}_and_${later_enzyme}.txt:$!";
open RATE,">>$path/run_results/$result_dir/characteristics_of_${filename}_by_${front_enzyme}_and_${later_enzyme}.txt" or die "cannot write to characteristics_of_${filename}_by_${front_enzyme}_and_${later_enzyme}.txt:$!"; 
open SEQFRAGINRANGE,">>$path/run_results/$result_dir/seq_of_frags_in_range_of_${filename}_by_${front_enzyme}_and_${later_enzyme}.txt" or die "cannot write to seq_of_frags_in_range_of_${filename}_by_${front_enzyme}_and_${later_enzyme}.txt:$!";
#put the scaffold names and their sequences into the array @scafds_and_seqs;
my @scafds_and_seqs;
foreach(<SCAFD>){
chomp;
push @scafds_and_seqs,$_;
}
#input the contents of enzyme.list into %enzyme_list;
my %enzyme_list=(
EcoRI => 'G|AATTC',
HinfI => 'G|ANTC',
AluI => 'AG|CT',
AvaII => 'G|GWCC',
BamHI => 'G|GATCC',
HindIII => 'A|AGCTT',
MseI => 'T|TAA',
MspI => 'C|CGG',
NdeI => 'CA|TATG',
PstI => 'CTGCA|G',
SalI => 'G|TCGAC',
XbaI => 'T|CTAGA',
XhoI => 'C|TCGAG',
NlaIII => 'NCATG|N',
);
my ($front_seq,$later_seq);
foreach(keys %enzyme_list){
if($_=~/$front_enzyme/i){
 $front_seq=$enzyme_list{$_};
}
if($_=~/$later_enzyme/i){
 $later_seq=$enzyme_list{$_};
}
}
#get the relative location of cutting of the recognation sequences; $front_enzyme_cut_loc and $later_enzyme_cut_loc are two relative locations we wanted;
my @front_enzyme_rec_seq=split //, $front_seq;
my @later_enzyme_rec_seq=split //, $later_seq;
my @front_enzyme_cut_seq;
my @later_enzyme_cut_seq;
my $front_enzyme_rec_start=0;
my $later_enzyme_rec_start=0;
my $front_enzyme_cut_loc;
my $later_enzyme_cut_loc;
foreach(@front_enzyme_rec_seq){
$front_enzyme_rec_start++;
if($_=~/\|/){
$front_enzyme_cut_loc=$front_enzyme_rec_start;
}
else{
push @front_enzyme_cut_seq,$_;
}
}
foreach(@later_enzyme_rec_seq){
$later_enzyme_rec_start++;
if($_=~/\|/){
$later_enzyme_cut_loc=$later_enzyme_rec_start;
}
else{
push @later_enzyme_cut_seq, $_;
}
}
#define the degenerate bases of front enzyme;
foreach(@front_enzyme_cut_seq){
if($_=~/N/){$_="A|T|C|G";}
if($_=~/R/){$_="A|G";}
if($_=~/Y/){$_="C|T";}
if($_=~/M/){$_="A|C";}
if($_=~/K/){$_="G|T";}
if($_=~/S/){$_="G|C";}
if($_=~/W/){$_="A|T";}
if($_=~/H/){$_="A|T|C";}
if($_=~/B/){$_="G|T|C";}
if($_=~/V/){$_="G|A|C";}
if($_=~/D/){$_="G|A|T";}
}
my($front_fir_base,$front_sec_base,$front_thi_base,$front_for_base,$front_fif_base,$front_six_base,$front_sev_base,$front_eig_base)=@front_enzyme_cut_seq;
my $num_of_front_enzyme_cut_seq=@front_enzyme_cut_seq;
#define the degenerate bases of later enzyme;
foreach(@later_enzyme_cut_seq){
if($_=~/N/){$_="A|T|C|G";}
if($_=~/R/){$_="A|G";}
if($_=~/Y/){$_="C|T";}
if($_=~/M/){$_="A|C";}
if($_=~/K/){$_="G|T";}
if($_=~/S/){$_="G|C";}
if($_=~/W/){$_="A|T";}
if($_=~/H/){$_="A|T|C";}
if($_=~/B/){$_="G|T|C";}
if($_=~/V/){$_="G|A|C";}
if($_=~/D/){$_="G|A|T";}
}
my($later_fir_base,$later_sec_base,$later_thi_base,$later_for_base,$later_fif_base,$later_six_base,$later_sev_base,$later_eig_base)=@later_enzyme_cut_seq;
my $num_of_later_enzyme_cut_seq=@later_enzyme_cut_seq;
#define the scalars and arrays used in the whole programme;
my $rate_overall=0;
my $rate_in_range=0;
my $seq_num=1;
my $lines_num=@scafds_and_seqs;
my $scafds_num=$lines_num/2;
my @length_of_fragment;
my @fragments_in_range;
my $overall_fragments_length=0;
my $PE100_overall_fragments_length=0;
my $PE100_overall_fragments_in_range_length=0;
my $overall_fragments_in_range_length=0;
my $overall_length_of_scfd=0;

########## CIRCLES START HERE!###########

#handle every scaffold and sequence of it; they are contained in the array @scafds_and_seqs;
while($seq_num<$lines_num){
my $scafd_num=$seq_num-1;
my $scafd_name=$scafds_and_seqs[$scafd_num];
my @bases=split //, $scafds_and_seqs[$seq_num];
my $length_of_scafd=@bases;
$seq_num=$seq_num+2;
my $num_compared=@bases-$num_of_front_enzyme_cut_seq-1;
#set the arrays which contain cut locations of front enzyme and later enzyme;
my @front_enzyme_locs;
my $front_enzyme_start=0;
my @later_enzyme_locs;
my $later_enzyme_start=0;
#find locs of recognition of front enzyme in this scaffold;the locs are meaningful in this array; not in the human context;
while($front_enzyme_start<$num_compared){
if($bases[$front_enzyme_start]=~/$front_fir_base/){
if($bases[$front_enzyme_start+1]=~/$front_sec_base/){
if($bases[$front_enzyme_start+2]=~/$front_thi_base/){
if($bases[$front_enzyme_start+3]=~/$front_for_base/){
unless ($front_fif_base){
push @front_enzyme_locs, $front_enzyme_start;
}
if($front_fif_base){
if($bases[$front_enzyme_start+4]=~/$front_fif_base/){
unless($front_six_base){
push @front_enzyme_locs, $front_enzyme_start;
}
if($front_six_base){
if($bases[$front_enzyme_start+5]=~/$front_six_base/){
unless($front_sev_base){
push @front_enzyme_locs,$front_enzyme_start;
}
if($front_sev_base){
if($bases[$front_enzyme_start+6]=~/$front_sev_base/){
unless($front_eig_base){
push @front_enzyme_locs,$front_enzyme_start;
}
if($bases[$front_enzyme_start+7]=~/$front_eig_base/){
push @front_enzyme_locs,$front_enzyme_start;
}
}
}
}
}
}
}
}
}
}
}
$front_enzyme_start++;
} 

#find locs of later enzyme in this scaffold;
while($later_enzyme_start<$num_compared){
if($bases[$later_enzyme_start]=~/$later_fir_base/){
if($bases[$later_enzyme_start+1]=~/$later_sec_base/){
if($bases[$later_enzyme_start+2]=~/$later_thi_base/){
if($bases[$later_enzyme_start+3]=~/$later_for_base/){
unless ($later_fif_base){
push @later_enzyme_locs, $later_enzyme_start;
}
if($later_fif_base){
if($bases[$later_enzyme_start+4]=~/$later_fif_base/){
unless($later_six_base){
push @later_enzyme_locs, $later_enzyme_start;
}
if($later_six_base){
if($bases[$later_enzyme_start+5]=~/$later_six_base/){
unless($later_sev_base){
push @later_enzyme_locs,$later_enzyme_start;
}
if($later_sev_base){
if($bases[$later_enzyme_start+6]=~/$later_sev_base/){
unless($later_eig_base){
push @later_enzyme_locs,$later_enzyme_start;
}
if($bases[$later_enzyme_start+7]=~/$later_eig_base/){
push @later_enzyme_locs,$later_enzyme_start;
}
}
}
}
}
}
}
}
}
}
}
$later_enzyme_start++;
} 
#integrate the front enzyme locs and later enzyme locs into a hash;
my %hash;
foreach(@front_enzyme_locs){
$hash{$_}=$front_enzyme;
}

foreach(@later_enzyme_locs){
$hash{$_}=$later_enzyme;
}
#sort the locs as the locations of them;
my @locs_sorted=sort{$a<=>$b} keys %hash;
my @locs_wanted_before;
my @enzyme_wanted_before;
my @locs_wanted_after;
my @enzyme_wanted_after;
my $locs_sorted_start=0;
while ($locs_sorted_start<@locs_sorted-1){
my $a=$locs_sorted_start;
my $b=$locs_sorted_start+1;
my $c=$locs_sorted[$a];
my $d=$locs_sorted[$b];
my $e=${hash}{$c};
my $f=${hash}{$d};
#$c and $d are locations in the array which starts from 0, to make it more convenient to read for human, $c+1 and $d+1 is given;
#filter the locs and put the fragments with different terminals into new arrays; these locs are also start with 0;
if( $e !~ $f){
push @locs_wanted_before, $c;
push @enzyme_wanted_before, $e;
push @locs_wanted_after, $d;
push @enzyme_wanted_after, $f;
#LOCS contains the locs of recognition ,not the cut locations; it is for human convenience;
}
$locs_sorted_start++;
}
#define the scalars which contains the length of fragments;
my $nums_of_fragments=@locs_wanted_before;
my @pairs_of_locs=(0..$nums_of_fragments-1);
my $length_of_fragments=0;
my $length_of_fragments_in_range=0;
my $PE100_length_of_fragments=0;
my $PE100_length_of_fragments_in_range=0;
my $count_all=0;
my $count_in_range=0;
my ($before_cut_coorr,$after_cut_coorr,$bf_ct_crr_hm_rd,$af_ct_crr_hm_rd);
foreach(@pairs_of_locs){
$count_all++;
my $before_coorr=${locs_wanted_before}[$_];
my $before_enzyme=${enzyme_wanted_before}[$_];
my $after_coorr=${locs_wanted_after}[$_];
my $after_enzyme=${enzyme_wanted_after}[$_];
my $seq_selected;
my $length_of_fragment;
my $bf_ct_crr;my $af_ct_crr;
#attentions here!!! the output sequences of fragments are in accordance with the results of illumian sequencing.For removing the sequencing primers, read1 and read2 can align with the fragments here directly!!! That is to say, fragments here contain the restriction overhangs!!!
if($front_enzyme=~/$before_enzyme/i){
$before_cut_coorr=$before_coorr+$front_enzyme_cut_loc-1;
$after_cut_coorr=$after_coorr+$num_of_later_enzyme_cut_seq-$later_enzyme_cut_loc;
$length_of_fragment=$after_cut_coorr-$before_cut_coorr+1;
$seq_selected=join "", @bases[$before_cut_coorr..$after_cut_coorr];
}
if($later_enzyme=~/$before_enzyme/i){
$before_cut_coorr=$before_coorr+$later_enzyme_cut_loc-1;
$after_cut_coorr=$after_coorr+$num_of_front_enzyme_cut_seq-$front_enzyme_cut_loc;
$length_of_fragment=$after_cut_coorr-$before_cut_coorr+1;
$seq_selected=join "", @bases[$before_cut_coorr..$after_cut_coorr];
}
$bf_ct_crr=$before_cut_coorr+1;
$af_ct_crr=$after_cut_coorr+1;
#the array "@length_of_fragment" contains the lengths of all fragments;
push @length_of_fragment, $length_of_fragment;
#we treat the fragments in range in the same way of all fragments as above;
if($length_of_fragment>($range1-1)){
	if($length_of_fragment<($range2+1)){
$count_in_range++;
print SEQFRAGINRANGE "$scafd_name-$count_all\n$seq_selected\n";
print LOCSINRG "$scafd_name-$count_all\t$before_enzyme\t$bf_ct_crr\t$after_enzyme\t$af_ct_crr\t$length_of_fragment\n";
$length_of_fragments_in_range+=$length_of_fragment;
	}}
#count all fragments length no matter in range or not!
$length_of_fragments+=$length_of_fragment;
}
my $reduced_rate=$length_of_fragments/$length_of_scafd;
my $fragment_in_range_rate_scafd=$length_of_fragments_in_range/$length_of_scafd;
$overall_fragments_length+=$length_of_fragments;
$overall_fragments_in_range_length+=$length_of_fragments_in_range;
$overall_length_of_scfd+=$length_of_scafd;
#print"$scafd_name is done\n"; 
}
my $overall_rate=$overall_fragments_length/$overall_length_of_scfd;
my $overall_rate_in_range=$overall_fragments_in_range_length/$overall_length_of_scfd;
#calculate the ratio of the number of fragments which are smaller than $range1; and the ratio of the number of fragments which between $range1 and $range2;
my $num_all_fragments=@length_of_fragment;
my $num_fragments_little_range1=0;
my $num_fragments_in_range=0;
foreach(@length_of_fragment){
if($_<$range1){$num_fragments_little_range1++;}
elsif($_<=$range2){$num_fragments_in_range++;}
}
my $num_rate_little_range1=$num_fragments_little_range1/$num_all_fragments;
my $num_rate_in_range=$num_fragments_in_range/$num_all_fragments;
my $ratio_value=$num_rate_in_range/$num_rate_little_range1;
printf RATE "AA\tthe reducing rate in range is\t%1.4f\n",$overall_rate_in_range;
printf RATE "BB\tthe ratio of fragments little than $range1 is\t%1.4f\n",$num_rate_little_range1;
printf RATE "CC\tthe ratio of fragments in range is\t%1.4f\n",$num_rate_in_range;
printf RATE "DD\tthe value of ratio to ratio is\t%.4f\n",$ratio_value;
print "$scafds_num chromsomes/scaffolds/contigs are successfully digested!\n";
