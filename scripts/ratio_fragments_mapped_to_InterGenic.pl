#!/usr/bin/perl -w
use strict;
my ($frontenzyme,$laterenzyme,$path,$ref,$gff,$result_dir);
while(@ARGV){
$_=shift @ARGV;
if($_=~/^-e1$/){$frontenzyme=shift @ARGV;}
elsif($_=~/^-e2$/){$laterenzyme=shift @ARGV;}
elsif($_=~/^-p$/){$path=shift @ARGV;}
elsif($_=~/^-gff$/){$gff=shift @ARGV;}
elsif($_=~/^-ref$/){$ref=shift @ARGV;}
elsif($_=~/^-dir$/){$result_dir=shift @ARGV;}
}
open GFF,"$path/gff/$gff" or die"cannot open $gff:$!";
open FGLCS,"$path/run_results/$result_dir/restrict_locs_in_range_of_${ref}_by_${frontenzyme}_and_${laterenzyme}.txt" or die"cannot open restrict_locs_in_range_of_${ref}_by_${frontenzyme}_and_${laterenzyme}.txt :$!";
open RSLT,">>$path/run_results/$result_dir/characteristics_of_${ref}_by_${frontenzyme}_and_${laterenzyme}.txt" or die"cannot write to characteristics_of_${ref}_by_${frontenzyme}_and_${laterenzyme}.txt:$!";
my($scfd_name_gene,$stt_pst_gene,$stp_pst_gene,$lgh_gene);
my($scfd_name_frgt,$stt_pst_frgt,$stp_pst_frgt,$lgh_frgt);
my(%genelocs,%frgtlocs);
my $scfd_name_gene_tmp="";
my $scfd_name_frgt_tmp="";
#handle locs_of_genes.txt
while(<GFF>){
chomp;
my @line=split /\t/,$_;
if($line[2]=~/mRNA/i){
($scfd_name_gene,$stt_pst_gene,$stp_pst_gene)=@line[0,3,4];
$lgh_gene=$stp_pst_gene-$stt_pst_gene+1;
unless($scfd_name_gene eq $scfd_name_gene_tmp){
if(exists $genelocs{$scfd_name_gene}){
push @{$genelocs{$scfd_name_gene}},($stt_pst_gene,$stp_pst_gene,$lgh_gene);
}else{
$scfd_name_gene_tmp=$scfd_name_gene;
@{$genelocs{$scfd_name_gene}}=($stt_pst_gene,$stp_pst_gene,$lgh_gene);
}
}
else{
push @{$genelocs{$scfd_name_gene}},($stt_pst_gene,$stp_pst_gene,$lgh_gene);
}
}
}
#handle locs_of_frgts.txt
while(<FGLCS>){
chomp;
my ($enzyme1,$enzyme2);
($scfd_name_frgt,$enzyme1,$stt_pst_frgt,$enzyme2,$stp_pst_frgt,$lgh_frgt)=split /\t/,$_;
$scfd_name_frgt=~s/>(\w+)-\d+/$1/;
unless($scfd_name_frgt eq $scfd_name_frgt_tmp){
if(exists $frgtlocs{$scfd_name_frgt}){
push @{$frgtlocs{$scfd_name_frgt}},($stt_pst_frgt,$stp_pst_frgt,$lgh_frgt);
}else{
$scfd_name_frgt_tmp=$scfd_name_frgt;
@{$frgtlocs{$scfd_name_frgt}}=($stt_pst_frgt,$stp_pst_frgt,$lgh_frgt);
}
}
else{
push @{$frgtlocs{$scfd_name_frgt}},($stt_pst_frgt,$stp_pst_frgt,$lgh_frgt);
}
}
my @scfds_frgt=sort keys %frgtlocs;
my $total_frgts=0;my $aln_frgts=0;
foreach my $scfd_frgt(@scfds_frgt){
my $scfd_gene=$scfd_frgt;
my @frgt_locs=@{$frgtlocs{$scfd_frgt}};
my $frgt_num=@frgt_locs/3;
my $stt_num_frgt=0;my $ix_stt_frgt=0;my $ix_stp_frgt=1;my $ix_lgh_frgt=2;
while($stt_num_frgt<$frgt_num){
$total_frgts++;
my $stt_frgt=$frgtlocs{$scfd_frgt}->[$ix_stt_frgt];
my $stp_frgt=$frgtlocs{$scfd_frgt}->[$ix_stp_frgt];
my $lghfrgt=$frgtlocs{$scfd_frgt}->[$ix_lgh_frgt];
if(exists $genelocs{$scfd_gene}){
my @gene_locs=@{$genelocs{$scfd_gene}};
my $gene_num=@gene_locs/3;
my $stt_num_gene=0;my $ix_stt_gene=0;my $ix_stp_gene=1;my $ix_lgh_gene=2;
my $judge1=0;
while($stt_num_gene<$gene_num){
my $stt_gene=$genelocs{$scfd_gene}->[$ix_stt_gene];
my $stp_gene=$genelocs{$scfd_gene}->[$ix_stp_gene];
my $lghgene=$genelocs{$scfd_gene}->[$ix_lgh_gene];
if($stt_frgt>=$stp_gene || $stp_frgt<=$stt_gene){
my $judge2=0;
$judge1=$judge1+$judge2;
}else{
my $judge2=1;
$judge1=$judge1+$judge2;
}
$stt_num_gene++;
$ix_stt_gene=$ix_stt_gene+3;
$ix_stp_gene=$ix_stp_gene+3;
$ix_lgh_gene=$ix_lgh_gene+3;
}
if($judge1==0){
$aln_frgts++;
}
}
else{
$aln_frgts++;
}
$stt_num_frgt++;
$ix_stt_frgt=$ix_stt_frgt+3;
$ix_stp_frgt=$ix_stp_frgt+3;
$ix_lgh_frgt=$ix_lgh_frgt+3;
}
}
my $rate_aln=$aln_frgts/$total_frgts;
print "there are $total_frgts frgts in total, there are $aln_frgts frgts which can aln intergenic regions, the rate is $rate_aln\n"; 
printf RSLT "FF\tthe intergenic mapping ratio is\t%1.4f\n",$rate_aln;
