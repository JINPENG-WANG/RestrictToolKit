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
#if($filename=~/locs_of_frgts_by_(\w+)_(\w+).txt/){
#$frontenzyme=$1;
#$laterenzyme=$2;
#}
open GFF,"$path/gff/$gff" or die"cannot open $gff:$!";
open FGLCS,"$path/run_results/$result_dir/restrict_locs_in_range_of_${ref}_by_${frontenzyme}_and_${laterenzyme}.txt" or die"cannot open restrict_locs_in_range_of_${ref}_by_${frontenzyme}_and_${laterenzyme}.txt:$!";
open RSLT,">>$path/run_results/$result_dir/characteristics_of_${ref}_by_${frontenzyme}_and_${laterenzyme}.txt" or die"cannot write to characteristics_of_${ref}_by_${frontenzyme}_and_${laterenzyme}.txt:$!";

my($scfd_name_cds,$stt_pst_cds,$stp_pst_cds,$lgh_cds);
my($scfd_name_frgt,$stt_pst_frgt,$stp_pst_frgt,$lgh_frgt);
my(%cdslocs,%frgtlocs);
my $scfd_name_cds_tmp="";
my $scfd_name_frgt_tmp="";

#handle locs_of_cds.txt
while(<GFF>){
chomp;
my @line=split /\t/, $_;
if($line[2]=~/CDS/){
($scfd_name_cds,$stt_pst_cds,$stp_pst_cds)=@line[0,3,4];
$lgh_cds=$stp_pst_cds-$stt_pst_cds+1;
unless($scfd_name_cds eq $scfd_name_cds_tmp){
if(exists $cdslocs{$scfd_name_cds}){
push @{$cdslocs{$scfd_name_cds}},($stt_pst_cds,$stp_pst_cds,$lgh_cds);
}else{
$scfd_name_cds_tmp=$scfd_name_cds;
@{$cdslocs{$scfd_name_cds}}=($stt_pst_cds,$stp_pst_cds,$lgh_cds);
}
}
else{
push @{$cdslocs{$scfd_name_cds}},($stt_pst_cds,$stp_pst_cds,$lgh_cds);
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
$scfd_name_frgt_tmp=$scfd_name_frgt;
@{$frgtlocs{$scfd_name_frgt}}=($stt_pst_frgt,$stp_pst_frgt,$lgh_frgt);
}
else{
push @{$frgtlocs{$scfd_name_frgt}},($stt_pst_frgt,$stp_pst_frgt,$lgh_frgt);
}
}
my @scfds_frgt=sort keys %frgtlocs;
#compare these two hashes.report every scaffold in the cdslocs and every cds located in this scaffold. 
my $total_frgts=0;my $aln_frgts=0;
foreach my $scfd_frgt(@scfds_frgt){
my $scfd_cds=$scfd_frgt;
my @frgt_locs=@{$frgtlocs{$scfd_frgt}};
my $frgt_num=@frgt_locs/3;
my $stt_num_frgt=0;my $ix_stt_frgt=0;my $ix_stp_frgt=1;my $ix_lgh_frgt=2;
while($stt_num_frgt<$frgt_num){
$total_frgts++;
my $stt_frgt=$frgtlocs{$scfd_frgt}->[$ix_stt_frgt];
my $stp_frgt=$frgtlocs{$scfd_frgt}->[$ix_stp_frgt];
my $lghfrgt=$frgtlocs{$scfd_frgt}->[$ix_lgh_frgt];
my $alntms=0;
if(exists $cdslocs{$scfd_cds}){
my @cds_locs=@{$cdslocs{$scfd_cds}};
my $cds_num=@cds_locs/3;
my $stt_num_cds=0;my $ix_stt_cds=0;my $ix_stp_cds=1;my $ix_lgh_cds=2;
while($stt_num_cds<$cds_num){
my $stt_cds=$cdslocs{$scfd_cds}->[$ix_stt_cds];
my $stp_cds=$cdslocs{$scfd_cds}->[$ix_stp_cds];
my $lghcds=$cdslocs{$scfd_cds}->[$ix_lgh_cds];
if($lghfrgt<=$lghcds){
if(($stt_cds<=$stt_frgt && $stt_frgt<=$stp_cds) || ($stt_cds<=$stp_frgt && $stp_frgt<=$stp_cds)){
$alntms++;
$aln_frgts++;
}}
else{
if(($stt_frgt<=$stt_cds && $stt_cds<=$stp_frgt) || ($stt_frgt<=$stp_cds && $stp_cds<=$stp_frgt)){
$alntms++;
$aln_frgts++;
}
}
$stt_num_cds++;
$ix_stt_cds=$ix_stt_cds+3;
$ix_stp_cds=$ix_stp_cds+3;
$ix_lgh_cds=$ix_lgh_cds+3;
}
my $frgt_order=$stt_num_frgt+1;
#print RSLT "$scfd_frgt\t$frgt_order\t$alntms\n";
}else{
my $frgt_order=$stt_num_frgt+1;
#print RSLT"$scfd_frgt\t$frgt_order\t$alntms\n";
}
$stt_num_frgt++;
$ix_stt_frgt=$ix_stt_frgt+3;
$ix_stp_frgt=$ix_stp_frgt+3;
$ix_lgh_frgt=$ix_lgh_frgt+3;
}
}
my $rate_aln=$aln_frgts/$total_frgts;
print "there are $total_frgts frgts in total, there are $aln_frgts frgts which can aln cdss, the rate is $rate_aln\n"; 
printf RSLT "EE\tthe cds mapping ratio is\t%1.4f\n",$rate_aln;
