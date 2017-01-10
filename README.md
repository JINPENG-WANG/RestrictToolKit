# RestrictToolKit
For simulating the restriction digest of reference genome.
RestrictToolKit Manual

Version:	1.0.1
Author:	Jinpeng Wang, Li Li, Xuedi Du, Guofan Zhang.
Institue:	Institue of Oceanology, Chinese Academy of Sciences
Last update:	2014-10-11

1. Introduction

RestrictToolKit is used for the simulation of double-enzyme GBS approach. GBS(Genotyping by sequencing) is a popular 
method for reducing the complexity of the whole genome. As sequencing the whole genome is expensive and sometimes
unnecessary, sequencing the representative part of the whole genome becomes practicable and attractive. The most important
thing of representative sequencing is choosing the right enzymes to digest the whole genome. The appropriate enzymes should
satisfy several criteria including: the recognition sites of these two enzymes are evenly dispersed on each single chromsome;
 the recognition sites of these two enzymes are not or rarely located at the repetitive region of the genome; the fragments 
produced by two-enzyme digestion are neither too long nor too short, they should be in the suitable range which is most
 effective for PCR amplification; the complexity reducing rate shoule be in the suitable range(1%-10%). To order to obtain 
the right enzymes pair, we developed this software: RestrictToolKit. RestrictToolKit is designed by considering the four
standards described above as main factors when evaluating some enzymes pair. 

RestrictToolKit is easy to use and does not rely on any other softwares. However, if you want the result to be more precise,
the "-b" option is recommended for it will evoke the blastn program and assess the repetitiveness of the digested DNA fragments.

There is one important thing that should get your attention: RestrictToolKit only evaluates the particular DNA fragments  
whose length are located in the Range. Because in the practice, not all DNA fragments are amplified during the PCR process. In our 
 experiments, the fragments whose length are located between 90bp-750bp are mostly amplified while other fragments are slightly
 amplified. So we set the default range is 90bp-750bp. In other cases, you may want to choose particular fragments by cutting the gel,
 you should change the default range to your cutting range. About how to change the setting, please use the "-h"/"-H" option for help.

Specifically speaking, we evaluate five mandatory properties and one optional property of the fragments in range produced by input 
restrict-enzyme pair. The weights of these properties are different. RestrictToolKit adopts a 10 points scoring system. 
That means every property has a full score of 10 points and according to the weight of each property, all property scores are 
transformed into the final score of the two-enzyme pair which also has the full score of 10 points.

The first mandatory property is the reducing rate of the fragments in range, the calculation formula is "y=-(200/9)x+92/9", the weight is 20%.
In the formula, the "x" refers to the reducing rate, the "y" refers to the score of this property.
The second mandatory property is the distribution of all fragments, the calculation formula is "y=(4/3)x+16/3"; the weight is 30%. 
In the formula, the "x" refers to the distribution value which is defined as "a/b". The "b" refers to the ratio of the fragments with length little than 90bp among all fragments while the "a" refers to the ratio of the fragments with length between 90bp-750bp among all fragments. 
The "y" refers to the score of this property.
The third mandatory property is the ratio of the fragments in range, the calculation formula is "y=20x-4"; the weight is 20%. In the formula, 
the "x" refers to the ratio, the "y" refers to the score.
The forth mandatory property is the ratio of the fragments which can align to the CDS region, the calculation formula is "y=20x"; 
the weight is 10%. In the formula, the "x" refers to the ratio, the "y" refers to the score.
The fifth mandatory property is the ratio of the fragments which can align to the Intergenic region, the calculation formula is "y=-20x+14"; 
the weight is 10%. In the formula, the "x" refers to the ratio, the "y" refers to the score.
The optional property is the ratio of the fragments which are uniquely aligned to the reference genome, the calculation formula is "y=(40/3)x-2"; the weight is 10%. In the formula, the "x" refers to the ratio, the "y" refers to the score.
If the "-b" option is evoked, the optional property is scored according to the formula. Else, the score of optional property is assigned 
to be 10 points.

2. Installation

2.1 Prerequisites

RestrictToolKit should be run on any standard UNIX-like enrivonment(Apple OS X, Linux, etc). RestrictToolKit is an 
independent pipeline and can be run without any additional external software. To obtain more exact results, however, RestrictToolKit
calls the "blastn" program of "blastall" if you use the "-b" option. So if you able the "-b" option, please install the legacy BLAST
of NCBI and make sure the command "blastall" can be executed under the "RestrictToolKit" directory. It is strongly recommended that you add  
the directory containing the command "blastall" into the environment variable $PATH, so RestrictToolKit can evoke "blastall" correctly.

2.2 Computer requirements

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

2.3 Build the software

First, you need to unpack the package "RestrictToolKit.tar" where you want to put. 

% tar -xvf RestrictToolKit.tar

A directory named "RestrictToolKit" will be created and you can change the directory name before the configuration (We will talk about that 
later). Once you have finished the configuration, you should not change the directory name else you will need a reconfiguration!
Under this directory, there are four subdirectories named "gff", "reference", "run_results", "scripts" seperately and an executable file 
"configuration" and an introductory file "README". 
The "gff" subdirectory contains the GFF file(s) you put. The "reference" subdirectory contains the REFERENCE file(s) you put. 
The "run_results" subdirectory will contain results. Once you start the "RestrictSimulation" pipeline, a subdirectory named 
after the reference filename and the enzyme-pair will be created as "REFERENCE_ENZYME1_ENZYME2". 
Please replace the REFERENCE with the reference filename you putin and the ENZYME1 with the enzyme 1 name and the ENZYME2 with 
the enzyme 2 name.
Under this subdirectory, several result files will be created including the final result file and several mediating files. 
The "scripts" subdirectory contains the perl scripts used to run the pipeline. 
The "configuration" file can recognize the location of the software and create an executable file "RestrictSimulation". The "configuration" 
file is writen in Perl, so you can run it as "perl configuration" or "./configuration". You should only run it under the 
"RestrictToolKit"(or other name you change to) directory:

% cd RestrictToolKit
% ./configuration

The "RestrictSimulation" file will be created. It is the only command file you will use of this software.
It runs several perl scripts by order and can evoke the "blastall" program.  

Before you run the pipeline, please move your reference file(s) into the subdirectory "reference" and gff file(s) into "gff".

One more thing you need to do if you use the "-b" option is to format the reference file using the "formatdb" command of the legacy BLAST.
You can format your reference file under the subdirectory "reference" of directory "RestrictToolKit" or moving all files after formatting to
the subdirectory "reference" if you format the reference elsewhere.  

3. What types of data does RestrictToolKit support?

3.1 Files contents format

RestrictToolKit is designed to process the reference file in the "fasta" format. It should be in the following format:
>ChromsomeName1
ATCGATCGATCGATCGATCGATCGATCGATCG
>ChromsomeName2
ATCGATCGATCGATCGATCGATCGATCGATCG
The chromesome name should be unique. Any kinds of space chracters are not allowed, like blank space, tab etc.

In some cases, the reference file is in the following format:
>ChromsomeName1
ATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCG
>ChromsomeName2
ATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCG
.....
Sorry, this format is not allowed. Please transform this format into the first format described above.

For the "GFF" file, it should be in the standard GFF format, here is an example:
seqid		source	type	start	end	score	strand	phase	attributes			
scaffold39120	GLEAN	mRNA	36514	50111	.	+	.	ID=OYG_10003188_10008591;	
scaffold39120	GLEAN	CDS	36514	36534	.	+	.	Parent=OYG_10003188_10008591;	
scaffold39120	GLEAN	CDS	43276	43353	.	+	.	Parent=OYG_10003188_10008591;	
scaffold39120	GLEAN	CDS	43766	43868	.	+	.	Parent=OYG_10003188_10008591;	
scaffold39120	GLEAN	CDS	44710	44741	.	+	.	Parent=OYG_10003188_10008591;	
scaffold39120	GLEAN	CDS	49875	50111	.	+	.	Parent=OYG_10003188_10008591;	
scaffold838	GLEAN	mRNA	36580	65627	.	-	.	ID=OYG_10003281_10026064;	

3.2 Files format

Make sure all the reference files and the gff files are in the Unix\Mac OS X format. This refers to the CR/LF format. 
DOS/Windows use CR+LF to indicate a newline. However, Unix/Linux use LF to indicate a newline. So if you copy a file from Windows system 
to Linux system, make sure the file is saved as in Unix/Linux format. If any files are in DOS/winows format, the software will 
not run correctly.

4. Running the pipeline

To run the pipeline, you just need to input the command:
"perl RestrictSimulation -ref REFERENCE -gff GFF -e1 ENZYME1 -e2 ENZYME2 -b(optional)". 

You can run the command in the "RestrictToolKit" directory or you can also copy it into your "bin" directory which is added into the
environment variable $PATH, so you can run it anywhere you like.

5. Pipeline components

The pipeline consists of five perl scripts.
1) survey_restrict_locs.pl. It will read the reference sequences and output the locations of the two-enzyme recognization sites 
and the fragments sequnces in range.

2) ratio_fragments_mapped_to_CDS.pl. It will read the GFF file of the reference and the locations file produced by "survey_restrict_locs.pl" and then calculate the ratio of the fragments in range which can mapped to the CDS region.

3) ratio_fragments_mapped_to_Intergenic.pl. It will read the GFF file of the reference and the locations file produced by "survey_restrict_locs.pl" and then calculate the ratio of the fragments in range which can mapped to the Intergenic region. 

4) extract_info_of_blastn.pl. If the "-b" option is evoked, this script will extract useful information from the alignment file
produced by blastn.The blastn result file will be named like "blastn_alignments_of_REFERENCE_by_ENZYME1_and_ENZYME2.txt".

5) alignments_stats.pl. If the "-b" option is evoked, this script will calculate the ratio of fragments in range which are uniquely
mapped to the reference.

6. What do the fields mean in RestrictToolKit output files?

1) In the "characteristict_of_REFERENCE_by_ENZYME1_and_ENZYME2.txt" file, there are 8 lines and 3 rows. The first column contains FLAG values.
They are used to distinguish different lines. You can just leave them alone. The second column contains several properties of these two
enzymes on this reference. The third column contains the corresponding values. You can estimate the performance of these particular 
two enzymes on this particular reference based on these values. The higher of the HH value, the better of these two enzymes on this 
reference.

2) In the "restrict_locs_in_range_of_REFERENCE_by_ENZYME1_and_ENZYME2.txt" file, there are 6 lines. The descriptions are as follows:
>Name of fragments	ENZYME-a	location of ENZYME-a	ENZYME-b	location of ENZYME-b	length of fragments
ENZYME-a can be ENZYME1 or ENZYME2, so as to ENZYME-b.

3) In the "seq_of_frags_in_range_of_REFERENCE_by_ENZYME1_and_ENZYME2.txt" file, there are many double-line. The first line of the 
double_line is the name of the fragment in the ">chromsomeName-[0-9]+" format. The second line of the double_line is the sequence
of this fragment.

4) In the "blastn_alignments_of_REFERENCE_by_ENZYME1_and_ENZYME2.txt" file, there are 12 columns. The descriptions are as following:
Query Name;Subject Name;Identity;Length of alignment;Num. of Mismatch;Num. of Gap;Start of querry alignment;Stop of querry alignment;
Start of subject alignment;Stop of subject alignment;Expectation;Score of alignment. For more information, please refer 
to http://blast.ncbi.nlm.nih.gov/Blast.cgi .

5) In the "selected_blastn_alignments_of_REFERENCE_by_ENZYME1_and_ENZYME2.txt" file, there are 2 columns. The first column contains
the names of the fragments in range. The second column contains the alignment(s) number. RestrictToolKit calculate the unique alignment
ratio of the fragments according to this file.

7. Other explanations.

In the first version of this software, we evoke the legacy BLAST program. However, we will evoke the new blast BLAST+ in the next
version in the near future. 

8. Report bugs.

We will appreciate that if you report any problems and bugs when using this software. Please send email to "restricttoolkit@163.com".
RestrictToolKit Team: Jinpeng Wang, Li Li, Xuedi Du, Guofan Zhang.


