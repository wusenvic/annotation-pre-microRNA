This is a new method for annotation of pre-microRNAs from total RNA-seq rawdata.
pre-requirement of software and files:
1. vsearch
2. pre-microRNA reference sequences in FASTA file
3. mature microRNA reference sequences in FASTA file.
   
FIrstly, preprocess rowdata with quality control for rawdata.
Secondely, go through "main_concept" file and perform sequence filter following the instruction.
Thirdly, using "anotation_pre_microRNA_based_on_kmer.py" to do a final annotation.
Lastly, combined single sample files to a matrix via the script of "combined*.py"  
