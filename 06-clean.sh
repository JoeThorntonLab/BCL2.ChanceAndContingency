#!/bin/sh

#program order: trim_galore, bbmerge, clumpify, seal, bbduk#
#This is the sixth step: remove the primer binding region of each library sequences


InDir="../05remove" # Sequence data that have the primer regions to be cleaned
MidDir="../06clean/mid" #Intermediate data during this process
OutDir="../06clean" #Folder to keep final cleaned sequence data that removed all primer binding regions


samplelistb="1_S1_R 2_S2_R 3_S3_R 4_S4_R 5_S5_R 6_S6_R 7_S7_R 8_S8_R 9_S9_R 10_S10_R 11_S11_R 12_S12_R 13_S13_R 14_S14_R 15_S15_R 16_S16_R 17_S17_R 18_S18_R 19_S19_R 20_S20_R 21_S21_R 22_S22_R 23_S23_R 24_S24_R 65_S65_R 66_S66_R 69_S69_R 70_S70_R 71_S71_R 72_S72_R"

samplelistc="25_S25_R 26_S26_R 27_S27_R 28_S28_R"

samplelistd="29_S29_R 30_S30_R 31_S31_R"

samplelistj="33_S33_R 34_S34_R 35_S35_R 36_S36_R"

samplelistk12="37_S37_R 38_S38_R"

samplelistk34="39_S39_R 40_S40_R"

samplelistl="41_S41_R 42_S42_R 43_S43_R 44_S44_R 45_S45_R 46_S46_R 47_S47_R 48_S48_R 49_S49_R 50_S50_R 51_S51_R 52_S52_R 53_S53_R 54_S54_R 55_S55_R 56_S56_R 57_S57_R 58_S58_R 59_S59_R 60_S60_R 67_S67_R 68_S68_R"

samplelistm="61_S61_R 62_S62_R 63_S63_R 64_S64_R"

A1F="GGTGGATCTGGTggcAGC"
B1R="GCAGCTGGACCTGCCTTAA"
B2F="ACGCACGAGTCCTCTTCAG"
B2R="GAGTTCGGTGGGGTCATGT"
B3F="TCAGGGACGGGGTGAACTG"
A3R="caacatgctccctcaatcgg"

C1R="CAGCAGGTCCAGCCCTTTC"
C2F="GTGGCACGTACATCTCCGT"
C2R="GTTCGGCGGTGCATTGGC"
C3F="TGTTAACTGGGGGCGCATC"

D1R="GCTGGACCTGCCTTGTCTC"
D2F="GCACGTACGTCGCCATTAC"
D2R="CTTTGGCGGCGCTTTGTGT"
D3F="CGGTGTTAATTGGGGGCGT"

J1R="GCTATCACTTAGGGGGTGCC"
J2F="CATGAGTTGGCGCGTGAC"
J2R="TGGGGCCGCATCGTATCC"
J3F="GGCCAACGAAATGTTCTCGG"

K121R="CTGCGTGGTGCGAAAGATAC"
K122F="CGAACAGGCAAAGAATTTGGTGA"
K122R="TGGGGACGTATTGTCGCTC"
K123F="AACTGTTTCGTGACGGCGT"

K341R="AGCTGGCCCAGCATTATCAC"
K342F="GCACGTACCTCTCCGCTTC"
K342R="TGCTTTCGGTGGGGCCCT"
K343F="CCGCGACGGTATCAATTGG"

L1R="GTACCTTCGGGAGCAGGC"
L2F="TACCGGCAGTCGCTGGAG"
L2R="TGGGGCAGGATTGTGACTC"
L3F="ATGTTTTCAGCGACGGCGT"

M1R="CAGGACCGGCACTTTCCC"
M2F="CTGCAAACACCCGCAGCA"
M2R="TCGGTGGAGCTCTTGCAGT"
M3F="AATTGGGGGCGTATCGTCG"


for s in $samplelistb

do 
#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1fr.fastq.gz out=$MidDir/$s.f1-l.fastq.gz k=18 literal=$A1F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f1-l.fastq.gz out=$OutDir/$s.F01.fastq k=19 literal=$B1R ktrim=r mm=f
../bbduk.sh in=$OutDir/$s.F01.fastq lhist=$OutDir/$s.F01.lhist.txt  #This is used to count the fragment length distributions in the library#


#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2fr.fastq.gz out=$MidDir/$s.f2-l.fastq.gz k=19 literal=$B2F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f2-l.fastq.gz out=$OutDir/$s.F02.fastq k=19 literal=$B2R ktrim=r mm=f
../bbduk.sh in=$OutDir/$s.F02.fastq lhist=$OutDir/$s.F02.lhist.txt 

#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3fr.fastq.gz out=$MidDir/$s.f3-l.fastq.gz k=19 literal=$B3F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f3-l.fastq.gz out=$OutDir/$s.F03.fastq k=20 literal=$A3R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F03.fastq lhist=$OutDir/$s.F03.lhist.txt 

done

for s in $samplelistc

do 
#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1fr.fastq.gz out=$MidDir/$s.f1-l.fastq.gz k=18 literal=$A1F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f1-l.fastq.gz out=$OutDir/$s.F01.fastq k=19 literal=$C1R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F01.fastq lhist=$OutDir/$s.F01.lhist.txt 

#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2fr.fastq.gz out=$MidDir/$s.f2-l.fastq.gz k=19 literal=$C2F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f2-l.fastq.gz out=$OutDir/$s.F02.fastq k=18 literal=$C2R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F02.fastq lhist=$OutDir/$s.F02.lhist.txt 


#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3fr.fastq.gz out=$MidDir/$s.f3-l.fastq.gz k=19 literal=$C3F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f3-l.fastq.gz out=$OutDir/$s.F03.fastq k=20 literal=$A3R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F03.fastq lhist=$OutDir/$s.F03.lhist.txt 


done

for s in $samplelistd

do 
#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1fr.fastq.gz out=$MidDir/$s.f1-l.fastq.gz k=18 literal=$A1F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f1-l.fastq.gz out=$OutDir/$s.F01.fastq k=19 literal=$D1R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F01.fastq lhist=$OutDir/$s.F01.lhist.txt 

#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2fr.fastq.gz out=$MidDir/$s.f2-l.fastq.gz k=19 literal=$D2F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f2-l.fastq.gz out=$OutDir/$s.F02.fastq k=19 literal=$D2R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F02.fastq lhist=$OutDir/$s.F02.lhist.txt 


#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3fr.fastq.gz out=$MidDir/$s.f3-l.fastq.gz k=19 literal=$D3F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f3-l.fastq.gz out=$OutDir/$s.F03.fastq k=20 literal=$A3R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F03.fastq lhist=$OutDir/$s.F03.lhist.txt 


done

for s in $samplelistj

do 
#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1fr.fastq.gz out=$MidDir/$s.f1-l.fastq.gz k=18 literal=$A1F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f1-l.fastq.gz out=$OutDir/$s.F01.fastq k=20 literal=$J1R ktrim=r mm=f
../bbduk.sh in=$OutDir/$s.F01.fastq lhist=$OutDir/$s.F01.lhist.txt 

#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2fr.fastq.gz out=$MidDir/$s.f2-l.fastq.gz k=18 literal=$J2F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f2-l.fastq.gz out=$OutDir/$s.F02.fastq k=18 literal=$J2R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F02.fastq lhist=$OutDir/$s.F02.lhist.txt 


#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3fr.fastq.gz out=$MidDir/$s.f3-l.fastq.gz k=20 literal=$J3F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f3-l.fastq.gz out=$OutDir/$s.F03.fastq k=20 literal=$A3R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F03.fastq lhist=$OutDir/$s.F03.lhist.txt 

done

for s in $samplelistk12

do 
#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1fr.fastq.gz out=$MidDir/$s.f1-l.fastq.gz k=18 literal=$A1F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f1-l.fastq.gz out=$OutDir/$s.F01.fastq k=20 literal=$K121R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F01.fastq lhist=$OutDir/$s.F01.lhist.txt 

#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2fr.fastq.gz out=$MidDir/$s.f2-l.fastq.gz k=23 literal=$K122F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f2-l.fastq.gz out=$OutDir/$s.F02.fastq k=19 literal=$K122R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F02.fastq lhist=$OutDir/$s.F02.lhist.txt 


#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3fr.fastq.gz out=$MidDir/$s.f3-l.fastq.gz k=19 literal=$K123F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f3-l.fastq.gz out=$OutDir/$s.F03.fastq k=20 literal=$A3R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F03.fastq lhist=$OutDir/$s.F03.lhist.txt 

done

for s in $samplelistk34

do 
#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1fr.fastq.gz out=$MidDir/$s.f1-l.fastq.gz k=18 literal=$A1F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f1-l.fastq.gz out=$OutDir/$s.F01.fastq k=20 literal=$K341R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F01.fastq lhist=$OutDir/$s.F01.lhist.txt 

#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2fr.fastq.gz out=$MidDir/$s.f2-l.fastq.gz k=19 literal=$K342F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f2-l.fastq.gz out=$OutDir/$s.F02.fastq k=18 literal=$K342R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F02.fastq lhist=$OutDir/$s.F02.lhist.txt 


#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3fr.fastq.gz out=$MidDir/$s.f3-l.fastq.gz k=19 literal=$K343F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f3-l.fastq.gz out=$OutDir/$s.F03.fastq k=20 literal=$A3R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F03.fastq lhist=$OutDir/$s.F03.lhist.txt 

done

for s in $samplelistl

do 

#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1fr.fastq.gz out=$MidDir/$s.f1-l.fastq.gz k=18 literal=$A1F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f1-l.fastq.gz out=$OutDir/$s.F01.fastq k=18 literal=$L1R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F01.fastq lhist=$OutDir/$s.F01.lhist.txt 

#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2fr.fastq.gz out=$MidDir/$s.f2-l.fastq.gz k=18 literal=$L2F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f2-l.fastq.gz out=$OutDir/$s.F02.fastq k=19 literal=$L2R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F02.fastq lhist=$OutDir/$s.F02.lhist.txt 


#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3fr.fastq.gz out=$MidDir/$s.f3-l.fastq.gz k=19 literal=$L3F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f3-l.fastq.gz out=$OutDir/$s.F03.fastq k=20 literal=$A3R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F03.fastq lhist=$OutDir/$s.F03.lhist.txt 

done



for s in $samplelistm

do 
#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1fr.fastq.gz out=$MidDir/$s.f1-l.fastq.gz k=18 literal=$A1F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f1-l.fastq.gz out=$OutDir/$s.F01.fastq k=18 literal=$M1R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F01.fastq lhist=$OutDir/$s.F01.lhist.txt 

#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2fr.fastq.gz out=$MidDir/$s.f2-l.fastq.gz k=18 literal=$M2F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f2-l.fastq.gz out=$OutDir/$s.F02.fastq k=19 literal=$M2R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F02.fastq lhist=$OutDir/$s.F02.lhist.txt 


#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3fr.fastq.gz out=$MidDir/$s.f3-l.fastq.gz k=19 literal=$M3F ktrim=l mm=f
../bbduk.sh in=$MidDir/$s.f3-l.fastq.gz out=$OutDir/$s.F03.fastq k=20 literal=$A3R ktrim=r mm=f 
../bbduk.sh in=$OutDir/$s.F03.fastq lhist=$OutDir/$s.F03.lhist.txt 

done




