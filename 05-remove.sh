#!/bin/sh

#program order: trim_galore, bbmerge, clumpify, seal, bbduk#
#This is the fifth step: remove sequences that don't have the primer regions 


InDir="../04seal" # Separated three fragments libraries based on target gene sequences
OutDir="../05remove" # Folder to keep sequences that have both primer aligned regions 


samplelistb="1_S1_R 2_S2_R 3_S3_R 4_S4_R 5_S5_R 6_S6_R 7_S7_R 8_S8_R 9_S9_R 10_S10_R 11_S11_R 12_S12_R 13_S13_R 14_S14_R 15_S15_R 16_S16_R 17_S17_R 18_S18_R 19_S19_R 20_S20_R 21_S21_R 22_S22_R 23_S23_R 24_S24_R 65_S65_R 66_S66_R 69_S69_R 70_S70_R 71_S71_R 72_S72_R"

samplelistc="25_S25_R 26_S26_R 27_S27_R 28_S28_R"

samplelistd="29_S29_R 30_S30_R 31_S31_R"

samplelistj="33_S33_R 34_S34_R 35_S35_R 36_S36_R"

samplelistk12="37_S37_R 38_S38_R"

samplelistk34="39_S39_R 40_S40_R"

samplelistl="41_S41_R 42_S42_R 43_S43_R 44_S44_R 45_S45_R 46_S46_R 47_S47_R 48_S48_R 49_S49_R 50_S50_R 51_S51_R 52_S52_R 53_S53_R 54_S54_R 55_S55_R 56_S56_R 57_S57_R 58_S58_R 59_S59_R 60_S60_R 67_S67_R 68_S68_R"

samplelistm="61_S61_R 62_S62_R 63_S63_R 64_S64_R"

A1F="GGTGGATCTGGTggcAGC"  #18
B1R="GCAGCTGGACCTGCCTTAA"  #19
B2F="ACGCACGAGTCCTCTTCAG"  #19
B2R="GAGTTCGGTGGGGTCATGT"  #19
B3F="TCAGGGACGGGGTGAACTG"  #19
A3R="caacatgctccctcaatcgg"  #20
  #
C1R="CAGCAGGTCCAGCCCTTTC"  #19
C2F="GTGGCACGTACATCTCCGT"  #19
C2R="GTTCGGCGGTGCATTGGC"  #18
C3F="TGTTAACTGGGGGCGCATC"  #19
  #
  #
D1R="GCTGGACCTGCCTTGTCTC"  #19
D2F="GCACGTACGTCGCCATTAC"  #19
D2R="CTTTGGCGGCGCTTTGTGT"  #19
D3F="CGGTGTTAATTGGGGGCGT"  #19
  #
  #
J1R="GCTATCACTTAGGGGGTGCC"  #20
J2F="CATGAGTTGGCGCGTGAC"  #18
J2R="TGGGGCCGCATCGTATCC"  #18
J3F="GGCCAACGAAATGTTCTCGG"  #20
  #
  #
K121R="CTGCGTGGTGCGAAAGATAC"  #20
K122F="CGAACAGGCAAAGAATTTGGTGA"  #23
K122R="TGGGGACGTATTGTCGCTC"  #19
K123F="AACTGTTTCGTGACGGCGT"  #19
  #
  #
K341R="AGCTGGCCCAGCATTATCAC"  #20
K342F="GCACGTACCTCTCCGCTTC"  #19
K342R="TGCTTTCGGTGGGGCCCT"  #18
K343F="CCGCGACGGTATCAATTGG"  #19
  #
  #
L1R="GTACCTTCGGGAGCAGGC"  #18
L2F="TACCGGCAGTCGCTGGAG"  #18
L2R="TGGGGCAGGATTGTGACTC"  #19
L3F="ATGTTTTCAGCGACGGCGT"  #19
  #
  #
M1R="CAGGACCGGCACTTTCCC"  #18
M2F="CTGCAAACACCCGCAGCA"  #18
M2R="TCGGTGGAGCTCTTGCAGT"  #19
M3F="AATTGGGGGCGTATCGTCG"  #19


#with mm=f, kmer based sorting is more strict may better for the next step


for s in $samplelistb

do 
#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1.fastq.gz outm=$OutDir/$s.f1f.fastq.gz out=$OutDir/$s.f1nof.fastq.gz stats=$OutDir/$s.f1f.stats.txt k=18 literal=$A1F mm=f 
../bbduk.sh in=$OutDir/$s.f1f.fastq.gz outm=$OutDir/$s.f1fr.fastq.gz out=$OutDir/$s.f1nofr.fastq.gz stats=$OutDir/$s.f1fr.stats.txt k=19 literal=$B1R mm=f

#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2.fastq.gz outm=$OutDir/$s.f2f.fastq.gz out=$OutDir/$s.f2nof.fastq.gz stats=$OutDir/$s.f2f.stats.txt k=19 literal=$B2F mm=f
../bbduk.sh in=$OutDir/$s.f2f.fastq.gz outm=$OutDir/$s.f2fr.fastq.gz out=$OutDir/$s.f2nofr.fastq.gz stats=$OutDir/$s.f2fr.stats.txt k=19 literal=$B2R mm=f


#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3.fastq.gz outm=$OutDir/$s.f3f.fastq.gz out=$OutDir/$s.f3nof.fastq.gz stats=$OutDir/$s.f3f.stats.txt k=19 literal=$B3F mm=f
../bbduk.sh in=$OutDir/$s.f3f.fastq.gz outm=$OutDir/$s.f3fr.fastq.gz out=$OutDir/$s.f3nofr.fastq.gz stats=$OutDir/$s.f3fr.stats.txt k=20 literal=$A3R mm=f

done

for s in $samplelistc

do 
#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1.fastq.gz outm=$OutDir/$s.f1f.fastq.gz out=$OutDir/$s.f1nof.fastq.gz stats=$OutDir/$s.f1f.stats.txt k=18 literal=$A1F mm=f 
../bbduk.sh in=$OutDir/$s.f1f.fastq.gz outm=$OutDir/$s.f1fr.fastq.gz out=$OutDir/$s.f1nofr.fastq.gz stats=$OutDir/$s.f1fr.stats.txt k=19 literal=$C1R mm=f

#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2.fastq.gz outm=$OutDir/$s.f2f.fastq.gz out=$OutDir/$s.f2nof.fastq.gz stats=$OutDir/$s.f2f.stats.txt k=19 literal=$C2F mm=f
../bbduk.sh in=$OutDir/$s.f2f.fastq.gz outm=$OutDir/$s.f2fr.fastq.gz out=$OutDir/$s.f2nofr.fastq.gz stats=$OutDir/$s.f2fr.stats.txt k=18 literal=$C2R mm=f


#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3.fastq.gz outm=$OutDir/$s.f3f.fastq.gz out=$OutDir/$s.f3nof.fastq.gz stats=$OutDir/$s.f3f.stats.txt k=19 literal=$C3F mm=f
../bbduk.sh in=$OutDir/$s.f3f.fastq.gz outm=$OutDir/$s.f3fr.fastq.gz out=$OutDir/$s.f3nofr.fastq.gz stats=$OutDir/$s.f3fr.stats.txt k=20 literal=$A3R mm=f

done

for s in $samplelistd

do 
#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1.fastq.gz outm=$OutDir/$s.f1f.fastq.gz out=$OutDir/$s.f1nof.fastq.gz stats=$OutDir/$s.f1f.stats.txt k=18 literal=$A1F mm=f 
../bbduk.sh in=$OutDir/$s.f1f.fastq.gz outm=$OutDir/$s.f1fr.fastq.gz out=$OutDir/$s.f1nofr.fastq.gz stats=$OutDir/$s.f1fr.stats.txt k=19 literal=$D1R mm=f

#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2.fastq.gz outm=$OutDir/$s.f2f.fastq.gz out=$OutDir/$s.f2nof.fastq.gz stats=$OutDir/$s.f2f.stats.txt k=19 literal=$D2F mm=f
../bbduk.sh in=$OutDir/$s.f2f.fastq.gz outm=$OutDir/$s.f2fr.fastq.gz out=$OutDir/$s.f2nofr.fastq.gz stats=$OutDir/$s.f2fr.stats.txt k=19 literal=$D2R mm=f


#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3.fastq.gz outm=$OutDir/$s.f3f.fastq.gz out=$OutDir/$s.f3nof.fastq.gz stats=$OutDir/$s.f3f.stats.txt k=19 literal=$D3F mm=f
../bbduk.sh in=$OutDir/$s.f3f.fastq.gz outm=$OutDir/$s.f3fr.fastq.gz out=$OutDir/$s.f3nofr.fastq.gz stats=$OutDir/$s.f3fr.stats.txt k=20 literal=$A3R mm=f


done

for s in $samplelistj

do 

#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1.fastq.gz outm=$OutDir/$s.f1f.fastq.gz out=$OutDir/$s.f1nof.fastq.gz stats=$OutDir/$s.f1f.stats.txt k=18 literal=$A1F mm=f 
../bbduk.sh in=$OutDir/$s.f1f.fastq.gz outm=$OutDir/$s.f1fr.fastq.gz out=$OutDir/$s.f1nofr.fastq.gz stats=$OutDir/$s.f1fr.stats.txt k=20 literal=$J1R mm=f

#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2.fastq.gz outm=$OutDir/$s.f2f.fastq.gz out=$OutDir/$s.f2nof.fastq.gz stats=$OutDir/$s.f2f.stats.txt k=18 literal=$J2F mm=f
../bbduk.sh in=$OutDir/$s.f2f.fastq.gz outm=$OutDir/$s.f2fr.fastq.gz out=$OutDir/$s.f2nofr.fastq.gz stats=$OutDir/$s.f2fr.stats.txt k=18 literal=$J2R mm=f


#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3.fastq.gz outm=$OutDir/$s.f3f.fastq.gz out=$OutDir/$s.f3nof.fastq.gz stats=$OutDir/$s.f3f.stats.txt k=20 literal=$J3F mm=f
../bbduk.sh in=$OutDir/$s.f3f.fastq.gz outm=$OutDir/$s.f3fr.fastq.gz out=$OutDir/$s.f3nofr.fastq.gz stats=$OutDir/$s.f3fr.stats.txt k=20 literal=$A3R mm=f

done

for s in $samplelistk12

do 
#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1.fastq.gz outm=$OutDir/$s.f1f.fastq.gz out=$OutDir/$s.f1nof.fastq.gz stats=$OutDir/$s.f1f.stats.txt k=18 literal=$A1F mm=f 
../bbduk.sh in=$OutDir/$s.f1f.fastq.gz outm=$OutDir/$s.f1fr.fastq.gz out=$OutDir/$s.f1nofr.fastq.gz stats=$OutDir/$s.f1fr.stats.txt k=20 literal=$K121R mm=f

#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2.fastq.gz outm=$OutDir/$s.f2f.fastq.gz out=$OutDir/$s.f2nof.fastq.gz stats=$OutDir/$s.f2f.stats.txt k=23 literal=$K122F mm=f
../bbduk.sh in=$OutDir/$s.f2f.fastq.gz outm=$OutDir/$s.f2fr.fastq.gz out=$OutDir/$s.f2nofr.fastq.gz stats=$OutDir/$s.f2fr.stats.txt k=19 literal=$K122R mm=f


#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3.fastq.gz outm=$OutDir/$s.f3f.fastq.gz out=$OutDir/$s.f3nof.fastq.gz stats=$OutDir/$s.f3f.stats.txt k=19 literal=$K123F mm=f
../bbduk.sh in=$OutDir/$s.f3f.fastq.gz outm=$OutDir/$s.f3fr.fastq.gz out=$OutDir/$s.f3nofr.fastq.gz stats=$OutDir/$s.f3fr.stats.txt k=20 literal=$A3R mm=f

done

for s in $samplelistk34

do 
#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1.fastq.gz outm=$OutDir/$s.f1f.fastq.gz out=$OutDir/$s.f1nof.fastq.gz stats=$OutDir/$s.f1f.stats.txt k=18 literal=$A1F mm=f 
../bbduk.sh in=$OutDir/$s.f1f.fastq.gz outm=$OutDir/$s.f1fr.fastq.gz out=$OutDir/$s.f1nofr.fastq.gz stats=$OutDir/$s.f1fr.stats.txt k=20 literal=$K341R mm=f

#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2.fastq.gz outm=$OutDir/$s.f2f.fastq.gz out=$OutDir/$s.f2nof.fastq.gz stats=$OutDir/$s.f2f.stats.txt k=19 literal=$K342F mm=f
../bbduk.sh in=$OutDir/$s.f2f.fastq.gz outm=$OutDir/$s.f2fr.fastq.gz out=$OutDir/$s.f2nofr.fastq.gz stats=$OutDir/$s.f2fr.stats.txt k=18 literal=$K342R mm=f


#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3.fastq.gz outm=$OutDir/$s.f3f.fastq.gz out=$OutDir/$s.f3nof.fastq.gz stats=$OutDir/$s.f3f.stats.txt k=19 literal=$K343F mm=f
../bbduk.sh in=$OutDir/$s.f3f.fastq.gz outm=$OutDir/$s.f3fr.fastq.gz out=$OutDir/$s.f3nofr.fastq.gz stats=$OutDir/$s.f3fr.stats.txt k=20 literal=$A3R mm=f

done

for s in $samplelistl

do 
#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1.fastq.gz outm=$OutDir/$s.f1f.fastq.gz out=$OutDir/$s.f1nof.fastq.gz stats=$OutDir/$s.f1f.stats.txt k=18 literal=$A1F mm=f
../bbduk.sh in=$OutDir/$s.f1f.fastq.gz outm=$OutDir/$s.f1fr.fastq.gz out=$OutDir/$s.f1nofr.fastq.gz stats=$OutDir/$s.f1fr.stats.txt k=18 literal=$L1R mm=f

#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2.fastq.gz outm=$OutDir/$s.f2f.fastq.gz out=$OutDir/$s.f2nof.fastq.gz stats=$OutDir/$s.f2f.stats.txt k=18 literal=$L2F mm=f
../bbduk.sh in=$OutDir/$s.f2f.fastq.gz outm=$OutDir/$s.f2fr.fastq.gz out=$OutDir/$s.f2nofr.fastq.gz stats=$OutDir/$s.f2fr.stats.txt k=19 literal=$L2R mm=f


#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3.fastq.gz outm=$OutDir/$s.f3f.fastq.gz out=$OutDir/$s.f3nof.fastq.gz stats=$OutDir/$s.f3f.stats.txt k=19 literal=$L3F mm=f
../bbduk.sh in=$OutDir/$s.f3f.fastq.gz outm=$OutDir/$s.f3fr.fastq.gz out=$OutDir/$s.f3nofr.fastq.gz stats=$OutDir/$s.f3fr.stats.txt k=20 literal=$A3R mm=f


done



for s in $samplelistm

do 
#This is used to trim fragment1 primers#
../bbduk.sh in=$InDir/$s.f1.fastq.gz outm=$OutDir/$s.f1f.fastq.gz out=$OutDir/$s.f1nof.fastq.gz stats=$OutDir/$s.f1f.stats.txt k=18 literal=$A1F mm=f
../bbduk.sh in=$OutDir/$s.f1f.fastq.gz outm=$OutDir/$s.f1fr.fastq.gz out=$OutDir/$s.f1nofr.fastq.gz stats=$OutDir/$s.f1fr.stats.txt k=18 literal=$M1R mm=f

#This is used to trim fragment2 primers#
../bbduk.sh in=$InDir/$s.f2.fastq.gz outm=$OutDir/$s.f2f.fastq.gz out=$OutDir/$s.f2nof.fastq.gz stats=$OutDir/$s.f2f.stats.txt k=18 literal=$M2F mm=f
../bbduk.sh in=$OutDir/$s.f2f.fastq.gz outm=$OutDir/$s.f2fr.fastq.gz out=$OutDir/$s.f2nofr.fastq.gz stats=$OutDir/$s.f2fr.stats.txt k=19 literal=$M2R


#This is used to trim fragment3 primers#
../bbduk.sh in=$InDir/$s.f3.fastq.gz outm=$OutDir/$s.f3f.fastq.gz out=$OutDir/$s.f3nof.fastq.gz stats=$OutDir/$s.f3f.stats.txt k=19 literal=$M3F mm=f
../bbduk.sh in=$OutDir/$s.f3f.fastq.gz outm=$OutDir/$s.f3fr.fastq.gz out=$OutDir/$s.f3nofr.fastq.gz stats=$OutDir/$s.f3fr.stats.txt k=20 literal=$A3R mm=f


done




