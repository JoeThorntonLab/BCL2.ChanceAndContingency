#!/bin/sh

#program order: trim_galore, bbmerge, clumpify, seal, bbduk#
#This is the first step: trimgalore, a quality based trim of sequence#

#users need to download Trim Galore from Babraham Bioinformatic: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/#

InDir="../fastq"  #Raw sequence data from MiSeq run. Users should change to their directory when use, same for all other folder and script paths.
OutDir="../01trim"  #Trimmed sequence data

samplelist="1_S1_R 2_S2_R 3_S3_R 4_S4_R 5_S5_R 6_S6_R 7_S7_R 8_S8_R 9_S9_R 10_S10_R 11_S11_R 12_S12_R 13_S13_R 14_S14_R 15_S15_R 16_S16_R 17_S17_R 18_S18_R 19_S19_R 20_S20_R 21_S21_R 22_S22_R 23_S23_R 24_S24_R 25_S25_R 26_S26_R 27_S27_R 28_S28_R 29_S29_R 30_S30_R 31_S31_R 33_S33_R 34_S34_R 35_S35_R 36_S36_R 37_S37_R 38_S38_R 39_S39_R 40_S40_R 41_S41_R 42_S42_R 43_S43_R 44_S44_R 45_S45_R 46_S46_R 47_S47_R 48_S48_R 49_S49_R 50_S50_R 51_S51_R 52_S52_R 53_S53_R 54_S54_R 55_S55_R 56_S56_R 57_S57_R 58_S58_R 59_S59_R 60_S60_R 61_S61_R 62_S62_R 63_S63_R 64_S64_R 65_S65_R 66_S66_R 67_S67_R 68_S68_R 69_S69_R 70_S70_R 71_S71_R 72_S72_R"

for s in $samplelist

do 

../trim_galore -paired --dont_gzip $InDir/$s.1.fastq.gz $InDir/$s.2.fastq.gz -o $OutDir 

done
