# antibiotic-resistance-process

The aim of this application is extract automatically the antibiotic resistance given a plate used in illumina sequencer machine.

execution:
python antibiotic_resistance.py --image images\plate.png

input: 
images/plate.png with a plate and ninety six samples

output:
images/outputXXX.png show the extracted samples, step bu step, given a plate and some requeriment like minimun rotation and scale

output/s1-A_0-0.png are the extracted resisance for each sample, the format is the next:
s: sample
next value is the row index
next value is the colmun index
next value is the absolute resistance founded
next value is the density of the resistance founded

ex: s4-A_122-0.23
is the sample 4-A with 122 pixels found as resistance with density of 23% for this well
