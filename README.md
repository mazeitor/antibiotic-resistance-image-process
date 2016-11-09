# antibiotic-resistance-process

The aim of this application is extract automatically the antibiotic resistance given a 96-well plate used in sequencer machine.

##execution:
python antibiotic_resistance.py --image images/plate.png

###input:
images/plate.png with a plate and ninety six samples

###output:
images/<platename>/putputXXX.png image with extracted samples
images/<platename>/<row>_<column>_<resistance>_<density>.png cropped image of extracted sample
images/<platename>/report.json json with extracted samples for each well
images/<platename>/log.txt log 

output/s1-A_0-0.png are the extracted resisance for each sample, the format is the next:
s: sample
next value is the row index
next value is the colmun index
next value is the absolute resistance founded
next value is the density of the resistance founded

ex: s4-A_122-0.23
is the sample 4-A with 122 pixels found as resistance with density of 23% for this well

##installing dependencies
###opencv
sudo apt-get install build-essential
sudo apt-get install cmake git libgtk2.0-dev pkg-config libavcodec-dev libavformat-dev libswscale-dev
sudo apt-get install python-opencv

###scilab
sudo apt-get install python-scipy

###python-tk
sudo apt-get install python-tk

###pip
sudo apt-get install python-pip

###matplotlib
pip install matplotlib
