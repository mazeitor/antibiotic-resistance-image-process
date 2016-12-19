# antibiotic-resistance-process

The aim of this application is extract automatically the antibiotic resistance given a 96-well plate used in sequencer machine.

##execution:
python antibiotic_resistance.py --image images/\<platename\>.png

###input:
images/plate.png with a plate and ninety six wells

###output:
images/\<platename\>/putputXXX.png image with extracted wells
images/<platename>/<row>_<column>_<resistance>_<density>.png cropped image of extracted well
images/<platename>/report.json json with extracted antibiotic resistance for each well
images/<platename>/log.txt log 

row: row index
column: colmun index
resistance: absolute resistance found
density: density of the resistance found

ex: 4-A_122-0.23
is the well 4-A, with 122 pixels found as resistance with density of 23%

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
