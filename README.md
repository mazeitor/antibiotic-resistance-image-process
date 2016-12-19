# antibiotic-resistance-process

The aim of this application is to extract automatically the antibiotic resistance given a 96-well plate used in a sequencer machine. The 96 wells are cultured with a bacteria combined with an antibiotic. The bacteria which is resistance to that single antibiotic, grow up and can recognised inside the well. The applications use computer vision methods to process this image and as a result construct a report for the collection of the 96 wells, and as a main feature give the density of the bacteria which have grown up inside each the well.


###execution:
python antibiotic_resistance.py --image images/\<platename\>.png

###input:
images/\<platename\>.png with a plate and ninety six wells

###output:
images/\<platename\>/outputXXX.png image with extracted wells
images/\<platename\>/\<row\>_\<column\>_\<resistance\>_\<density\>.png cropped image of extracted well
images/\<platename\>/report.json json with extracted antibiotic resistance for each well
images/\<platename\>/log.txt log 

Description of the schema:
row: row index
column: colmun index
resistance: absolute resistance found
density: density of the resistance found

report example:
```
   "7-J":{  
      "density":0.17,
      "column":"A",
      "resistance":122,
      "total":706,
      "row":"4"
   },
```
output images example:
```  
4-A_122-0.23, is the well 4-A, with 122 pixels found as resistance with density of 17%
```
output log example:
```
customizing scale well: found False, num wells 93, min radius value 18, max radius value 23
customizing scale well: found False, num wells 96, min radius value 18, max radius value 24
customizing grid matching: found False, num wells recognized 96
Succesfully processed plate, found 96 wells
```

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
