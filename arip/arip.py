'''
Copyright 2016, Oriol Mazariegos Canellas <oriol.mazariegos@gmail.com> 
 
This file is part of the ARIP application.
ARIP is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
ARIP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more denode_bs.
You should have received a copy of the GNU General Public License
along with ARIP.  If not, see <http://www.gnu.org/licenses/>.
'''

'''
Created on 01 August 2016
@author: oriol mazariegos
@copyright: Copyright 2016, ARIP
@credits: ["Oriol Mazariegos"]
@license: GPLv3
@version: 1.0.0
@maintainer: oriol mazariegos
@email: oriol.mazariegos@gmail.com
@status: beta
'''

import argparse
import cv2
import math
import numpy as np
import json
import os
import shutil
from numpy.random import randn
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform

def sorting(wells, SORTINGINDEX):
	'''
	@brief: sorting wells structure by index -> rows or columns
	@param wells: numpyarray as matrix structure to manage wells
	@param SORTINGINDEX: index to specify the space to sort
	@return: wells sorted by index
	'''
	i = np.lexsort((wells[:,1], wells[:,SORTINGINDEX])) ###get sorted indexes
	return wells[i]

def clustering(wells, INDEX, THRESHOLD):
	'''
	@brief: clustering wells structure by index -> rows or columns
	@param wells: numpyarray as matrix structure to manage wells
	@param INDEX: index to specify the space to label
	@param THRESHOLD: 
	@return: wells labeled by index
	'''
	label = 0
	idx = 0
	idx += 1
	while idx < len(wells):
		x = wells[idx][INDEX]
		if (x-wells[idx-1][INDEX])>THRESHOLD:
			label += 1
		wells[idx][INDEX+3] = label
		idx += 1
	return wells

def indexing(wells, INDEX):
	'''
	@brief: change cluster id for consecutives id's in rows and columns.
	@param wells: numpyarray as matrix structure to manage wells
	@param INDEX: index to specify the space to label
	@return: wells labeled by index
	'''
	idx = 0
	dic = {}
	for i in np.unique(wells[:,INDEX]):
		dic[i] = idx
		idx += 1
	idx = 0
	while idx < len(wells):
		key = wells[idx][INDEX]
		newkey = dic[key]
		wells[idx][INDEX] = newkey
		idx += 1	
	return wells
	
def removing(wells, INDEX, NUMBER):
	'''
	@brief: remove associates elements with labels less than total number of ROWS or COLUMNS
	@param wells: numpyarray as matrix structure to manage wells
	@param INDEX: index to specify the space to remove
	@param NUMBER: 
	@return: wells cleaned by index
	'''
	histogram = {}
	for points in wells:
		key = points[INDEX]	
		if key in histogram:
			histogram[key] += 1
		else:
			histogram[key] = 1

	cleaned = np.array([])
	cleaned.resize(0,5)
	idx = 0
	while idx < len(wells):
		key = wells[idx][INDEX]
		if histogram[key] >= NUMBER:
			s=cleaned.shape
			cleaned.resize(s[0]+1,5)
			cleaned[s[0]]=wells[idx]
		idx += 1
	return cleaned
			
def antibioticextraction(image, radius):
	'''
	@brief: well segmentation and extraction the antitiobic resistance for that specific well
	@param image: image with a well
	@param radius: radius of the well
	@return: well segmented image, resistance antitiotic metric and total evaluated pixels
	'''
	frame=0
	total = int(np.pi*((radius-frame)*(radius-frame)))
	resistance=0
	x = image.shape[0]
	y = image.shape[1]
	i=0
	while i < x:
		j=0
		while j < y:
			dist = math.hypot(i - y/2, j - x/2)
			if dist >= (radius-frame):	
				image.itemset((i,j),255)
				image.itemset((i,j),255)
				image.itemset((i,j),255)	
			if dist < (radius-frame) and image[i,j] == 0:
				resistance+=1 
			j+=1
		i+=1
	return image,resistance,total

def segmentation(image,wells,radius):
	LABELROWS=[]
	LABELCOLUMNS=[]

	for lr in range(SHAPE[1]):
		LABELROWS.append(str(lr))
	for lc in range(SHAPE[0]):
		LABELCOLUMNS.append(str(chr(lc+65)))

        croppedwells=[]

	for well in wells:
		y = int(well[0])
		x = int(well[1])
		labelY = int(well[3])
		labelX = int(well[4])

		dimension = np.shape(image)
		x1=x-radius
		x2=x+radius
		if x1<0:
			x1=1
		if x2>dimension[0]:
			x2=dimension[0]-1
		y1=y-radius
		y2=y+radius
		if y1<0:
			y1=1
		if y2>dimension[1]:
			y2=dimension[1]-1
	
		cropped = image[x1:x2,y1:y2]
		croppedgray = cv2.cvtColor(cropped, cv2.COLOR_BGR2GRAY)
		
		# Otsu's thresholding after Gaussian filtering
		blur = cv2.GaussianBlur(croppedgray,(5,5),0)		
		threshold,img = cv2.threshold(blur,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
		#here we have img segmented, we can crop the circle
		img,resistance,total = antibioticextraction(img,radius)
		croppedwells.append({"image":img,"row":LABELROWS[labelX],"column":LABELCOLUMNS[labelY],"resistance":resistance,"total":total})
	return croppedwells
	
def quality(wells):
	'''
	@brief: check if all wells have been found or if are missing or more than ninety six
	@param wells: numpyarray as matrix structure to manage wells
	@return: true or false
	'''
	##count wells per each column, should be 8 per column, total 12 columns.
	error = True

	columns = dict()
	rows = dict()

	for i in range(SHAPE[0]):
		columns[i] = 0
	for i in range(SHAPE[1]):
		rows[i] = 0

	for well in wells:
		labelY = int(well[3])
		labelX = int(well[4])

		columns[labelY] = columns[labelY]+1
		rows[labelX] = rows[labelX]+1

	for column in columns:
		if column != SHAPE[1]:
			return False
	for row in rows:
		if row != SHAPE[0]:
			return False

	return error

def normalizingradius(wells,normalizingerror):
        radiusavg = int(np.mean(wells, axis=0)[2])-normalizingerror
	wells[:,2] = radiusavg
	return wells

def paint(wells,output,output_name,platename):
	'''
	@brief: paint in an output image the wells found combinationing the input image
	@param wells: numpyarray as matrix structure to manage wells
	@param output: output image with original image
	@param output_name: output name for the image to write
	'''
	path = "output/{0}".format(platename)

	# ensure at least some wells were found
	if wells is not None:
		# convert the (x, y) coordinates and radius of the wells to integers
		wells = np.round(wells[:,0:3]).astype("int")

		# loop over the (x, y) coordinates and radius of the wells
		for (x, y, r) in wells:
			# draw the circle in the output image, then draw a rectangle
			# corresponding to the center of the circle
			cv2.circle(output, (x, y), r, (0, 255, 0), 4)
			cv2.rectangle(output, (x - 5, y - 5), (x + 5, y + 5), (0, 128, 255), -1)
	 
		# write to file
		filename = "{0}/{1}.jpg".format(path,output_name)
		cv2.imwrite(filename, output, [int(cv2.IMWRITE_JPEG_QUALITY), 90])

def paintcoord(well_x, well_y, radius, output, output_name):
	'''
	@brief: given a rows and columns coordinates paint them in an output image combinationing the input image
	@param well_x: set of wells rows coordinates
	@param well_y: set of wells columns coordinates
	@param radius: radius of the well  
	@param output: output image with original image
	@param output_name: output name for the image to write
	'''
	# loop over the (x, y) coordinates and radius of the wells
	for x in well_x.tolist()[0]:
		x = int(x)
		for y in well_y.tolist()[0]:
			y = int(y)
			# draw the circle in the output image, then draw a rectangle
			# corresponding to the center of the circle
			r = radius
			cv2.circle(output, (x, y), r, (0, 255, 0), 4)
			cv2.rectangle(output, (x - 5, y - 5), (x + 5, y + 5), (0, 128, 255), -1)

	# write to file
	filename = "images/{0}.jpg".format(output_name)
	cv2.imwrite(filename, output, [int(cv2.IMWRITE_JPEG_QUALITY), 90])
	
def write(wells, platename):
	'''
	@brief: write each well found as a file
	@param wells: numpyarray as matrix structure to manage wells
	@param platename: platename
	@param log: list of log
	'''
	path = "output/{0}".format(platename)

	data=dict()
	for object in wells:
		cropped = object["image"]
		row = object["row"]
		column = object["column"]
		resistance = object["resistance"]
		total = object["total"]
		density = round(float(resistance)/float(total),2)
		filename = "{0}/{1}-{2}_{3}-{4}.jpg".format(path,row,column,resistance,density)
		cv2.imwrite(filename, cropped, [int(cv2.IMWRITE_JPEG_QUALITY), 90])

		object["density"] = density
		object.pop("image",None)
		data[row+"-"+column]=object
	
	filename = "{0}/{1}.json".format(path,"report")
	jsonfile = open(filename, 'wb')
	json.dump(data,jsonfile)


def writeLog(platename,log):
	path = "output/{0}".format(platename)
	filename = "{0}/{1}.txt".format(path,"log")
	logfile = open(filename, 'wb')
	for line in log:
		logfile.write(line)
		logfile.write("\n")
	logfile.close()
	


def execution1(image, outputs, minRadius, maxRadius, clusterthreshold, platename):
       	###circle detection
	try: 
		##convert image to grayscale
		gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
		gray = cv2.medianBlur(gray,3)

		##detect wells in the image
		wells = cv2.HoughCircles(
						gray,                                           ## image
						cv2.cv.CV_HOUGH_GRADIENT,       ## method for detecting wells
						1,                                                     ## canny filter
						(minRadius+maxRadius)/2,
						param1=30,                                      ## 30
						param2=20,                                      ## 15
						minRadius=minRadius,            ## 20 ## 18 relaxed
						maxRadius=maxRadius             ## 25 ## 27 relaxed
					  )

		LABELTHRESHOLD_INDEX = clusterthreshold
		
		##preparing well data structure
		wells = wells[0,:]
		wells = np.insert(wells, 3, 0, axis=1) ## add dimensionality
		wells = np.insert(wells, 4, 0, axis=1) ## add dimensionality
	
		##writting wells with original
		paint(wells,outputs["output1"],"output1",platename)

		##column 0,1 speficy the COORDINATES and COLUMN 2,3 identify ROW and COLUMN label
		##sorting and clustering objects by rows
		wells = sorting(wells, ROW_INDEX)
		wells = clustering(wells, ROW_INDEX, LABELTHRESHOLD_INDEX)

		##sorting and clustering objects by columns
		wells = sorting(wells, COLUMN_INDEX)
		wells = clustering(wells, COLUMN_INDEX, LABELTHRESHOLD_INDEX)

		error = quality(wells)
		return error,len(wells),wells
	except:
		return False,0,None

def execution2(image, outputs, wells, platename):
	###cleaning wells
	##remove associates elements with labels less than total number of ROWS or COLUMNS
	wells = removing(wells, ROW_INDEX+3, NUM_LABELS_IN_ROWS)
	wells = removing(wells, COLUMN_INDEX+3, NUM_LABELS_IN_COLUMNS)

	##writting wells with original
	paint(wells,outputs["output2"],"output2",platename)

	##indexing cluster id for consecutives id's in rows and columns
	wells = indexing(wells,ROW_INDEX+3)
	wells = indexing(wells,COLUMN_INDEX+3)

	##here we can check if we have 96 samples then we analyse the image otherwise we analyse but warning with message
	error = quality(wells)
	return error, len(wells), wells


def execution3(image, outputs, normalizingerror, wells, platename):
        ###segmentation wells

        ##getting an average of the radius
        radiusavg = int(np.mean(wells, axis=0)[2])-normalizingerror
	wells = normalizingradius(wells,normalizingerror)
        ##writting wells with original
        paint(wells,outputs["output3"],"output3", platename)

        ##given a matrix of samples and an average radius, aply an otsu segmentation and get only the wells, not bounding box
        wells = segmentation(image,wells,radiusavg)
       
	return wells

def process(args):
	log = []

	input_path = args["image"]
	base, platename = os.path.split(input_path)
	platename, extension = os.path.splitext(platename)
	minRadius = 18
	maxRadius = 23 ##scale wells recognizing 
	normalizingerror = 4 
	clusterthreshold = 2 ##cluster wells recognizing
	if args["minRadius"] is not None:
		minRadius = int(args["minRadius"])
	if args["maxRadius"] is not None:
		maxRadius = int(args["maxRadius"])
	if args["normError"] is not None:
		normalizingerror = int(args["normError"])
	if args["threshold"] is not None:
		clusterthreshold = int(args["threshold"])
	if args["shape"] is not None:
		global SHAPE
		SHAPE = args["shape"]

	##preparing folder for outputs
        path = "output/{0}".format(platename)
	if not os.path.exists(path):
		os.makedirs(path)
	else:
	        shutil.rmtree(path) #removes all the subdirectories!
	        os.makedirs(path)

	##load the image, clone it for output
	image = cv2.imread(input_path)
	outputs = dict()
	outputs["output1"]=image.copy()
	outputs["output2"]=image.copy()
	outputs["output3"]=image.copy()

	##dynamic algorithm to find the result that converge, in this case we take account the maxradius of a well to be robust in scale and threshold to recognize a well is inside a grid (96-well plate is a grid with 8 row and 12 columns)
	MAXWELLS = SHAPE[0] * SHAPE[1]
	
	numwells = 0
	thresholditerations = 5	#5
	clusterthreshold = 5	#5
	while numwells != MAXWELLS and thresholditerations>0:
		iterations = 10 
		maxRadius = 23 

		while numwells < MAXWELLS and iterations>0:
			error, numwells, wells = execution1(image, outputs, minRadius, maxRadius, clusterthreshold, platename)
			log.append("customizing scale well: found {0}, num wells {1}, min radius value {2}, max radius value {3}, clusterthreshold {4}".format(error, numwells, minRadius, maxRadius,clusterthreshold))
			maxRadius = maxRadius + 1
			iterations = iterations - 1
			
			if numwells>=MAXWELLS:
				error, numwells, wells = execution2(image, outputs, wells, platename)
				log.append("customizing grid matching: found {0}, num wells recognized {1}".format(error, numwells))
		
		clusterthreshold = clusterthreshold - 1
		thresholditerations = thresholditerations - 1

	if numwells == MAXWELLS:
		wells = execution3(image, outputs, normalizingerror, wells, platename)
		log.append("Succesfully processed plate, found 96 wells")
	else:
		log.append("No processed plate, not found 96 wells")
		wells = []
		

	if len(wells):
        	##write the results in a separated file
		write(wells,platename)

	writeLog(platename,log)


##initializing variables
NUM_LABELS_IN_ROWS = 4
NUM_LABELS_IN_COLUMNS = 3
ROW_INDEX = 0
COLUMN_INDEX = 1
SHAPE = (12,8)


if __name__ == '__main__':
	# construct the argument parser and parse the arguments
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--image", required = True, help = "Path to the image")
	ap.add_argument("--minRadius", required = False, help = "Min radius for hough circles method")
	ap.add_argument("--maxRadius", required = False, help = "Max radius for hough circles method")
	ap.add_argument("--normError", required = False, help = "Error value to normalize the radius of a circle to be evaluated")
	ap.add_argument("--threshold", required = False, help = "Threshold to construct groups of rows and colums")
	ap.add_argument("--shape", required = False, help = "Shape of the grid to be recognized")

	args = vars(ap.parse_args())

	process(args)
