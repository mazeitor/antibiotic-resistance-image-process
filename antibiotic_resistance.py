# import the necessary packages
import numpy as np
import argparse
import cv2
import pdb

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
from numpy.random import randn
import math

def sorting(circles, SORTINGINDEX):
	i = np.lexsort((circles[:,1], circles[:,SORTINGINDEX])) ###get sorted indexes
	return circles[i]

def labeling(circles, INDEX, THRESHOLD):
	label = 0
	idx = 0
	idx += 1
	while idx < len(circles):
		x = circles[idx][INDEX]
		if (x-circles[idx-1][INDEX])>THRESHOLD:
			label += 1
		circles[idx][INDEX+3] = label
		idx += 1
	return circles

def changingLabels(circles, INDEX):
	idx = 0
	dic = {}
	for i in np.unique(circles[:,INDEX]):
		dic[i] = idx
		idx += 1
	idx = 0
	while idx < len(circles):
		key = circles[idx][INDEX]
		newkey = dic[key]
		circles[idx][INDEX] = newkey
		idx += 1	
	return circles
	
def removing(circles, INDEX, NUMBER):
    histogram = {}
    for points in circles:
		key = points[INDEX]	
		if key in histogram:
			histogram[key] += 1
		else:
			histogram[key] = 1

    cleaned = np.array([])
    cleaned.resize(0,5)
    idx = 0
    while idx < len(circles):
		key = circles[idx][INDEX]
		if histogram[key] >= NUMBER:
			s=cleaned.shape
			cleaned.resize(s[0]+1,5)
			cleaned[s[0]]=circles[idx]
		idx += 1
    return cleaned

def normalizingCoordinates(circle):
	for circle in circles:
		r = circle[3]
		c = circle[4]
		x[c,r] = circle[0]
		y[r,c] = circle[1]

	xavg = np.mean(x, axis=0)
	yavg = np.mean(y, axis=0)
	
	return xavg,yavg

			
def segmentedCircle(image, radius):
	frame = 5
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
		
def getSegmentedCircles(image, rows, columns, radius):
	
	LABELROWS=["1","2","3","4","5","6","7","8"]
	LABELCOLUMNS=["A","B","C","D","E","F","G","H","I","J","K","l"]

	croppedcircles=[]
	
	labelX = 0
	for x in rows.tolist()[0]:
		labelY = 0
		x = int(x)
		for y in columns.tolist()[0]:	
			y=int(y)			
			cropped = image[x-radius:x+radius,y-radius:y+radius]
			croppedgray = cv2.cvtColor(cropped, cv2.COLOR_BGR2GRAY)
			
			# Otsu's thresholding after Gaussian filtering
			blur = cv2.GaussianBlur(croppedgray,(5,5),0)
			threshold,img = cv2.threshold(blur,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
			##here we have img segmented, we can crop the circle
			img,resistance,total = segmentedCircle(img,radius)
		
			croppedcircles.append({"image":img,"row":LABELROWS[labelX],"column":LABELCOLUMNS[labelY],"resistance":resistance,"total":total})
			
			labelY += 1
		labelX += 1
	
	return croppedcircles
	
def paint(circles,output,output_name):
	# ensure at least some circles were found
	if circles is not None:
		# convert the (x, y) coordinates and radius of the circles to integers
		circles = np.round(circles[:,0:3]).astype("int")

		# loop over the (x, y) coordinates and radius of the circles
		for (x, y, r) in circles:
			# draw the circle in the output image, then draw a rectangle
			# corresponding to the center of the circle
			cv2.circle(output, (x, y), r, (0, 255, 0), 4)
			cv2.rectangle(output, (x - 5, y - 5), (x + 5, y + 5), (0, 128, 255), -1)
	 
		# show the output image
		# cv2.imshow("output", np.hstack([image, output]))
		# cv2.waitKey(0)
		filename = "images\{0}.jpg".format(output_name)
		cv2.imwrite(filename, output, [int(cv2.IMWRITE_JPEG_QUALITY), 90])

def paintcoord(circle_x, circle_y, radius, output, output_name):
	# loop over the (x, y) coordinates and radius of the circles
	for x in circle_x.tolist()[0]:
		x = int(x)
		for y in circle_y.tolist()[0]:
			y = int(y)
			# draw the circle in the output image, then draw a rectangle
			# corresponding to the center of the circle
			r = circleavg
			cv2.circle(output, (x, y), r, (0, 255, 0), 4)
			cv2.rectangle(output, (x - 5, y - 5), (x + 5, y + 5), (0, 128, 255), -1)

	# show the output image
	cv2.imshow("output", np.hstack([image, output]))
	cv2.waitKey(0)
	filename = "images\{0}.jpg".format(output_name)
	cv2.imwrite(filename, output, [int(cv2.IMWRITE_JPEG_QUALITY), 90])
	
def write(circles):
	for object in circles:
		cropped = object["image"]
		row = object["row"]
		column = object["column"]
		resistance = object["resistance"]
		total = object["total"]
		density = round(float(resistance)/float(total),2)
		filename = "output\s{0}-{1}_{2}-{3}.jpg".format(row,column,resistance,density)
		cv2.imwrite(filename, cropped, [int(cv2.IMWRITE_JPEG_QUALITY), 90])
	
	
# construct the argument parser and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--image", required = True, help = "Path to the image")
ap.add_argument("--minRadius", required = False, help = "Min radius for hough circles method")
ap.add_argument("--maxRadius", required = False, help = "Max radius for hough circles method")
args = vars(ap.parse_args())

input_path = args["image"]
minRadius = 18
maxRadius = 27
if args["minRadius"] is not None:
	minRadius = int(args["minRadius"])
if args["maxRadius"] is not None:
	maxRadius = int(args["maxRadius"])

# load the image, clone it for output, and then convert it to grayscale
image = cv2.imread(input_path)
output1 = image.copy()
output2 = image.copy()
output3 = image.copy()


gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
gray = cv2.medianBlur(gray,3)

# detect circles in the image
circles = cv2.HoughCircles(
				gray, 						## image
				cv2.cv.CV_HOUGH_GRADIENT, 	## method for detecting circles
				1, 							## canny filter
				15,							## minimun distance between circles detected
				param1=30,					## 30
				param2=15,					## 15
				minRadius=minRadius, 		## 20 ## 18 relaxed
				maxRadius=maxRadius 		## 25 ## 27 relaxed
			  )


		  
#initializing variables
NUM_LABELS_IN_ROWS = 4
NUM_LABELS_IN_COLUMNS = 3
LABELTHRESHOLD_INDEX = 3
ROW_INDEX = 0
COLUMN_INDEX = 1

##preparing circle data structure
circles = circles[0,:]
circles = np.insert(circles, 3, 0, axis=1) ## add dimensionality
circles = np.insert(circles, 4, 0, axis=1) ## add dimensionality

##writting circles with original
paint(circles,output1,"output1")

##column 0,1 speficy the COORDINATES and COLUMN 2,3 identify ROW and COLUMN label
##sorting and labeling objects by rows
circles = sorting(circles, ROW_INDEX) 
circles = labeling(circles, ROW_INDEX, LABELTHRESHOLD_INDEX) 

##sorting and labeling objects by columns
circles = sorting(circles, COLUMN_INDEX)
circles = labeling(circles, COLUMN_INDEX, LABELTHRESHOLD_INDEX) 
 
##remove associates elements with labels less than total number of ROWS or COLUMNS
circles = removing(circles, ROW_INDEX+3, NUM_LABELS_IN_ROWS) 
circles = removing(circles, COLUMN_INDEX+3, NUM_LABELS_IN_COLUMNS) 

##writting circles with original
paint(circles,output2,"output2")

##here we can check if we have 96 samples then we analyse the image otherwise we analyse but warning with message

##changing labels for consecutives labels in rows and columns
circles = changingLabels(circles,ROW_INDEX+3)
circles = changingLabels(circles,COLUMN_INDEX+3)

##normalizing coordinates, 8 rows and 12 columns
x = np.matrix(np.arange(8*12).reshape((8, 12)))
y = np.matrix(np.arange(8*12).reshape((12, 8)))
columns,rows = normalizingCoordinates(circles)
circleavg = int(np.mean(circles, axis=0)[2])-3	

##writting normalized circles with original
paintcoord(x,y,circleavg,output3,"output3")

##given a matrix of samples and an average radius, aply a otsu segmentation and get only the circles, not bounding box
circles = getSegmentedCircles(image, rows, columns, circleavg)

##write the results in a separated file
write(circles)	