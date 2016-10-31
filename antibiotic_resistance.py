import argparse
import cv2
import math
import numpy as np
from numpy.random import randn
from matplotlib import pyplot as plt
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

def segmentation(image,wells,radius):
        LABELROWS=["1","2","3","4","5","6","7","8"]
        LABELCOLUMNS=["A","B","C","D","E","F","G","H","I","J","K","l"]

        croppedwells=[]

	#index = np.lexsort((wells[:,4],wells[:,3]))
	#wells = wells[index]
	for well in wells:
		y = int(well[0])
		x = int(well[1])
		labelY = int(well[3])
		labelX = int(well[4])

		cropped = image[x-radius:x+radius,y-radius:y+radius]
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
	@brief: check if there all wells have been found or if are missing or more than ninety six
	@param wells: numpyarray as matrix structure to manage wells
	@return: true or false
	'''
	return len(wells) == 96
	
def paint(wells,output,output_name):
	'''
	@brief: paint in an output image the wells found combinationing the input image
	@param wells: numpyarray as matrix structure to manage wells
	@param output: output image with original image
	@param output_name: output name for the image to write
	'''
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
		filename = "images/{0}.jpg".format(output_name)
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
	
def write(wells):
	'''
	@brief: write each well found as a file
	@param wells: numpyarray as matrix structure to manage wells
	'''
	for object in wells:
		cropped = object["image"]
		row = object["row"]
		column = object["column"]
		resistance = object["resistance"]
		total = object["total"]
		density = round(float(resistance)/float(total),2)
		filename = "output/{0}-{1}_{2}-{3}.jpg".format(row,column,resistance,density)
		cv2.imwrite(filename, cropped, [int(cv2.IMWRITE_JPEG_QUALITY), 90])
	
if __name__ == '__main__':
	# construct the argument parser and parse the arguments
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--image", required = True, help = "Path to the image")
	ap.add_argument("--minRadius", required = False, help = "Min radius for hough circles method")
	ap.add_argument("--maxRadius", required = False, help = "Max radius for hough circles method")
	ap.add_argument("--normError", required = False, help = "Error value to normalize the radius of a circle to be evaluated")
	ap.add_argument("--threshold", required = False, help = "Threshold to construct groups of rows and colums")

	args = vars(ap.parse_args())

	input_path = args["image"]
	minRadius = 18
	maxRadius = 27
	normalizingerror = 3
	labelthreshold = 3
	if args["minRadius"] is not None:
		minRadius = int(args["minRadius"])
	if args["maxRadius"] is not None:
		maxRadius = int(args["maxRadius"])
	if args["normError"] is not None:
		normalizingerror = int(args["normError"])
	if args["threshold"] is not None:
		labelthreshold = int(args["threshold"])

	##load the image, clone it for output
	image = cv2.imread(input_path)
	outputs = dict()
	outputs["output1"]=image.copy()
	outputs["output2"]=image.copy()
	outputs["output3"]=image.copy()

	##convert image to grayscale
	gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
	gray = cv2.medianBlur(gray,3)

	##detect wells in the image
	wells = cv2.HoughCircles(
					gray, 						## image
					cv2.cv.CV_HOUGH_GRADIENT, 	## method for detecting wells
					1, 							## canny filter
					15,							## minimun distance between wells detected
					param1=30,					## 30
					param2=15,					## 15
					minRadius=minRadius, 		## 20 ## 18 relaxed
					maxRadius=maxRadius 		## 25 ## 27 relaxed
				  )
		  
	##initializing variables
	NUM_LABELS_IN_ROWS = 4
	NUM_LABELS_IN_COLUMNS = 3
	LABELTHRESHOLD_INDEX = labelthreshold
	ROW_INDEX = 0
	COLUMN_INDEX = 1

	##preparing well data structure
	wells = wells[0,:]
	wells = np.insert(wells, 3, 0, axis=1) ## add dimensionality
	wells = np.insert(wells, 4, 0, axis=1) ## add dimensionality

	##writting wells with original
	paint(wells,outputs["output1"],"output1")

	##column 0,1 speficy the COORDINATES and COLUMN 2,3 identify ROW and COLUMN label
	##sorting and clustering objects by rows
	wells = sorting(wells, ROW_INDEX) 
	wells = clustering(wells, ROW_INDEX, LABELTHRESHOLD_INDEX) 

	##sorting and clustering objects by columns
	wells = sorting(wells, COLUMN_INDEX)
	wells = clustering(wells, COLUMN_INDEX, LABELTHRESHOLD_INDEX) 
	 
	error = quality(wells)
	print error,len(wells)

	##remove associates elements with labels less than total number of ROWS or COLUMNS
	wells = removing(wells, ROW_INDEX+3, NUM_LABELS_IN_ROWS) 
	wells = removing(wells, COLUMN_INDEX+3, NUM_LABELS_IN_COLUMNS) 

	##writting wells with original
	paint(wells,outputs["output2"],"output2")

	##indexing cluster id for consecutives id's in rows and columns
	wells = indexing(wells,ROW_INDEX+3)
	wells = indexing(wells,COLUMN_INDEX+3)

	##here we can check if we have 96 samples then we analyse the image otherwise we analyse but warning with message
	error = quality(wells)
	print error, len(wells)

	#print error , type(wells), wells
	##getting an average of the radius
	radiusavg = int(np.mean(wells, axis=0)[2])-normalizingerror	
	
	##writting normalized wells with original
	#paintcoord(x,y,radiusavg,outputs["output3"],"output3")

	##given a matrix of samples and an average radius, aply a otsu segmentation and get only the wells, not bounding box
	wells = wellSegmentation(image,wells,radiusavg)
	
	##write the results in a separated file
	write(wells)	
