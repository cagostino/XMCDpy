#lets gooooo. Author: Christopher Agostino
import numpy as np
import os
import osgeo
import scipy
import gdal
import matplotlib.pyplot as plt
from math import *
from matplotlib.backend_bases import MouseEvent
import pyfits
def assign_images(lst = None):
	new_array =[]
	file_list =[]
	for fil in os.listdir("."):
		if fil.endswith(".fits"):
			new_array.append(pyfits.open(fil)[0].data)
			file_list.append(fil)
		elif fil.endswith(".tif"):
			new =gdal.Open(fil)

			new_array.append(scipy.array(new.GetRasterBand(1).ReadAsArray()))
			file_list.append(fil)
	return [file_list, new_array]
list_images = assign_images()[1]
list_files = assign_images()[0]
photon_energies = []
for i in list_files:
	if 'fits' in i:

		i = i.replace(".fits","")
	if 'tif' in i:
		i = i.replace("_L.tif","")
	photon_energies.append(float(i))
def pick_domain_wall_point(image):
	plt.imshow(image,cmap='Greys_r')
	coords = plt.ginput(1)
	coords = coords[0]
	return np.array(coords) 
#point = pick_domain_wall_point(list_images[23])
def calculate_nth_average_intensity(image,coordinate, size, number):
	intensities_right = []
	intensities_left = []
	output_array = []
	count_x = 0
	while count_x < size:
		intensities_right.append(image[coordinate[1]][coordinate[0]+count_x+size*(number-1)])
		count_y = 1
		while count_y < size:
			intensities_right.append(image[coordinate[1]+count_y][coordinate[0]+count_x+size*(number-1)])
			count_y+=1
		count_x +=1 
	count_x =0
	while count_x > -size:
                intensities_left.append(image[coordinate[1]][coordinate[0]+count_x-size*(number-1)])
                count_y = 1
                while count_y < size:
                        intensities_left.append(image[coordinate[1]+count_y][coordinate[0]+count_x-size*(number-1)])
                        count_y+=1
                count_x -=1
	return [np.sum(intensities_left)/(abs(count_x)*count_y)-np.sum(intensities_right)/(abs(count_x)*count_y), np.sum(intensities_left)/(abs(count_x)*count_y), np.sum(intensities_right)/(abs(count_x)*count_y) ]
#point1=[153,774]
point =[int(i) for i in list(raw_input("enter an x and y coordinate: ").split(','))]
pixel_list = [5,8,10,15,20]
def producedata(image_list, pixel_list,number):
	for pixel in pixel_list:
		for num in range(number):
			intens =[]
			list_intensities= []
			intensleft = []
			intensright = []
			for i in range(len(image_list)):
				list_intensities.append(calculate_nth_average_intensity(image_list[i],point,pixel,num))
			for j in range(len(list_intensities)):
				intens.append(list_intensities[j][0])
				intensleft.append(list_intensities[j][1])
				intensright.append(list_intensities[j][2])
			fileref= open(str(point[1])+"_"+str(point[0])+"_" + str(num) +"_area_"+str(pixel) +".csv" ,'w')
			for k in range(len(photon_energies)):
				fileref.write(str(photon_energies[k]) + " " + str(intens[k]) + " " + str(intensleft[k]) +" " +str(intensright[k]) +"\n")
			fileref.close()

producedata(list_images,pixel_list,5	)