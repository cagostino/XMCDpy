#lets gooooo. Author: Christopher Agostino
import numpy as np
import os
import osgeo
import scipy
import gdal
import matplotlib.pyplot as plt
from lmfit.models import LorentzianModel
from lmfit.models import GaussianModel
from lmfit.models import VoigtModel
#import pyfits
import gc
plt.rc('text',usetex=True)
#make these plots look dank
plt.rc('font',family='serif')
l = os.listdir("left/")
r =os.listdir("right/")
l.sort()
r.sort()
from numpy import trapz
from scipy.integrate import simps
def read_dat(filename):
	fil = np.loadtxt(filename, skiprows=1)
	return fil.transpose()
def assign_images(lst = None):
	new_array_r =[]
	file_list_r =[]
	new_array_l =[]
	file_list_l =[]
	for fil in l:
		if fil.endswith(".fits"):
			new_array_l.append(pyfits.open('left/'+fil)[0].data)
			file_list_l.append(fil)
		elif fil.endswith(".tif"):
			new =gdal.Open('left/'+fil)

			new_array_l.append(scipy.array(new.GetRasterBand(1).ReadAsArray()))
			file_list_l.append(fil)
	for fil in r:
		if fil.endswith(".fits"):
			new_array_r.append(pyfits.open('right/'+fil)[0].data)
			file_list_r.append(fil)
		elif fil.endswith(".tif"):
			new =gdal.Open('right/'+fil)

			new_array_r.append(scipy.array(new.GetRasterBand(1).ReadAsArray()))
			file_list_r.append(fil)	
	return [file_list_l,file_list_r, new_array_l,new_array_r]
list_images_l = assign_images()[2]
list_images_r = assign_images()[3]

list_files_l = assign_images()[0]
list_files_r = assign_images()[1]

photon_energies = []
for i in list_files_l:
	if 'fits' in i:
		#if '_L' in i:
		#	i = i.replace("_L.fits","")
		if '_R' in i:
			i = i.replace("_R.fits","")
		else:
		
			i=i.replace(".fits","")
		if "_" in i:
			i=i.replace("_",'.')
	if 'tif' in i:
		#if '_L' in i:
		#	i = i.replace("_L.tif","")
		if '_R' in i:
			i = i.replace("_R.tif","")	
		else:
			i=i.replace(".tif",'')
		if "_" in i:
			i=i.replace("_",'.')


	photon_energies.append(float(i))
gc.collect()

def pick_domain_wall_point(image):
	plt.imshow(image,cmap='Greys_r')
	coords = plt.ginput(1)
	coords = coords[0]
	return np.array(coords) 
#point = pick_domain_wall_point(list_images[23])
def calculate_nth_average_intensity(image_left,image_right,coordinate, size, number,direct=1,det=0,phot_e=0):
	intensities_right = []
	intensities_left = []
	output_array = []
	count_x = 0
	norm_l =[]
	x, y = coordinate
	norm_r =[]
	#determine best normalization given a point
	
	diff_1 =np.array( [])
	diff_2 = np.array([])
	if det:	
		for a in range(coordinate[0]-5,coordinate[0]+5):
			left_1 =[]
			right_1 =[]
			left_2 =[]
			right_2 =[]
			for j in range(size+1):
				for i in range(size):
					left_1.append(image_left[y+i][a-j])
					right_1.append(image_left[y+i][a+j])
					
					left_2.append(image_right[y+i][a-j])
					right_2.append(image_right[y+i][a+j])

			left_1, right_1 = np.array(left_1), np.array(right_1)
			left_1, right_1 = np.mean(left_1), np.mean(right_1)
			m0_1 = (left_1+right_1)/2
			left_2, right_2 = np.array(left_2), np.array(right_2)
			left_2, right_2 = np.mean(left_2), np.mean(right_2)
			m0_2 = (left_2+right_2)/2
			diff1 = m0_1-image_left[y][a]
			diff_1 = np.append(diff_1,diff1)
			diff2 = m0_2-image_right[y][a]
			diff_2 = np.append(diff_2, diff2)
		choice = np.where( abs(diff_1) + abs(diff_2) == min( abs(diff_1)+abs(diff_2)) )[0][0]
		choice -=5
		#coordinate[0] += choice
		
		x+=choice
		if phot_e == 778.0:
			print choice,x
	norm_l =[]
	norm_r =[]
	for j in range(-size,size+1):
		for i in range(size):
			norm_l.append(image_left[y+i][x+j])
			norm_r.append(image_right[y+i][x+j])
	norm_l, norm_r = np.array(norm_l),np.array(norm_r)
	norm_l =np.mean(norm_l)
	norm_r = np.mean(norm_r)
	back = (norm_l+norm_r)
	while count_x < size:
		intensities_right.append(image_right[y][x+direct*count_x+size*(number-1)])
		count_y = 1
		while count_y < size:
			intensities_right.append(image_right[y+count_y][x+direct*count_x+direct*size*(number-1)])
			count_y+=1
		count_x +=1 
	count_x =0
	while count_x<size:
                intensities_left.append(image_left[y][x+direct*count_x+direct*size*(number-1)])
                count_y = 1
                while count_y < size:
                        intensities_left.append(image_left[y+count_y][x+direct*count_x+direct*size*(number-1)])
      	                count_y+=1
                count_x +=1
	intens_left=np.sum(intensities_left)/abs(count_x*count_y)#-norm_l
	
	intens_right = np.sum(intensities_right)/(abs(count_x)*count_y)#-norm_r
	xmcd= intens_left-intens_right
	intens=(intens_left-norm_l) -(intens_right-norm_r)
	sum_l_r= intens_right+intens_left
	I_plus = sum_l_r/2. + intens/2.
	I_minus = sum_l_r/2. -intens/2.
	return [intens,I_plus, I_minus,intens_left, intens_right,norm_l, norm_r,intens_left-norm_l, intens_right-norm_r] 

#point1=[153,774]
#point =[int(i) for i in list(raw_input("enter an x and y coordinate: ").split(','))]

def lorentzfit(photon_e, xmcd):
	
	mod = LorentzianModel()
	pars = mod.guess(xmcd, x=photon_e)
	out  = mod.fit(xmcd, pars, x=photon_e)
	return out.best_fit
def gauss(photon_e, xmcd):
	mod= GaussianModel()
	pars = mod.guess(xmcd, x=photon_e)
	out=mod.fit(xmcd, pars,x=photon_e)
	return out.best_fit
def voigt(photon_e, xmcd):
	mod =VoigtModel()
	pars = mod.guess(xmcd, x= photon_e)
	out=mod.fit(xmcd, pars, x=photon_e)
	return out.best_fit
def ls_ratio(a3,a2):
	ls=(2*(a3+a2))/(3*(a3-2*a2))
	return ls
pixel_list = [5,8,10,15]
def producedata(image_list_l,image_list_r, pixel_list,number,point,direct=1,det=0):
	totxmcd=[]	
	for pixel in pixel_list:
		for num in range(1,number+1):
			intens =[]
			list_intensities= []
			intensleft = []
			print num
			intensright = []
			back_l = []
			back_r =[]
			norm_l = []
			norm_r =[]
			combo_l=[]
			combo_r=[]
			for i in range(len(image_list_r)):
				if photon_energies[i] == 778.0: print photon_energies[i]
				list_intensities.append(calculate_nth_average_intensity(image_list_l[i],image_list_r[i],point,pixel,num,direct=direct,det=det,phot_e = photon_energies[i]))
			for j in range(len(list_intensities)):
				intens.append(list_intensities[j][0])
				intensleft.append(list_intensities[j][1])
				intensright.append(list_intensities[j][2])
				back_l.append(list_intensities[j][3])
				back_r.append(list_intensities[j][4])
				norm_l.append(list_intensities[j][5])
				norm_r.append(list_intensities[j][6])
				combo_l.append(list_intensities[j][7])
				combo_r.append(list_intensities[j][8])
			#turning them into absorptions
			combo_l = np.array(combo_l)
			combo_l-= np.median(combo_l[0:9])
			combo_r= np.array(combo_r)
			combo_r-= np.median(combo_r[0:9])			
			xmcd = combo_r-combo_l
			totxmcd.append(xmcd)
			abs_l = -np.log(back_l)
			abs_r = -np.log(back_r)
			#writing them to files
			tes,xtes = np.array(xmcd[0:30]), np.array(photon_energies[0:30])
			tes,xtes= np.append(tes,xmcd[-10:]), np.append(xtes, photon_energies[-10:])
			m, yint = np.polyfit(xtes, tes,1)
		
			xmcd = xmcd- np.array(photon_energies) *m-yint #trying to do linear fit
	
			'''mod = LorentzianModel()
			photon_e = np.array(photon_energies)
			pars= mod.guess(xmcd, x=photon_e)
			out= mod.fit(xmcd, pars, x=photon_e)
			lorentzfit_top = out.best_fit'''
			#doing a lorentzian fit on them for the double peaks
			lorentzfit_top = lorentzfit(np.array(photon_energies)[-45:], xmcd[-45:])
			lorentzfit_bot = lorentzfit(np.array(photon_energies[:-40]), xmcd[:-40])

			# integ to l/s with trapz
			integ_lorentz_top_trapz = trapz(lorentzfit_top, x= np.array(photon_energies)[-45:])
			integ_lorentz_bot_trapz = trapz(lorentzfit_bot, x=np.array(photon_energies[:-40]) )
			# integ to l/s with simps
			integ_lorentz_top_simps = simps(lorentzfit_top, x= np.array(photon_energies)[-45:])
			integ_lorentz_bot_simps = simps(lorentzfit_bot, x=np.array(photon_energies[:-40]) )
			
			lorentzfit_bot = np.append(lorentzfit_bot, np.zeros(40))
			lorentzfit_top = np.append(np.zeros(len(photon_energies)-45),lorentzfit_top)
			lorentz_tot = lorentzfit_top + lorentzfit_bot

			ls_trapz_lor = ls_ratio(integ_lorentz_top_trapz,integ_lorentz_bot_trapz)
			ls_simps_lor = ls_ratio(integ_lorentz_top_simps, integ_lorentz_bot_simps)

			gaussfit_top = gauss(np.array(photon_energies)[-45:], xmcd[-45:])
			gaussfit_bot = gauss(np.array(photon_energies[:-40]), xmcd[:-40])
			# integ to l/s with trapz
			integ_gauss_top_trapz = trapz(gaussfit_top, x= np.array(photon_energies)[-45:])
			integ_gauss_bot_trapz = trapz(gaussfit_bot, x=np.array(photon_energies[:-40]) )
			# integ to l/s with simps
			integ_gauss_top_simps = simps(gaussfit_top, x= np.array(photon_energies)[-45:])
			integ_gauss_bot_simps = simps(gaussfit_bot, x=np.array(photon_energies[:-40]) )
			
			gaussfit_bot = np.append(gaussfit_bot, np.zeros(40))
			gaussfit_top = np.append(np.zeros(len(photon_energies)-45),gaussfit_top)
			gaussfit_tot = gaussfit_top + gaussfit_bot
			
			ls_trapz_gauss = ls_ratio(integ_gauss_top_trapz,integ_gauss_bot_trapz)
			ls_simps_gauss = ls_ratio(integ_gauss_top_simps, integ_gauss_bot_simps)

			
			voigtfit_top = voigt(np.array(photon_energies)[-45:], xmcd[-45:])
			voigtfit_bot = voigt(np.array(photon_energies[:-40]), xmcd[:-40])

			# integ to l/s with trapz
			integ_voigt_top_trapz = trapz(voigtfit_top, x= np.array(photon_energies)[-45:])
			integ_voigt_bot_trapz = trapz(voigtfit_bot, x=np.array(photon_energies[:-40]) )
			# integ to l/s with simps
			integ_voigt_top_simps = simps(voigtfit_top, x= np.array(photon_energies)[-45:])
			integ_voigt_bot_simps = simps(voigtfit_bot, x=np.array(photon_energies[:-40]))
			

			voigtfit_bot = np.append(voigtfit_bot, np.zeros(40))
			voigtfit_top = np.append(np.zeros(len(photon_energies)-45),voigtfit_top)
			voigtfit_tot = voigtfit_top + voigtfit_bot
		 	
			ls_trapz_voigt = ls_ratio(integ_voigt_top_trapz,integ_voigt_bot_trapz)
			ls_simps_voigt = ls_ratio(integ_voigt_top_simps, integ_voigt_bot_simps)
				
			if direct==1:
				dirstr = 'r'
			elif direct==-1:
				dirstr='l'
			fold = str(point[0])+"_"+str(point[1]) 
			if fold not in os.listdir("."):
				os.mkdir(fold)
			fileref= open(fold+"/"+str(point[0])+"_"+str(point[1])+"_" + str(num) +"_area_"+str(pixel) +'dir_'+dirstr+".csv" ,'w')
			fileref.write('eV intens trans_l trans_r XMCD combo_l combo_r abs_l abs_r background_l background_r lorentzfit gaussfit voigtfit ls_lor_simps ls_lor_trapz ls_gauss_simps ls_gauss_trapz ls_voigt_simps ls_voigt_trapz '+'\n')
			
			for k in range(len(photon_energies)):
				fileref.write(str(photon_energies[k]) + " " + str(intens[k]) + " " + str(back_l[k]) +" " +str(back_r[k])+" " 
					+str(xmcd[k]) +" "+str(combo_l[k])+" "+str(combo_r[k]) +" " +str(abs_l[k]) + " " +str(abs_r[k])+" "+str(norm_l[k]) + " " +str(norm_r[k]) +" " +str(lorentz_tot[k])+" " +str(gaussfit_tot[k])+" " +str(voigtfit_tot[k])+" " +str(ls_simps_lor)+" " +str(ls_trapz_lor)+" " +str(ls_simps_gauss)+" " +str(ls_trapz_gauss)+" " +str(ls_simps_voigt)+" " +str(ls_trapz_voigt)+ "\n")
			'''
			plt.plot(photon_energies, combo_r,label='absorption')

			plt.close()
			plt.plot(photon_energies, combo_l,label='absorption')
			plt.title('combined left', size=20)
			plt.xlabel('Phot on Energies [eV]',size =20)
			plt.ylabel('Intensity',size=20)
			#plt.legend(loc=2)
			plt.savefig(str(point[0]) + "_" +str(point[1]) + "_" + str(num) +"_area_" +str(pixel)+'combo_l.png')#,dpi=500)
			plt.close()

			plt.plot(photon_energies, abs_r,label='absorption')
			plt.xlabel('Photon Energies [eV]',size =20)
			plt.title('absorption right', size=20)
			plt.ylabel('Intensity',size=20)
			#plt.legend(loc=2)
			plt.savefig(str(point[0]) + "_" +str(point[1]) + "_" + str(num) +"_area_" +str(pixel)+'absorption_r.png')#,dpi=500)
			plt.close()
			plt.plot(photon_energies, abs_l,label='absorption')
			plt.title('absorption left', size=20)
			plt.xlabel('Photon Energies [eV]',size =20)
			plt.ylabel('Intensity',size=20)
			#plt.legend(loc=2)
			plt.savefig(str(point[0]) + "_" +str(point[1]) + "_" + str(num) +"_area_" +str(pixel)+'absorption_l.png')#,dpi=500)
			plt.close()
			'''
			'''
	
			plt.plot(photon_energies, xmcd,'o',label='absorption')
			plt.plot(photon_energies, lorentz_tot,label='absorption')
			plt.title('xmcd', size=20)
			plt.xlabel('Photon Energies [eV]',size =20)
			plt.ylabel('Intensity',size=20)
			#plt.legend(loc=2)
			plt.savefig(str(point[0]) + "_" +str(point[1]) + "_" + str(num) +"_area_" +str(pixel)+'xmcd.png')#,dpi=500)
			plt.close()
			'''
			'''
			plt.plot(photon_energies, back_l,label='absorption')
			plt.title('raw left', size=20)
			plt.xlabel('Photon Energies [eV]',size =20)
			plt.ylabel('Intensity',size=20)
			#plt.legend(loc=2)
			plt.savefig(str(point[0]) + "_" +str(point[1]) + "_" + str(num) +"_area_" +str(pixel)+'raw_l.png')#,dpi=500)
			#if num == 2 and pixel == 8:   plt.show()	
			plt.close()
			plt.plot(photon_energies, back_r,label='absorption')
			plt.title('raw right', size=20)
			plt.xlabel('Photon Energies [eV]',size =20)
			plt.ylabel('Intensity',size=20)
			#plt.legend(loc=2)
			plt.savefig(str(point[0]) + "_" +str(point[1]) + "_" + str(num) +"_area_" +str(pixel)+'raw_r.png')#,dpi=500)
			plt.close()	

			plt.plot(photon_energies, norm_l,label='absorption')
			plt.title('backround left', size=20)
			plt.xlabel('Photon Energies [eV]',size =20)
			plt.ylabel('Intensity',size=20)
			#plt.legend(loc=2)
			plt.savefig(str(point[0]) + "_" +str(point[1]) + "_" + str(num) +"_area_" +str(pixel)+'norm_l.png')#,dpi=500)

			plt.close()
			plt.plot(photon _energies, norm_r,label='absorption')
			plt.title('background right', size=20)
			plt.xlabel('Photon Energies [eV]',size =20)
			plt.ylabel('Intensity',size=20)
			#plt.legend(loc=2)
			plt.savefig(str(point[0]) + "_" +str(point[1]) + "_" + str(num) +"_area_" +str(pixel)+'norm_r.png')#,dpi=500)
			plt.close()	
					'''
			fileref.close()
			del back_l
			del back_r
			del combo_l
			del combo_r
			del abs_r
			del abs_l
			del norm_l
			del norm_r
			''' '''
	
	return totxmcd

#coord=[749, 543]
#producedata(list_images_l, list_images_r, pixel_list, 4, coord)
#right

xmcdr_1 = producedata(list_images_l,list_images_r, pixel_list, 4, [341,620],direct=1,det=0)

xmcdr_2 =producedata(list_images_l,list_images_r, pixel_list, 4, [652,501],direct=-1,det=0)
xmcdr_3 =producedata(list_images_l,list_images_r, pixel_list, 4, [532,518],direct=1,det=0)
xmcdr_4 =producedata(list_images_l,list_images_r, pixel_list, 4, [833,620],direct=1,det=0)
'''
producedata(list_images_l,list_images_r, pixel_list, 4, [334,525],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [334,525],direct=-1,det=0)

producedata(list_images_l,list_images_r, pixel_list, 4, [334,530],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [334,530],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [335,535],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [335,535],direct=-1,det=0)



producedata(list_images_l,list_images_r, pixel_list, 4, [334,540],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [334,540],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [334,545],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [334,545],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [334,550],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [334,550],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [334,555],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [334,555],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [334,560],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [335,560],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [338,565],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [338,565],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [334,570],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [334,570],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [337,575],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [337,575],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [341,590],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [341,590],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [341,595],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [341,595],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [341,600],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [341,600],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [341,605],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [341,605],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [341,610],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [341,610],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [341,615],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [341,615],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [342,620],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [342,620],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [342,625],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [342,625],direct=-1,det=0)
'''



producedata(list_images_l,list_images_r, pixel_list, 4, [535,520],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [513,540],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [523,558],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [514,602],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [414,643],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [458,493],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [465,527],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [424,707],direct=1,det=0)

producedata(list_images_l,list_images_r, pixel_list, 4, [535,520],direct=-1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [513,540],direct=-1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [523,558],direct=-1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [514,602],direct=-1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [414,643],direct=-1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [458,493],direct=-1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [465,527],direct=-1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [424,707],direct=-1,det=0)

#left
xmcdl_1 = producedata(list_images_l,list_images_r, pixel_list, 4, [341,620],direct=-1,det=0)
xmcdl_2=producedata(list_images_l,list_images_r, pixel_list, 4, [652,501],direct=1,det=0)
xmcdl_4 =producedata(list_images_l,list_images_r, pixel_list, 4, [829,620],direct=-1,det=0)

producedata(list_images_l,list_images_r, pixel_list, 4, [534,526],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [534,526],direct=-1,det=0)

producedata(list_images_l,list_images_r, pixel_list, 4, [530,542],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [530,542],direct=-1,det=0)

producedata(list_images_l,list_images_r, pixel_list, 4, [517,582],direct=-1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [517,582],direct=1,det=0)

producedata(list_images_l,list_images_r, pixel_list, 4, [524,634],direct=-1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [524,634],direct=1,det=0)

producedata(list_images_l,list_images_r, pixel_list, 4, [485,609],direct=-1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [485,609],direct=1,det=0)

producedata(list_images_l,list_images_r, pixel_list, 4, [488,580],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [488,580],direct=-1,det=0)

producedata(list_images_l,list_images_r, pixel_list, 4, [474,539],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [474,539],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [455,505],direct=-1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [455,505],direct=1,det=0)

producedata(list_images_l,list_images_r, pixel_list, 4, [726,722],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [726,722],direct=-1,det=0)

producedata(list_images_l,list_images_r, pixel_list, 4, [248,592],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [248,592],direct=-1,det=0)


producedata(list_images_l,list_images_r, pixel_list, 4, [345,591],direct=1,det=0)
producedata(list_images_l,list_images_r, pixel_list, 4, [345,591],direct=-1,det=0)
#ycoords = np.arange(710, 800,5)
#xcoords= np.ones(len(ycoords))*348
#for i in range(len(xcoords)):
	#producedata(list_images_l,list_images_r,pixel_list,15,[752, 537])
#	producedata(list_images_l,list_images_r,pixel_list,15,[xcoords[i], ycoords[i]])
#	print 'point #',i
