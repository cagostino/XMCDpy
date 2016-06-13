#lets gooooo. Authors: Christopher Agostino & MacCallum Robertson
import sys
import numpy as np
import os
import scipy
import gdal
import matplotlib.pyplot as plt
import Image


from PySide import QtCore, QtGui
from PySide.QtGui import *
from PySide.QtCore import QTimer, SIGNAL, SLOT
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
from mpl_toolkits.mplot3d import Axes3D

from lmfit.models import LorentzianModel
from lmfit.models import GaussianModel
from lmfit.models import VoigtModel
#import pyfits
import gc
from numpy import trapz
from scipy.integrate import simps
plt.rc('text',usetex=True)
#make these plots look dank
plt.rc('font',family='serif')
l=[]
r=[]
l = os.listdir("left/")
r =os.listdir("right/")
l.sort()
r.sort()

#GUI Part
app =QApplication.instance()

if not app:
    app = QApplication(sys.argv)

   
      
#################First Window####################   
w = QWidget()

#Window#
w.resize(650,625)
w.setWindowTitle('Local XMCD Spectro-Microscopy (LoXS)')




def show1(w,energy,energy2):
    
    #Image#
    imagelabel=QLabel(w)
    image=QImage('Capture2.PNG')
    imagelabel.setPixmap(QPixmap.fromImage(image))
    #image.move(300,100)
    scrollArea=QScrollArea(w)
    scrollArea.setBackgroundRole(QPalette.Dark)
    scrollArea.setWidget(imagelabel)
    scrollArea.move(350,10)
    scrollArea.resize(280,230)
    scrollArea.setVisible(False)
        
    #Image#
    imagelabel2=QLabel(w)
    image2=QImage('capture3.PNG')
    imagelabel2.setPixmap(QPixmap.fromImage(image2))
    #image.move(300,100)
    scrollArea2=QScrollArea(w)
    scrollArea2.setBackgroundRole(QPalette.Dark)
    scrollArea2.setWidget(imagelabel2)
    scrollArea2.move(350,270)
    scrollArea2.resize(280,230)
    scrollArea2.setVisible(False)

    scrollArea.setVisible(True)
    scrollArea2.setVisible(True)
    


def show2(w,energy,energy2):
    
    #Image#
    imagelabel3=QLabel(w)
    image3=QImage('capture.PNG')
    imagelabel3.setPixmap(QPixmap.fromImage(image3))

    scrollArea3=QScrollArea(w)
    scrollArea3.setBackgroundRole(QPalette.Dark)
    scrollArea3.setWidget(imagelabel3)
    scrollArea3.move(350,10)
    scrollArea3.resize(280,230)
    scrollArea3.setVisible(False)
        
    #Image#
    imagelabel4=QLabel(w)
    image4=QImage('capture.PNG')
    imagelabel4.setPixmap(QPixmap.fromImage(image4))
    
    scrollArea4=QScrollArea(w)
    scrollArea4.setBackgroundRole(QPalette.Dark)
    scrollArea4.setWidget(imagelabel4)
    scrollArea4.move(350,270)
    scrollArea4.resize(280,230)
    scrollArea4.setVisible(False)

    scrollArea3.setVisible(True)
    scrollArea4.setVisible(True)
    
#Drop Down Material Menu#
material=QComboBox(w)
material.resize(150,20)

material.addItem('Metal')
material.addItem('Cobalt')
material.addItem('Nickel')
material.addItem('Iron')
material.move(30,35)



#Buttons#
button2= QPushButton('Show',w)
button2.setToolTip('Show L3 Images')
button2.resize(button2.sizeHint())
button2.move(375,530)
 
button3= QPushButton('Show 2',w)
button3.setToolTip('Show L3 Images')
button3.resize(button3.sizeHint())
button3.move(500,530)


button4= QPushButton('Run',w)
button4.setToolTip('Run the data collection script')
button4.resize(button4.sizeHint())
button4.move(500,575)

button5= QPushButton('Skip',w)
button5.setToolTip('Skip to data analysis section ')
button5.resize(button5.sizeHint())
button5.move(375,575)


#Labels#
label1 = QLabel(w)
label1.setText('Local XMCD Spectro-Microscopy (LoXS) ')
label1.move(30, 10)

#mb=QMessageBox(w)
#mb.setText("hi")
#mb.exec_()

label2 = QLabel(w)
label2.setText('Select Coordinates (Center of DW)')
label2.move(30, 115)

label3 = QLabel(w)
label3.setText('Select Area Size (Pixels)')
label3.move(30, 185)

label4 = QLabel(w)
label4.setText('                                                ')
label4.move(350, 245)

label5 = QLabel(w)
label5.setText('                                                 ')
label5.move(350,505)

label3 = QLabel(w)
label3.setText('Select Number of Areas')
label3.move(30, 245)

#Textboxes#
textbox1 = QLineEdit(w)
textbox1.move(30, 135)
textbox1.resize(200, 20)
textbox1.setText('Enter X')

textbox2 = QLineEdit(w)
textbox2.move(30, 160)
textbox2.resize(200, 20)
textbox2.setText('Enter Y')

textbox3 = QLineEdit(w)
textbox3.move(185, 35)
textbox3.resize(50, 20)
textbox3.setText('L3 (eV)')

textbox4 = QLineEdit(w)
textbox4.move(240, 35)
textbox4.resize(50, 20)
textbox4.setText('L2 (eV)')

textbox5 = QLineEdit('Import Polarization 1 (.tif)', w)
textbox5.setText('File Name For First Polarization')
textbox5.resize(200,20)
textbox5.move(30, 60)

textbox6= QLineEdit(w)
textbox6.setText('File Name For Second Polarization')
textbox6.resize(200,20)
textbox6.move(30,90)

textbox7= QLineEdit(w)
textbox7.setText('')
textbox7.resize(30,20)
textbox7.move(30,265)

store_count= QLineEdit(w)
store_count.setText('')
store_count.resize(30,20)
store_count.move(30,265)
store_count.hide()


#Check Boxes
Check1=QCheckBox("Transmission intensity",w)
Check1.move(30,290)
Check1.setChecked(False)

Check2=QCheckBox("Absorption intensity",w)
Check2.move(30,320)
Check2.setChecked(False)

Check3=QCheckBox("Local Background",w)
Check3.move(30,350)
Check3.setChecked(False)

Check4=QCheckBox("Background Subtracted",w)
Check4.move(30,380)
Check4.setChecked(False)

Check5=QCheckBox("XMCD",w)
Check5.move(30,410)
Check5.setChecked(False)

Check6=QCheckBox("Voigt Fit",w)
Check6.move(30,440)
Check6.setChecked(False)

Check7=QCheckBox("Lorentzian Fit",w)
Check7.move(30,470)
Check7.setChecked(False)

Check8=QCheckBox("Gaussian Fit",w)
Check8.move(30,500)
Check8.setChecked(False)

Check9=QCheckBox("L/S Grid",w)
Check9.move(30,530)
Check9.setChecked(False)

Check10=QCheckBox("L/S vs. Area",w)
Check10.move(30,560)
Check10.setChecked(False)

Check11=QCheckBox("5 x 5 ",w)
Check11.move(30,200)
Check11.setChecked(False)

Check12=QCheckBox("8 x 8",w)
Check12.move(100,200)
Check12.setChecked(False)

Check13=QCheckBox("10 x 10",w)
Check13.move(30,220)
Check13.setChecked(False)

Check14=QCheckBox("15 x 15",w)
Check14.move(100,220)
Check14.setChecked(False)

image=QImage()
image.isNull()
#image.move(300,300)
#image.load('filename')


#Actions
def on_click_button2():
    #c1=closest(float(textbox3.text()),photon_energies)
    #c2=closest(float(textbox4.text()),photon_energies)
    #c1=str(c1)
    #c2=str(c2)
    #if '.' in c1:
    #    c1 = c1.replace(".","_")
    #if '.' in c2:
    #    c2 = c2.replace(".","_")
    #c1=c1+'_R.tif' 
    #c2=c2+'_R.tif'   
    
    show1(w,3,3)
    
    
def on_click_button3():
    #c3=closest(float(textbox3.text()),photon_energies)
    #c4=closest(float(textbox4.text()),photon_energies)

    #c3=str(c3)
    #c4=str(c4)
    #if '.' in c3:
    #    c3 = c3.replace(".","_")
    #if '.' in c4:
    #    c4 = c4.replace(".","_")
    #c3=c3+'_L.tif' 
    #c4=c4+'_L.tif'
    #show2(w,c3,c4)
    show2(w,3,3)
    
    
def on_click_button4():
    xcoord=float(textbox1.text())
    ycoord=float(textbox2.text())
    area_num=int(textbox7.text())
    #check()
    #count=0
    #if Check1.isChecked() == True:
    #    count+=1
    #if Check2.isChecked() == True:
    #    count+=1
    #if Check3.isChecked() == True:
    #    count+=1
    #if Check4.isChecked() == True:
    #    count+=1
    #if Check5.isChecked() == True:
    #    count+=1
    #if Check6.isChecked() == True:
    #    count+=1
    #if Check7.isChecked() == True:
    #    count+=1
    #if Check8.isChecked() == True:
    #    count+=1
    #if Check9.isChecked() == True:
    #    count+=1
    #if Check10.isChecked() == True:
    #    count+=1
    #store_count.setText(count)
    #x,y=600,600
    #if count<=2:
    #    x,y=800,500
    #if count>=3 and count<=4:
    #    x,y=1000,500
    #if count>=5 and count<=6:
    #    x,y=1200,500
    #if count>=7 and count<=8:
    #    x,y=1400,500
    #if count>=9 and count<=10:
    #    x,y=1600,500

    #w2.resize(x,y)
    
    #store_count.setText(str(x))
    
    producedata(list_images_l,list_images_r,area_num,[xcoord,ycoord],direct=1,det=0)
    #count=0
 
    w.hide()
    w2.show()
    
def on_click_button5():
    w.hide()
    w2.show()

#ComboBox activation
def on_activated2(text):
	text_str = str(text)
	textbox3.setText(str(Metals[text_str][0]))
	textbox4.setText(str(Metals[text_str][1]))
	label4.setText('First Polarization @' + str(Metals[text_str][0])+ 'eV')
	label5.setText('Second Polarization @' +str(Metals[text_str][1])+ 'eV')


#Metal L edges#
Metals = {
                  'Cobalt':[778, 792],	
                  'Iron':[708, 720],
                  'Nickel':[820, 853]	
                 
}

material.activated[str].connect(on_activated2)

button2.clicked.connect(on_click_button2)
button3.clicked.connect(on_click_button3)
button4.clicked.connect(on_click_button4)
button5.clicked.connect(on_click_button5)
##############################################################
##########################Second Window#######################
w2 = QWidget()

w2.resize(300,500)
w2.setWindowTitle('LoXS Data Viewer')

new='on'


textboxS = QLineEdit(w2)
textboxS.move(30, 185)
textboxS.resize(100, 20)
textboxS.setText(new)
textboxS.hide()

#def show_recent
#def show_new
def showB():
    labelB.show()
    labelB1.show()
    labelB2.show()
    labelB3.show()
    textboxB1.show()
    textboxB2.show()
    
def hideB():
    labelB.hide()
    labelB1.hide()
    labelB2.hide()
    labelB3.hide()
    textboxB1.hide()
    textboxB2.hide()

labelB = QLabel(w2)
labelB.setText('Recent File Parameters')
labelB.move(30, 135)
labelB.show() 
       
labelB1 = QLabel(w2)
labelB1.setText('Input Area Number')
labelB1.move(30, 160)
labelB1.show()
    
labelB2 = QLabel(w2)
labelB2.setText('Input Pixel Value')
labelB2.move(30, 220)
labelB2.show()

labelB3 = QLabel(w2)
labelB3.setText('Selected File:')
labelB3.move(30, 385)
labelB3.show()

    
textboxB1= QLineEdit(w2)
textboxB1.setText('')
textboxB1.resize(80,20)
textboxB1.move(30,190)
textboxB1.show()
    
textboxB2= QLineEdit(w2)
textboxB2.setText('')
textboxB2.resize(80,20)
textboxB2.move(30,250)
textboxB2.show()
           
def showC():
    labelC.show()
    labelC1.show()
    labelC2.show()
    labelC3.show()
    labelC4.show()
    textboxC1.show()
    textboxC2.show()
    textboxC3.show()
    textboxC4.show()
def hideC():
    labelC.hide()
    labelC1.hide()
    labelC2.hide()
    labelC3.hide()
    labelC4.hide()
    textboxC1.hide()
    textboxC2.hide()
    textboxC3.hide()
    textboxC4.hide()
    
labelC = QLabel(w2)
labelC.setText('New File Parameters')
labelC.move(30, 135)
labelC.hide()

labelC1 = QLabel(w2)
labelC1.setText('Input Area Number')
labelC1.move(30, 240)
labelC1.hide()
    
labelC2 = QLabel(w2)
labelC2.setText('Input Pixel Value')
labelC2.move(30, 290)
labelC2.hide()
    
labelC3 = QLabel(w2)
labelC3.setText('Input (X,Y) Coordinates')
labelC3.move(30, 165)
labelC3.hide()
    
labelC4 = QLabel(w2)
labelC4.setText('Selected File:')
labelC4.move(30, 385)
labelC4.hide()
    
    
textboxC1= QLineEdit(w2)
textboxC1.setText('')
textboxC1.resize(80,20)
textboxC1.move(30,260)
textboxC1.hide()
    
textboxC2= QLineEdit(w2)
textboxC2.setText('')
textboxC2.resize(80,20)
textboxC2.move(30,310)
textboxC2.hide()

textboxC3 = QLineEdit(w2)
textboxC3.move(30, 185)
textboxC3.resize(100, 20)
textboxC3.setText('Enter X')
textboxC3.hide()
    
    
    
textboxC4 = QLineEdit(w2)
textboxC4.move(30, 215)
textboxC4.resize(100, 20)
textboxC4.setText('Enter Y')
textboxC4.hide()
  
labelD = QLabel(w2)
labelD.setText('                                                                                            ')
labelD.move(30, 400)  

def plots(x,y,x_axis,y_axis,title,name):
    plt.plot(x, y)
    plt.xlabel(x_axis,fontsize=20)
    plt.ylabel(y_axis,fontsize=20)
    plt.title(title,fontsize=22)
    plt.grid(True)
    plt.savefig(name+".png")
    plt.close()
def plots2(x,y,y2,x_axis,y_axis,info1,info2,title,name):
    left=plt.plot(x, y)
    right=plt.plot(x,y2)
    plt.xlabel(x_axis,fontsize=20)
    plt.ylabel(y_axis,fontsize=20)
    plt.title(title,fontsize=22)
    plt.grid(True)
    plt.savefig(name+".png")
    plt.close()
    



store_value=QLineEdit(w2)
store_value.setText('')
store_value.hide()

store_value2=QLineEdit(w2)
store_value2.setText('')
store_value2.hide()
def makefile():
    
    if textboxS.text()=='on':
        xcoord=textbox1.text()
        ycoord=textbox2.text()
        xcoord=float(xcoord)
        ycoord=float(ycoord)
        fold=str(xcoord)+"_"+str(ycoord)
        filename=2
        #directions = {1:'r',-1:'l'}
        #dirstr = directions[direct]
        dirstr=1
        pixel=textboxB2.text()
        num=textboxB1.text()
        filename=fold+"/"+str(xcoord)+"_"+str(ycoord)+"_" + str(num) +"_area_"+str(pixel) +'dir_'+dirstr+".csv"
    else:
        xcoord=textboxC3.text()
        ycoord=textboxC4.text()
        xcoord=float(xcoord)
        ycoord=float(ycoord)
        fold=str(xcoord)+"_"+str(ycoord)
        dirstr='r'
        pixel=textboxC2.text()
        num=textboxC1.text()
        filename=fold+"/"+str(xcoord)+"_"+str(ycoord)+"_" + str(num) +"_area_"+str(pixel) +'dir_'+dirstr+".csv"    
    
    
    labelD.setText(filename)
    labelD.show()
    
    store_value.setText(filename)

def show_im_A():
    w2.resize(1100,600)
    imageA=QImage('transmission.PNG')
    pixmapA=QPixmap.fromImage(imageA)
    pixmapA=pixmapA.scaled(700,700,QtCore.Qt.KeepAspectRatio)
    imagelabelA=QLabel(w2)
    imagelabelA.setPixmap(pixmapA)
    imagelabelA.move(300,30)
    #imagelabelA.resize(100,100)
    imagelabelA.setScaledContents(True)
    imagelabelA.show()
def show_im_B():
    w2.resize(1100,600)
    imageB=QImage('abs.PNG')
    pixmapB=QPixmap.fromImage(imageB)
    pixmapB=pixmapB.scaled(700,700,QtCore.Qt.KeepAspectRatio)
    imagelabelB=QLabel(w2)
    imagelabelB.setPixmap(pixmapB)
    imagelabelB.move(300,30)
    #imagelabelA.resize(100,100)
    imagelabelB.setScaledContents(True)
    imagelabelB.show()
def show_im_C():
    w2.resize(1100,600)
    imageC=QImage('background.PNG')
    pixmapC=QPixmap.fromImage(imageC)
    pixmapC=pixmapC.scaled(700,700,QtCore.Qt.KeepAspectRatio)
    imagelabelC=QLabel(w2)
    imagelabelC.setPixmap(pixmapC)
    imagelabelC.move(300,30)
    #imagelabelA.resize(100,100)
    imagelabelC.setScaledContents(True)
    imagelabelC.show()
def show_im_D():
    w2.resize(1100,600)
    imageD=QImage('bs.PNG')
    pixmapD=QPixmap.fromImage(imageD)
    pixmapD=pixmapD.scaled(700,700,QtCore.Qt.KeepAspectRatio)
    imagelabelD=QLabel(w2)
    imagelabelD.setPixmap(pixmapD)
    imagelabelD.move(300,30)
    #imagelabelA.resize(100,100)
    imagelabelD.setScaledContents(True)
    imagelabelD.show()
def show_im_E():
    w2.resize(1100,600)
    imageE=QImage('xmcd.PNG')
    pixmapE=QPixmap.fromImage(imageE)
    pixmapE=pixmapE.scaled(700,700, QtCore.Qt.KeepAspectRatio)
    imagelabelE=QLabel(w2)
    imagelabelE.setPixmap(pixmapE)
    imagelabelE.move(300,30)
    #imagelabelA.resize(100,100)
    imagelabelE.setScaledContents(True)
    imagelabelE.show()
def show_im_F():
    w2.resize(1100,600)
    imageF=QImage('voigt.PNG')
    pixmapF=QPixmap.fromImage(imageF)
    pixmapF=pixmapF.scaled(700,700,QtCore.Qt.KeepAspectRatio)
    imagelabelF=QLabel(w2)
    imagelabelF.setPixmap(pixmapF)
    imagelabelF.move(300,30)
    #imagelabelA.resize(100,100)
    imagelabelF.setScaledContents(True)
    imagelabelF.show()
def show_im_G():
    w2.resize(1100,600)
    imageG=QImage('lorentzian.PNG')
    pixmapG=QPixmap.fromImage(imageG)
    pixmapG=pixmapG.scaled(700,700,QtCore.Qt.KeepAspectRatio)
    imagelabelG=QLabel(w2)
    imagelabelG.setPixmap(pixmapG)
    imagelabelG.move(300,30)
    #imagelabelA.resize(100,100)
    imagelabelG.setScaledContents(True)
    imagelabelG.show()
def show_im_H():
    w2.resize(1100,600)
    imageH=QImage('gaussian.PNG')
    pixmapH=QPixmap.fromImage(imageH)
    pixmapH=pixmapH.scaled(700,700,QtCore.Qt.KeepAspectRatio)
    imagelabelH=QLabel(w2)
    imagelabelH.setPixmap(pixmapH)
    imagelabelH.move(300,30)
    #imagelabelA.resize(100,100)
    imagelabelH.setScaledContents(True)
    imagelabelH.show()
     
        
#    if Check9.isChecked() == True:
#        imagelabel3=QLabel(w)
#        image3=QImage('capture.PNG')
#        imagelabel3.setPixmap(QPixmap.fromImage(image3))
        
        
#    if Check10.isChecked() == True:
#        imagelabel3=QLabel(w)
#        image3=QImage('capture.PNG')
#        imagelabel3.setPixmap(QPixmap.fromImage(image3))

#def on_click_buttonA():
    #print 10
    #show_plots()
    
    
    
    


def on_click_buttonB():
    
    hideC()
    showB()
    textboxS.setText('on')
    #count=textboxS.text()
    #count=float(count)
    #if count%2==1:
    #    showC()
    #if count%2==0:
    #    hideC()
    #print count
    
    #count=count+1
    #textboxS.setText(str(count))
    
def on_click_buttonC():
    hideB()
    showC()
    textboxS.setText('off')
    #count=textboxS.text()
    #count=float(count)
    #if count%2==1:
    #    showC()
    #if count%2==0:
    #    hideC()
    #print count
    
    #count=count+1
    
def on_click_buttonD():
    makefile()
    filename=store_value.text()
    read_dat(filename)
    #show_plots()
    
        

#buttonA= QPushButton('Push To Plot',w2)
#buttonA.setToolTip('Display Selected Plots')
#buttonA.resize(buttonA.sizeHint())
#buttonA.move(30,430)

buttonB= QPushButton('Use Recent Files',w2)
buttonB.setToolTip('Select files just generated by LoXS')
buttonB.resize(buttonB.sizeHint())
buttonB.move(30,50)

buttonC= QPushButton('Use New Files',w2)
buttonC.setToolTip('Select already existing files ')
buttonC.resize(buttonC.sizeHint())
buttonC.move(30,90)

buttonD= QPushButton('Retrieve File',w2)
buttonD.setToolTip('Retrieve file corresponding to the parameters above')
buttonD.resize(buttonD.sizeHint())
buttonD.move(30,340)

options=QComboBox(w2)
options.resize(150,20)
options.move(30,430)
options.hide()

def on_activated3(text):
	text_str = str(text)
	if text_str=='Transmission':
	    show_im_A()
	if text_str=='Absorption':
	    show_im_B()
	if text_str=='Background':
	    show_im_C() 
	if text_str=='Background Subtraction':
	    show_im_D()
	if text_str=='XMCD':
	    show_im_E()   
	if text_str=='Voigt':
	    show_im_F()    
	if text_str=='Lorentzian':
	    show_im_G()    
	if text_str=='Gaussian':
	    show_im_H()   
	     


options.activated[str].connect(on_activated3)
#buttonA.clicked.connect(on_click_buttonA)
buttonB.clicked.connect(on_click_buttonB)
buttonC.clicked.connect(on_click_buttonC)
buttonD.clicked.connect(on_click_buttonD)
############################################
def closest(num,lst):
    closest_num=min(photon_energies, key=lambda x:abs(x-num))
    return closest_num

def read_dat(filename):
    fil = np.loadtxt(filename, skiprows=1)
    fil_t=fil.transpose()
    
    ##Transmission##
    if Check1.isChecked() == True:
        plots2(fil_t[0],fil_t[2],fil_t[3],"Photon Energy", "Intensity","Polarization 1","Polarization 2", "Transmitted Intensity", "transmission")
        options.addItem('Transmission')
    ##Absorption##
    if Check2.isChecked() == True:
        plots2(fil_t[0],fil_t[7],fil_t[8],"Photon Energy", "Intensity","Polarization 1","Polarization 2", "Absorption Intensity", "abs")
        options.addItem('Absorption')
    ##background##
    if Check3.isChecked() == True:
        plots2(fil_t[0],fil_t[9],fil_t[10],"Photon Energy", "Intensity","Polarization 1","Polarization 2", "Background Intensity", "background")
        options.addItem('Background')
    ##combo##
    if Check4.isChecked() == True:
        plots2(fil_t[0],fil_t[5],fil_t[6],"Photon Energy", "Intensity","Polarization 1","Polarization 2","Background Subtracted", "bs")
        options.addItem('Background Subtraction')
    ##XMCD##
    if Check5.isChecked() == True:
        plots(fil_t[0],fil_t[4],"Photon Energy", "Intensity", "XMCD", "xmcd")
        options.addItem('XMCD')
    ##Voigt##
    if Check6.isChecked()==True:
        plots2(fil_t[0],fil_t[4],fil_t[13],"Photon Energy", "Intensity","XMCD","Voigt", "Voigt Fit", "voigt")
        options.addItem('Voigt')
    ##Lorenztian##
    if Check7.isChecked()==True:
        plots2(fil_t[0],fil_t[4],fil_t[11],"Photon Energy", "Intensity","XMCD","Lorentzian", "Lorentzian Fit", "lorentzian")
        ptions.addItem('Lorentzian')   
    ##Gaussian##
    if Check8.isChecked()==True:
        plots2(fil_t[0],fil_t[4],fil_t[12],"Photon Energy", "Intensity", "XMCD","Gaussian","Gaussian Fit", "gaussian")
        ptions.addItem('Gaussian')
    options.show()


	
	
def assign_images(lst = None):
	new_array_r =[]
	file_list_r =[]
	new_array_l =[]
	file_list_l =[]
	for fil in l:
		if fil.endswith(".tif"):
			new =gdal.Open('left/'+fil)
			new_array_l.append(scipy.array(new.GetRasterBand(1).ReadAsArray()))
			file_list_l.append(fil)
	for fil in r:
		if fil.endswith(".tif"):
			new =gdal.Open('right/'+fil)
			new_array_r.append(scipy.array(new.GetRasterBand(1).ReadAsArray()))
			file_list_r.append(fil)	
	return [file_list_l,file_list_r, new_array_l,new_array_r]
list_files_l,list_files_r , list_images_l, list_images_r = assign_images()
photon_energies = []
for i in list_files_l:
	if 'tif' in i:
		if '_R' in i:
			i = i.replace("_R.tif","")
		if "_L" in i:
			i=i.replace("_L.tif",'')
                if "_" in i:
			i=i.replace("_",'.')	
		else:
			i=i.replace(".tif",'')
		

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
	#determines best normalization given a certain point, not usually necessary
	diff_1 = np.array([])
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
	#now for actual xmcd calculation
	norm_l =[]
	norm_r =[]
	#gets the left and right normalizations
	for j in range(-size,size+1): 
		for i in range(size):
			norm_l.append(image_left[y+i][x+j])
			norm_r.append(image_right[y+i][x+j])
	#turn them into numpy arrays idk why
	norm_l, norm_r = np.array(norm_l),np.array(norm_r) 
	#set them equal to their mean to subtract from the XMCD later on
	norm_l =np.mean(norm_l)
	norm_r = np.mean(norm_r)
	
	back = (norm_l+norm_r)
	#overall background
	#now get the actual xray absorption curves for each
	# does the right polarized images
	while count_x < size:
		intensities_right.append(image_right[y][x+direct*count_x+size*(number-1)])
		intensities_left.append(image_left[y][x+direct*count_x+direct*size*(number-1)])
		count_y = 1
		while count_y < size:
			intensities_right.append(image_right[y+count_y][x+direct*count_x+direct*size*(number-1)])
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
	#performs a single peak lorentzian fit on data
	mod = LorentzianModel()
	pars = mod.guess(xmcd, x=photon_e)
	out  = mod.fit(xmcd, pars, x=photon_e)
	return out.best_fit
def gauss(photon_e, xmcd):
	#performs single peak gaussian fit on data
	mod= GaussianModel()
	pars = mod.guess(xmcd, x=photon_e)
	out=mod.fit(xmcd, pars,x=photon_e)
	return out.best_fit
def voigt(photon_e, xmcd):
	#performs single peak voigt fit on data
	mod =VoigtModel()
	pars = mod.guess(xmcd, x= photon_e)
	out=mod.fit(xmcd, pars, x=photon_e)
	return out.best_fit
def ls_ratio(a3,a2):
	#using sum rules, calculates the L/S ratio
	ls=(2*(a3+a2))/(3*(a3-2*a2))
	return ls


#def check():
#    pixel_list=[]
#    if Check11.isChecked() == True:
#        pixel_list.append(5)
    
#    if Check12.isChecked() == True:
#        pixel_list.append(8)
    
#    if Check13.isChecked() == True:
#        pixel_list.append(10)
    
#    if Check14.isChecked() == True: 
#        pixel_list.append(15)  
#    print pixel_list 

   
def producedata(image_list_l,image_list_r,number,point,direct=1,det=0):
	#produces lots of xmcd curves for a given point using different sizes for averaging and different distances 
	pixel_list=[]
        if Check11.isChecked() == True:
            pixel_list.append(5)
    
        if Check12.isChecked() == True:
            pixel_list.append(8)
    
        if Check13.isChecked() == True:
            pixel_list.append(10)
    
        if Check14.isChecked() == True: 
            pixel_list.append(15)  
	
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
			#assign all useful info from calculations
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
			#SUBTRACTING  A LINEAR FIT OF THE START AND END PART BECAUSE THEY ARE SUPPOSED TO BE ZERO
			tes,xtes = np.array(xmcd[0:30]), np.array(photon_energies[0:30])
			tes,xtes= np.append(tes,xmcd[-10:]), np.append(xtes, photon_energies[-10:])
			m, yint = np.polyfit(xtes, tes,1)
			#ZEROED XMCD
			xmcd = xmcd- np.array(photon_energies) *m-yint #trying to do linear fit
	
			
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
			directions = {1:'r',-1:'l'}
			dirstr = directions[direct]
			fold = str(point[0])+"_"+str(point[1]) 
			if fold not in os.listdir("."):
				os.mkdir(fold)
			fileref= open(fold+"/"+str(point[0])+"_"+str(point[1])+"_" + str(num) +"_area_"+str(pixel) +'dir_'+dirstr+".csv" ,'w')
			fileref.write('eV intens trans_l trans_r XMCD combo_l combo_r abs_l abs_r background_l background_r lorentzfit gaussfit voigtfit ls_lor_simps ls_lor_trapz ls_gauss_simps ls_gauss_trapz ls_voigt_simps ls_voigt_trapz '+'\n')
			
			for k in range(len(photon_energies)):
				fileref.write(str(photon_energies[k]) + " " + str(intens[k]) + " " + str(back_l[k]) +" " +str(back_r[k])+" " 
					+str(xmcd[k]) +" "+str(combo_l[k])+" "+str(combo_r[k]) +" " +str(abs_l[k]) + " " +str(abs_r[k])+" "+str(norm_l[k]) + " " +str(norm_r[k]) +" " +str(lorentz_tot[k])+" " +str(gaussfit_tot[k])+" " +str(voigtfit_tot[k])+" " +str(ls_simps_lor)+" " +str(ls_trapz_lor)+" " +str(ls_simps_gauss)+" " +str(ls_trapz_gauss)+" " +str(ls_simps_voigt)+" " +str(ls_trapz_voigt)+ "\n")
			fileref.close()
			del back_l
			del back_r
			del combo_l
			del combo_r
			del abs_r
			del abs_l
			del norm_l
			del norm_r
	print totxmcd
	
w.show()
#w.hide()
w2.hide()

sys.exit(app.exec_())
	
#producedata(list_images_l,list_images_r, pixel_list,5,[500,500],direct=1,det=0)