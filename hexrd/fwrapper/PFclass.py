import Symclass as sym
import numpy as np
from fortranlib import lambert as lam
from math import atan2, acos
import time

class PF:

	''' this is the pole figure class which deal with
	    any sort of computation or i/o of pole figures '''

	def __init__(self, npoles, nlam, hkl):

		print("Initializing pole figure data...")
		tstart = time.time()

		'''number of pole figures'''
		self.npoles 	 = npoles 				# number of input pole figures

		'''crystallographic pole directions; shape is 3 x npoles'''
		self.hkl 		 = hkl

		'''parameter for lambert equal area projection'''
		self.nlam		 = nlam 					# this is the lambert projection sampling paramter
		self.mm 		 = self.npoles * (2 * self.nlam + 1)**2

		'''the lambert coordinates of the regular grid'''
		lamx 			 = np.array([x / self.nlam for x in np.arange(-self.nlam, self.nlam+1,dtype=np.float64)],dtype=np.float64)
		lamy 			 = np.array([x / self.nlam for x in np.arange(-self.nlam, self.nlam+1,dtype=np.float64)],dtype=np.float64)
		xylist 			 = np.array([[x, y] for x in lamx for y in lamy],dtype=np.float64)

		'''converted to unit directions on the northern hemisphere'''
		ierr 			 = 0
		self.xyzlist 		 = np.array([lam.lambert2dsquareforwarddouble(xy, ierr) for xy in xylist],dtype=np.float64)

		'''the x,y and z components are stored as three different arrays
		no specific reason, just the way things are done in MTEX'''

		'''x-component of direction cosines of sampled points'''
		self.x 			 = self.xyzlist[0:,0]
		'''y-component of direction cosines of sampled points'''
		self.y 			 = self.xyzlist[0:,1]
		'''z-component of direction cosines of sampled points'''
		self.z 			 = self.xyzlist[0:,2]

		'''linear ids of the sampling points. this is used to determine which row
		of the forward matrix is interpolated to'''

		self.idlam 		 = np.array([np.rint(x * self.nlam + self.nlam) * (2 * self.nlam + 1) + \
									np.rint(y * self.nlam + self.nlam) for x in lamx for y in lamy],dtype=np.int32)

		'''intensities at sampled points'''
		self.intensity 	 = np.ones(self.x.shape[0],dtype=np.float64)

		'''spherical angles at the sampled points
		not really required but calculated nonetheless'''

		#azimuthal angle
		self.azimuth 	 = np.array([atan2(y,x) for x in self.x for y in self.y])
		# polar angle
		self.polar 		 = np.array([acos(z) for z in self.z])

		tstop = time.time()
		trun = tstop - tstart
		print ("Execution time [secs] = ", "%8.2f" % trun, "\n")

	def addnoise(self):
		pass

	def WritePFData(self):
		pass

	def ImportPFData(self):
		pass

	def normalize(self):
		pass
