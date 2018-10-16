import Symclass as sym
import numpy as np
import constants as const
import time

class ODF:

	''' this is the oreintation distribution function class
	    contains all functions and variables required to perform
	    calculation as well as i/o of the data.'''

	def __init__(self, space, ncub, crystal_pgnum):

		print("Initializing ODF data...")
		tstart = time.time()

		''' since we are using a discrete method of representation, we
            represent the data at discrete grid popints in the case of
            cubochoric grid and a discrete nodal points in the case of
            the rodrigues parametrization.'''

		if space.lower() == 'cubochoric':
			print("Setting space to Cubochoric coordinates")
		elif space.lower() == 'rodrigues':
			print ("Setting space to Rodrigues parametrization")
		else:
			raise Exception("only Cubochoric and Rodrigues representation of SO(3) is allowed")

		self.pgnum = crystal_pgnum

		self.space = space  	#'cubochoric' or 'rodrigues'

		''' these variables are for the cubochoric representation only
		 number of points in semi-edge of cube; determines the resolution of the ODF'''
		self.ncub     = ncub
		self.nn 	  = (2 * self.ncub + 1)**3
		# integer value of the grid points at these coordinates
		self.gridpts  = np.array([[x,y,z] for x in np.arange(-self.ncub,self.ncub+1,dtype=np.float64) \
								for y in np.arange(-self.ncub,self.ncub+1,dtype=np.float64) \
								for z in np.arange(-self.ncub, self.ncub + 1, dtype=np.float64)], \
								dtype=np.float64)
		c = const.Constants()

		# cubochoric coordinates stored as numpy array
		self.cugrid  = np.array([[x/c.Lambert3D.ap,y/c.Lambert3D.ap,z/c.Lambert3D.ap] \
								for x in np.arange(-self.ncub,self.ncub+1,dtype=np.float64) \
								for y in np.arange(-self.ncub,self.ncub+1,dtype=np.float64) \
								for z in np.arange(-self.ncub, self.ncub + 1, dtype=np.float64)], \
								dtype=np.float64)


		''' these variables are for the rodrigues representation only
		to be appended as the implementation details become clearer'''
		nodes   = [] # node points in rodrigues space saved as list

		''' symmetry initialization for both cryustal and sample '''
		self.crystal_symmetry = sym.Symmetry(self.pgnum) # defaults to cubic crystal symmetry
		self.sample_symmetry  = sym.Symmetry(1)  # defaults to triclinic crystal symmetry

		tstop = time.time()
		trun = tstop - tstart
		print("Execution time [secs] = ", "%8.2f" % trun, "\n")

	def normalize(Self):
		pass

	def TextureIndex(self):
		pass

	def ImportODF(self):
		pass

	def WriteODF(self):
		pass