from fortranlib import symmetry as sym
from fortranlib import so3 as so3

class Symmetry:
	''' this is the symmetry class and handels all the "ROTATIONAL"
	    symmetries for any point group. Note that general space group
	    symmetry handling will be added in the future '''

	def __init__(self, pgnum):

		[FZtype, FZorder] = so3.getfztypeandorder(pgnum)

		[Nqsym, Pm] = sym.symmetry_init(pgnum)

		self.pgnum 		= pgnum
		self.FZtype 	= FZtype
		self.FZorder 	= FZorder
		self.Nqsym 		= Nqsym
		self.Pm 		= Pm