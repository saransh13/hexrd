import time
from fortranlib import crystal

class Crystal:

    ''' this class deals with all the matrices needed for the crystallographic
    computations such as structure matrix and metric tensors in direct and
    reciprocal space. As with previous routines, we will use the fortran code
    and simply write the simplest possible wrapper to call the already tried
    and tested routines'''

    def __init__(self):
        [self.dmt, self.rmt, self.dsm, self.rsm, self.trigmat]
