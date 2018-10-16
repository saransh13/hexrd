from numpy import pi as cPi

class Constants:
    def __init__(self):
        self.Lambert2D = self.Lambert2D()
        self.Lambert3D = self.Lambert3D()

    class Lambert2D:
        def __init__(self):
            self.ap = 1.0

    class Lambert3D:
        def __init__(self):
            self.ap = cPi**(2.0/3.0)



