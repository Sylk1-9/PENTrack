import magpylib as magpy
import numpy as np
from scipy.spatial.transform import Rotation as R

print("Script visited : SF2.py in in/")


## funktion to generate piece of circle as numpy array
class SF_field:
    def __init__(self, phi, phi_off , x, x_off, r, cur, k):
        self.phi = phi   # open angle of saddle coil in degree
        self.phi_off = phi_off # angular offset of coil from center
        self.x = x  # length of saddle coil in mm
        self.x_off = x_off
        self.r = r  # radius of saddle coil in mm
        self.cur = cur  # current through saddle coil in A
        self.k = k  # resolution of saddle coil curves

    def circ(self, phi, x, x_off, r):  # generate points on a circle of radius r at height z
        return np.array([[r * np.sin(phii / 360 * 2 * np.pi), r * np.cos(phii / 360 * 2 * np.pi), x + x_off] for phii in phi])

    def sf_field(self):

        ## create arrays copntaining all the angles for wires around the cylinder
        angles = np.linspace(-self.phi/2 + self.phi_off, self.phi/2 + self.phi_off, self.k)
        angles2 = np.linspace(self.phi/2 + self.phi_off, -self.phi/2 + self.phi_off, self.k)
        angles3 = np.linspace(self.phi/2 + 180 + self.phi_off, -self.phi/2 + 180 + self.phi_off, self.k)
        angles4 = np.linspace(-self.phi/2 + 180 + self.phi_off, self.phi/2 + 180 + self.phi_off, self.k)

        # round sadlle coil paths
        rpath1 = self.circ(angles2, -self.x/2, self.x_off, self.r)
        rpath2 = self.circ(angles, self.x/2, self.x_off, self.r )
        rpath3 = self.circ(angles4, -self.x/2, self.x_off, self.r)
        rpath4 = self.circ(angles3, self.x/2, self.x_off, self.r)
        # start and end positions of coil rungs

        startrung1 = rpath1[-1]
        endrung1 = rpath2[0]
        startrung2 = rpath2 [-1]
        endrung2 = rpath1[0]
        endrung3 = rpath4[0]
        startrung3 = rpath3[-1]
        endrung4 = rpath3[0]
        startrung4 = rpath4[-1]


        sf_objects = [
            mag.current.Line(self.cur, [startrung1, endrung1]),
            mag.current.Line(self.cur, rpath1),
            mag.current.Line(self.cur, rpath2),
            mag.current.Line(self.cur, [startrung2, endrung2]),
            mag.current.Line(self.cur, rpath3),
            mag.current.Line(self.cur, rpath4),
            mag.current.Line(self.cur, [ startrung3, endrung3]),
            mag.current.Line(self.cur, [startrung4, endrung4])
        ]
        spinflipper = mag.Collection(sf_objects, style_label="sf1_c1")

        return spinflipper


def buildSource(t=0, winding=100):
    '''
    The axis is NOT defined as for PENTrack. The x direction is upward, the the coil planes are perpendicular to z.
    The z axies is the symetry one for tSPECT cylinder. x and z coordinate are swaped below, in exportMesh function.
    '''
    # C1
    phi_C1 = 120
    phi_off_C1 = 0
    length_C1 = 60
    x_off_C1 = 975
    r_C1 = 31.5
    cur_C1 = 22
    k_C1 = 100
    f_C1 = 38314864


    # C2
    phi_C2 = 120
    phi_off_C2 = 90
    length_C2 = 60
    x_off_C2 = 985
    r_C2 = 31.5
    cur_C2 = 22
    k_C2 = 100
    f_C2 = 38314864
    
    spinflipper_C1 = SF_field(phi_C1, phi_off_C1, length_C1, x_off_C1, r_C1, cur_C1*np.sin(f_C1*t), k_C1).sf_field()
    spinflipper_C2 = SF_field(phi_C2, phi_off_C2, length_C2, x_off_C2, r_C2, cur_C2*np.cos(f_C2*t), k_C2).sf_field()

    spinflipper = [spinflipper_C1, spinflipper_C2]
    spinflipper = mag.Collection(spinflipper, style_label="sf" + "str(SF)")
    

    return spinflipper
