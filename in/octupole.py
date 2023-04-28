import magpylib as magpy
import numpy as np
from scipy.spatial.transform import Rotation as R

print("Script visited : octupole.py in in/")
def BuildMagnet():
    z = 2*np.sin(np.deg2rad(11.25)/2)*54
    y = 84.1 - 54
    x = 57.5
    gap = 0 # in [deg]
    a = 11.25 - gap 
    N = 32 # number of Segments
    ang=2*np.pi/N
    k = 3
    Br = 1100 # magnetic remanence 
    #mu_0 = 1.25663706212e-06 # magnatic constant
    M = Br #/mu_0
    octupole = magpy.Collection()
    for i in np.arange(0,24,1):
        for deg in np.arange(0,360, 11.25 ):
            magnetization = (M*np.cos((k+1)*np.deg2rad(deg)), M*np.sin((k+1)*np.deg2rad(deg)), 0)
            dimension = (54, 84.1, x, -a/2, a/2)
            position = (x/2+i*x + 570, 0, 0)
            orientation = R.from_rotvec((0,90, 0), degrees=True)
            color = 'rgb(%i, %i, %i)' % tuple(np.random.randint(0, 255, 3))
            segment = magpy.magnet.CylinderSegment(magnetization, dimension, position, orientation, style={'color':color})
            octupole.add(segment.rotate_from_angax(deg, axis='x',anchor=0))
            # print(i, deg)

    return octupole
