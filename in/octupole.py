import magpylib as magpy
import numpy as np

print("Script visited : octupole.py in in/")
def BuildMagnet():

    x = 2*np.sin(np.deg2rad(11.25)/2)*54 #339.29200658769764/32 #8.283127022 - 5.318535781
    y = 84.1 - 54
    z = 57.5

    # N = 32 # number of Segments
    # ang=2*np.pi/N

    k = 4

    Br = 1100 # magnetic remanence 
    #mu_0 = 1.25663706212e-06 # magnatic constant
    M = Br #/mu_0
    #m = np.sqrt(M**2/3)

    # dimension = (x, y, z)
    dimension = (x, y, z*24)

    octupole = magpy.Collection()
    # for i in np.arange(0, 24, 1):
        # for deg in np.arange(0, 360, 11.25 ):
        #     magnetization = (M*np.cos((k+1)*np.deg2rad(deg)), M*np.sin((k+1)*np.deg2rad(deg)),0)
        #     position = (0, 54 + y/2, z/2+i*57.5 + 570)
        #     cube = magpy.magnet.Cuboid(magnetization=magnetization, dimension=dimension, position=position)
        #     octupole.add(cube.rotate_from_angax(deg, axis='z', anchor=0))

    for deg in np.arange(0, 360, 11.25 ):
        rad = np.deg2rad(deg)
        magnetization = (M*np.cos((k+1)*rad), M*np.sin((k+1)*rad), 0)
        position = (0, 54 + y/2, z*24/2 + 570)
        cube = magpy.magnet.Cuboid(magnetization=magnetization, dimension=dimension, position=position)
        octupole.add(cube.rotate_from_angax(deg, axis='z', anchor=0))

            
    return octupole
