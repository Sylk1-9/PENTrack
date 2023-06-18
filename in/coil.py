import magpylib as magpy
import numpy as np
from scipy.spatial.transform import Rotation as R

print("Script visited : coil.py in in/")
# def BuildMagnet():
#     '''
#     the axis is defined as for PENTrack. The z direction is upward, the the coil planes are perpendicular to x
#     x is the axe of symetry of the tSPECT cylinder.
#     '''
#     I_main = 330000 # 700000 
#     I_c3 = 450000 #350000
#     I_c5 = -200000#150000
#     I_AHC = 255000 #255000
    
#     # coils = ['coil1', 'coil2a', 'coil2b', 'coil2c', 'coil4a', 'coil4b', 'coil4c', 'coil6', 'coil7', 'coil8', 'coil9', 'coil10', 'coil11', 'coil3', 'coil5', 'coil12', 'coil13']
    
#     J = [256.527680856013*I_main, 146.853146853147*I_main, 189.810189810189*I_main, 256.463597424696*I_main, 146.901709401710*I_main, 
#          189.255189255189*I_main ,256.410256410256*I_main, 256.502028155571*I_main, 256.410256410257*I_main
#          , 256.410256410256*I_main, 256.391089849006*I_main, 256.410256410256*I_main, 256.324900133156*I_main
#          , 256.410256410256*I_c3, 256.410256410256*I_c5, 3.04761904761905*I_AHC, -3.04761904761905*I_AHC]
    
#     width = [4.8, 16, 7, 26.33, 17.6, 12.6, 34.8, 4.8, 7.8, 9, 1.8, 19.8, 27, 4.2, 4.2, 21, 21]
#     height = [286.6, 71.5, 71.5, 71.5, 46.8, 46.8, 46.8, 139.7, 178.1, 315.9, 668.9, 302.9, 150.2, 4.55, 4.55, 50, 50]
#     r = [129, 254.8, 270.8, 277.8, 242.8, 260.4, 273, 140.8, 150, 150, 192.8, 120, 120, 131.1, 131.3, 374, 374]
#     x = [-429, -137, -137, -137, 67.3, 67.3, 67.3, 139.6, 295.4, 474.1, 990.6, 1900, 2202.9, -117.8, 113.2, 1020, 1565]
#     s = 100
    
#     coil = magpy.Collection()
#     for i in np.arange(len(J)):
#         for j in np.linspace(x[i], x[i] + height[i], s):
#             for k in np.arange(r[i], r[i] + width[i], s):
#                 current = J[i]*width[i]*height[i]/(1000**2*s)
#                 orientation = R.from_rotvec((0,90,0), degrees=True)
#                 coil.add(magpy.current.Loop(current=current, diameter=k*2, position=(j,0,0), orientation=orientation))

#     return coil


def buildSource(t=0, winding=100):
    '''
    The axis is NOT defined as for PENTrack. The x direction is upward, the the coil planes are perpendicular to z.
    The z axies is the symetry one for tSPECT cylinder. x and z coordinate are swaped below, in exportMesh function.
    '''
    I_main = 330000 # 700000 
    I_c3 = 450000 #350000
    I_c5 = -200000#150000
    I_AHC = 255000 #255000
    
    J = [256.527680856013*I_main, 146.853146853147*I_main, 189.810189810189*I_main, 256.463597424696*I_main, 146.901709401710*I_main, 
         189.255189255189*I_main ,256.410256410256*I_main, 256.502028155571*I_main, 256.410256410257*I_main
         , 256.410256410256*I_main, 256.391089849006*I_main, 256.410256410256*I_main, 256.324900133156*I_main
         , 256.410256410256*I_c3, 256.410256410256*I_c5, 3.04761904761905*I_AHC, -3.04761904761905*I_AHC]
    
    width = [4.8, 16, 7, 26.33, 17.6, 12.6, 34.8, 4.8, 7.8, 9, 1.8, 19.8, 27, 4.2, 4.2, 21, 21]
    height = [286.6, 71.5, 71.5, 71.5, 46.8, 46.8, 46.8, 139.7, 178.1, 315.9, 668.9, 302.9, 150.2, 4.55, 4.55, 50, 50]
    r = [129, 254.8, 270.8, 277.8, 242.8, 260.4, 273, 140.8, 150, 150, 192.8, 120, 120, 131.1, 131.3, 374, 374]
    x = [-429, -137, -137, -137, 67.3, 67.3, 67.3, 139.6, 295.4, 474.1, 990.6, 1900, 2202.9, -117.8, 113.2, 1020, 1565]
    # s = 100
    
    coil = magpy.Collection()
    for i in np.arange(len(J)):
        color = 'rgb(%i, %i, %i)' % tuple(np.random.randint(0, 255, 3))
        for j in np.linspace(x[i], x[i] + height[i], winding):
            for k in np.arange(r[i], r[i] + width[i], winding):
                current = J[i]*width[i]*height[i]/(1000**2*winding)
                orientation = None# R.from_rotvec((0,90,0), degrees=True)
                # orientation = R.from_rotvec((0,90,0), degrees=True)
                # coil.add(magpy.current.Loop(current=current, diameter=k*2, position=(j,0,0), orientation=orientation))
                coil.add(magpy.current.Loop(current=current, diameter=k*2, position=(0,0,j), orientation=orientation,style={'color':color}))
                # rotate to get octupole along x axis.
                coil.rotate_from_angax(angle=90, anchor=(0,0,0), axis='y')
    

    return coil
