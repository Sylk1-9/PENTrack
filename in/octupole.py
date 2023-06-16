import magpylib as magpy
import numpy as np
from scipy.spatial.transform import Rotation as R

print("Script visited : octupole.py in in/")

# def BuildMagnet():
#     z = 2*np.sin(np.deg2rad(11.25)/2)*54
#     y = 84.1 - 54
#     x = 57.5
#     gap = 0 # in [deg]
#     a = 11.25 - gap 
#     N = 32 # number of Segments
#     ang=2*np.pi/N
#     k = 3
#     Br = 1100 # magnetic remanence 
#     #mu_0 = 1.25663706212e-06 # magnatic constant
#     M = Br #/mu_0
#     octupole = magpy.Collection()
#     for i in np.arange(0,24,1):
#         for deg in np.arange(0,360, 11.25 ):
#             magnetization = (M*np.cos((k+1)*np.deg2rad(deg)), M*np.sin((k+1)*np.deg2rad(deg)), 0)
#             dimension = (54, 84.1, x, -a/2, a/2)
#             position = (x/2+i*x + 570, 0, 0)
#             orientation = R.from_rotvec((0,90, 0), degrees=True)
#             color = 'rgb(%i, %i, %i)' % tuple(np.random.randint(0, 255, 3))
#             segment = magpy.magnet.CylinderSegment(magnetization, dimension, position, orientation, style={'color':color})
#             octupole.add(segment.rotate_from_angax(deg, axis='x',anchor=0))
#             # print(i, deg)

#     return octupole

def buildTrapez(magnetization=(1, 0, 0), dimension=(1, 3, 3, 4), position=(0, 0, 0), color="blue"):

    (a, b, h, l) = dimension
    v = []
    for z in [-l/2, l/2]:
        for x in [0, h]:
            ys = [-a/2, a/2] if (x==0) else [-b/2, b/2]
            for y in ys:
                v.append([x, y, z])

    v = np.array(v)

    vs = list(np.zeros((12, 3, 3)))
    vs[0]  = (v[3], v[0], v[1])
    vs[1]  = (v[2], v[0], v[3])
    vs[2]  = (v[1], v[0], v[4])
    vs[3]  = (v[4], v[0], v[2])
    vs[4]  = (v[2], v[6], v[4])
    vs[5]  = (v[3], v[6], v[2])
    vs[6]  = (v[7], v[6], v[3])
    vs[7]  = (v[4], v[6], v[7])
    vs[8]  = (v[7], v[5], v[4])
    vs[9]  = (v[4], v[5], v[1])
    vs[10] = (v[3], v[5], v[7])
    vs[11] = (v[1], v[5], v[3])

    trapez = magpy.Collection(style_label="Trapezoid")

    for ind, vert in enumerate(vs):
        trapez.add(
            magpy.misc.Triangle(
                magnetization=magnetization,
                vertices=vert,
                style_label=f"Triangle_{ind+1:02d}",
                style={'color':color} # 'rgb(%i, %i, %i)' % tuple(np.random.randint(0, 255, 3))}
            )
        )

    return trapez.move(position)



def buildSource(t=0, nseg=24, stype="cub", delta=0, nring=24):
    # nseg = number of segement along x.
    # nring = number of ring along x for the octupole (24 by default)
    # nseg=1 + nring=1  will display only 1 ring.
    # nseg=1 + nring=24 will display one long ring, the length of the octupole (equivalent to 24 join segments)
    z = 57.5*(nring/nseg)  # heigh of segment
    M = 1100 # Magnetization
    N = 32 # number of segment in circle, around z
    phi = 360/N # angle between segments
    r1 = 54
    r2 = 84.1
    area = np.pi*(r2*r2 - r1*r1)/N
    
    octupole = magpy.Collection()
    for i in np.arange(0, nseg):
        for deg in np.arange(0, 360, phi):
            rad = np.deg2rad(deg)
            magnetization = (M*np.cos(4*rad), M*np.sin(4*rad), 0)
            orientation = None # R.from_rotvec((0, 90, 0), degrees=True)
            color = 'rgb(%i, %i, %i)' % tuple(np.random.randint(0, 255, 3))

            if(stype == "cub"):
                y = 2*54*np.sin(np.deg2rad(phi)/2)
                x = area/y #84.1 - 54
                position = (54 + x/2, 0, 570 + z/2 + i*z)
                dimension = (x, y, z)
                # print("cube area ", x*y)
                segment = magpy.magnet.Cuboid(magnetization, dimension, position, orientation, {'color':color})
                segment.rotate_from_angax(deg, axis='z', anchor=0)
            elif(stype == "cyl"):
                position = (0, 0, 570 + z/2 + i*z)
                # dimension = (r1, r2, z, deg-phi/2, deg+phi/2)
                dimension = (r1, r2, z, -phi/2, phi/2)
                # print("cylinder area", np.pi*(r2*r2 - r1*r1)/N)
                segment = magpy.magnet.CylinderSegment(magnetization, dimension, position, orientation,{'color':color})
                segment.rotate_from_angax(deg, axis='z', anchor=0)
                
            elif(stype == "tra"):
                position = (54, 0, 570 + z/2 + i*z)
                a = 2*r1*np.sin(np.deg2rad(phi)/2) # base down
                aa = -1
                bb = -2*r1
                cc = area/a*r1*2
                k = np.sin(np.deg2rad(phi)/2) # cst
                h = (-bb - (bb*bb - 4*aa*cc)**.5)/(2*aa) # height
                b = 2*(r1+h)*np.sin(np.deg2rad(phi)/2) # base up
                l = z
                # print("trapez area", (a+b)/2 * h)
                dimension = (a, b, h, l)
                segment = buildTrapez(magnetization, dimension, position, color=color)
                segment.rotate_from_angax(deg, axis='z', anchor=0)
                # posxy = np.array(segment.position[:2])
                # print(np.sum(posxy**2))
            else:
                return

            segment.rotate_from_angax(np.random.rand()*delta, np.random.rand(3)*2-1, anchor=None)
            octupole.add(segment)
    return octupole

