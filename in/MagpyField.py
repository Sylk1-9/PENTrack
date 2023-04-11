import magpylib as magpy
import numpy as np

print("Script visited : MagpyField in /in")
def BField(magnet, x, y, z):
    # print("Function visited : BField in MagpyField.py")

    # B = 1.0/1000 * magpy.getB(magnet, position=1000*np.array([x, y, z]).T, sumup=True)
    B = magpy.getB(magnet, 1000*np.array([x, y, z]), sumup=True)
    # print("magnet = ", magnet)
    # print("result should be = ", B)
    # return [x, y, z]
    return list(B/1000.0)
