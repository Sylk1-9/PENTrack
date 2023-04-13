# import magpylib as magpy
# import numpy as np

print("Script visited : -- -MagpyField in /in")
# take mm, output mT
def BField(magnet, xyz):
    # print(xyz)
    return magnet.getB(xyz)
