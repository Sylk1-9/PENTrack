# import magpylib as magpy
# import numpy as np

# print("Script visited : PythonField in /in")
# take mm, output mT
def BField(source, xyz):
    # print(xyz)
    return source.getB(xyz)
