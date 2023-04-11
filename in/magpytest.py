import magpylib as magpy

print("Script visited : magpytest.py in in/")
def BuildMagnet():
    print("Function visited : BuildMagnet in magpytest.py")
    cube = magpy.magnet.Cuboid(magnetization=(1000,1000,1000), dimension=(10,10,10))
    # print("in create magnet")
    return cube
