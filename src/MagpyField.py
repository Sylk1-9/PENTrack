import magpylib as magpy

print("Script visited : MagpyField in /src")
def getB(magnet, x, y, z, t=0):
    sensor = magpy.Sensor(position=(x,y,z), style_size=1.8)
    B = sensor.getB(magnet)
    # print("magnet = ", magnet)
    # print("result should be = ", B)
    # return [x, y, z]
    return list(B)
