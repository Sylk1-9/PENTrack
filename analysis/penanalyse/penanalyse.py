import os as os
import glob as glob
import h5py as h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv
# from pyvistaqt import BackgroundPlotter
import pyvistaqt as pvqt
import vtkplotlib as vpl
import stl as stl
from stl.mesh import Mesh
import magpylib as magpy
import magtspect as magt
from colorsys import hls_to_rgb
import numpy as np 


def get_distinct_colors(n):
    colors = []
    for i in np.arange(0., 360., 360. / n):
        h = i / 360.
        l = (50 + np.random.rand() * 10) / 100.
        s = (90 + np.random.rand() * 10) / 100.
        colors.append(hls_to_rgb(h, l, s))
    return colors


plt.ion()
plt.show()

keydic = {
        'neutrontrack':'nt', 'neutronhit':'nh', 'neutronend':'ne', 'neutronsnapshot':'ns', 'neutronspin':'nS',
        'protontrack':'pt', 'protonhit':'ph', 'protonend':'pe', 'protonsnapshot':'ps', 'protonspin':'pS',
        'electrontrack':'et', 'electronhit':'eh', 'electronend':'ee', 'electronsnapshot':'es', 'electronspin':'eS',
}

magpy.defaults.display.style.current.arrow.width = 4


# stl_file_paths = glob.glob(workingdir + "in/CAD/*.stl") + glob.glob(workingdir + "in/CAD/store/*.stl")
# class SetVisibilityCallback:
#     """Helper callback to keep a reference to the actor being modified."""

#     def __init__(self, actors):
#         self.actors = actors

#     def __call__(self, state):
#         if isinstance(self.actors, (list, tuple, set, dict)):
#             for actor in self.actors:
#                 actor.SetVisibility(state)
#         else:
#             self.actors.SetVisibility(state)



class MyCustomRoutine:
    def __init__(self, mesh):
        self.output = mesh  # Expected PyVista mesh type
        # default parameters
        self.kwargs = {
            'radius': 0.5,
            'theta_resolution': 30,
            'phi_resolution': 30,
        }

    def __call__(self, param, value):
        self.kwargs[param] = value
        self.update()

    def update(self):
        # This is where you call your simulation
        result = pv.Sphere(**self.kwargs)
        self.output.copy_from(result)
        return

            
class data:    

    def __init__(self, dfile="000000000001", dfilext=".h5", keys=('nt'), loadstl=True, ptdir="/home/sly/Work/Physics/Neutrons/tSPECT/PENtrack/pentrack/"):
        self.ptdir = ptdir
        self.keys = keys
        self.dfile = dfile
        self.dfilext = dfilext
        self.df = {}
        self.stlfiles = []
        self.stlactors = []
        # self.fig = fig
        # self.plotter = plotter
        self.loadata()
        if loadstl:
            self.loadstl()
        self.axes = None
        self.plotter = None
            
    def loadata(self):
        if (("h5" in self.dfilext) or ("hdf5" in self.dfilext))  :
            fh5 = h5py.File(self.ptdir + "out/" + self.dfile + ".h5",'r')
            for key in fh5.keys():
                if (key not in 'config') and (key in keydic):
                    print("loading " + key + " in dataframe df")
                    
                    self.df[keydic[key]] = pd.DataFrame(fh5[key][:])
                else:
                    self.df["config"] = fh5["config"]

        elif (("out" in self.dfilext) or ("txt" in self.dfilext)):
            for key in keydic:
                outpath = self.ptdir + "out/" + self.dfile + key + ".out"
                if os.path.exists(outpath):
                    outpathc = outpath
                    print("loading " + key + " in dataframe df")
                    self.df[keydic[key]] = pd.read_csv(outpathc, sep=" ")
        
            # spinlist = []
            # data_iter = pd.read_hdf(workingdir + "out/saved/" + config_name + ".h5", key='neutronspin', columns=['x', 'y', 'z', 'Sx', 'Sy', 'Sz'], iterator=True, chunksize=100, where=None, start=50000000)
            # it = 0
            # for chunk in data_iter:
            #    if (it > 10):
            #       break
    
            #    it = it+1
            #    #train cnn on chunk here
            #    spinlist.append(chunk[['x', 'y', 'z', 'Sx', 'Sy', 'Sz']])
            #    # print(chunk[['x', 'y', 'z', 'Sx', 'Sy', 'Sz']])
            #    # print(chunk)
            #    print(spinlist)

        elif ".root" in self.dfile:
            print("not implemented")
            return

    def loadstl(self, config=None):
        if config is None:
            config = self.df["config"]
            
        geolist = config["GEOMETRY"]
        for geo in geolist:
            strgeo = str(geo[1]).split("\\t")
            if "ignored" not in strgeo[1]:
                self.stlfiles.append(self.ptdir + "in/" + strgeo[1])

    def loadallstl(self):
        stl_file_paths = glob.glob(self.ptdir + "in/STLs/*.stl") + glob.glob(self.ptdir + "in/STLs/store/*.stl")
        for geo in geolist:
            strgeo = str(geo[1]).split("\\t")
            if "ignored" not in strgeo[1]:
                self.allstlfiles.append(self.ptdir + "in/" + strgeo[1])


    def initplotter(self):
        # axes = pv.Axes(show_actor=True, actor_scale=2.0, line_width=5)
        plotter = pvqt.BackgroundPlotter()
        plotter.show_bounds(color="black")
        return plotter
        
    # def toggle_stl(flag):
    #     for actor in self.stlactors:
    #         actor.SetVisibility(flag)
        
    # def toggle_actor(self, flag):
    #     for actor in self.stlactors:
    #         actor.SetVisibility(flag)
    def opacity_stl(self, opacity):
        for actor in self.stlactors:
            actor.GetProperty().SetOpacity(opacity)
                
    def plotstl(self, plotter=None, opacity=1, stlfiles=None):

        # self.axes = pv.Axes(show_actor=True, actor_scale=2.0, line_width=5) if self.axes is None else self.axes
        self.plotter = self.initplotter() if self.plotter is None else self.plotter
        self.stlactors = []
        stlfiles = self.stlfiles if stlfiles is None else stlfiles

        size = 30
        pos = 2*size
        colors = get_distinct_colors(len(stlfiles))

        
        for stlfile, color in zip(stlfiles, colors):
            print(stlfile)
            # color= tuple(np.random.randint(0, 255, 3)) + (1,)
            mesh = pv.read(stlfile).scale([1000.0, 1000.0, 1000.0], inplace=False)
            # engine = MyCustomRoutine(mesh)
            # mesh = mesh.rotate_y(-90, point=axes.origin, inplace=False)
            actor = self.plotter.add_mesh(mesh, opacity=opacity, color=color)
            
            self.stlactors.append(actor)
            
            # Make a separate callback for each widget
            callback = SetVisibilityCallback(actor)
            self.plotter.add_checkbox_button_widget(
                callback,
                value=True,
                position=(5.0, pos),
                size=size,
                border_size=1,
                color_on=color,
                color_off='grey',
                background_color='grey',
            )
            pos = pos + size + (size // 10)
        
        # self.plotter.add_checkbox_button_widget(self.toggle_actor, size=30, position=(10, 10), color_on="blue", value=True)
        self.plotter.add_slider_widget(self.opacity_stl, (0.0, 1.0), pointa=(0.04, 0.6), pointb=(0.04, 0.9),  slider_width=0.02, tube_width=0.005, color="black")

                
        return self.plotter, self.stlactors

    # def toggle_mag(flag):
    #     for actor in self.magactors:
    #         self.magactors[actor].SetVisibility(flag)

    def opacity_mag(self, opacity):
        for actor in self.magactors:
            actor.GetProperty().SetOpacity(opacity)

            
    def plotmag(self, opacity=1):

        # magpy.defaults.display.style.base.opacity = opacity
        
        # self.axes = if self.axes is None else self.axes
        self.plotter = self.initplotter() if self.plotter is None else self.plotter
        self.stlactors = []
      
        co = magt.coil(wx=20, wr=50).source.rotate_from_angax(90, axis='y', anchor=0)
        octu = magt.octupole(nzseg=1, nring=24, stype="cyl").source.rotate_from_angax(90, axis='y', anchor=0)
        comp = magt.octupole(nzseg=1, nring=5, stype="cyl", comp=True).source.rotate_from_angax(90, axis='y', anchor=0)
        sf1 = magt.spinflipper(isf=1, kseg=20).source#.rotate_from_angax(-90, axis='y', anchor=0)
        sf2 = magt.spinflipper(isf=2, kseg=20).source#.rotate_from_angax(-90, axis='y', anchor=0)

        # coll.set_children_styles(magnetization_color_south="blue")
        
        size = 30
        pos = 2*size
        sources = [co, octu, comp, sf1, sf2]
        prevactors = self.plotter.actors.copy()
        colors = get_distinct_colors(len(sources))
        for source, color in zip(sources, colors):
            # color = 'rgb(%i, %i, %i)' % tuple(np.random.randint(0, 255, 3))
            # color= tuple(np.random.randint(0, 255, 3),) # + (1,)
            # source.set_children_styles(color='rgb(%i, %i, %i)' % color)
            source.set_children_styles(color=color)
            source.set_children_styles(opacity=opacity)
             
            oldvactors = self.plotter.actors.copy()
            self.plotter = source.show(backend='pyvista', canvas=self.plotter, return_fig=True, style_magnetization_show=False)
            actor = [v for k,v in self.plotter.actors.items() if k not in oldvactors]

            # Make a separate callback for each widget
            callback = SetVisibilityCallback(actor)
            self.plotter.add_checkbox_button_widget(
                callback,
                value=True,
                position=(pos, 5.0),
                size=size,
                border_size=1,
                color_on=color,
                color_off='grey',
                background_color='grey',
            )
            pos = pos + size + (size // 10)

        self.magactors = [v for k,v in self.plotter.actors.items() if k not in prevactors]
        # self.plotter.add_checkbox_button_widget(self.toggle_actor, size=30, position=(40, 10), color_on="green", value=True)

        self.plotter.add_slider_widget(self.opacity_mag, (0, 1), pointa=(0.6, 0.03), pointb=(0.9, 0.03), slider_width=0.02, tube_width=0.005, color="black")

        return self.plotter, self.magactors



da = data()
pl, magactors = da.plotmag(opacity=0.9)
pl, stlactors = da.plotstl(opacity=0.3)

# for actor in stlactors:
#     actor.GetProperty().SetOpacity(0.1)

####################

# magpy.defaults.display.style.magnet.magnetization.color.mode = 'bicolor'
# magpy.defaults.display.style.base.path.line.width = 10
# magpy.defaults.display.style.markers.path.line.width

# axes = pv.Axes(show_actor=True, actor_scale=1.0, line_width=2)

#######################

# pl = pv.Plotter()

# co = magt.coil(wx=200, wr=50).source
# octu = magt.octupole(nzseg=1, stype="cyl").source
# comp = magt.octupole(nzseg=1, stype="cyl", comp=True).source
# sf1 = magt.spinflipper(isf=1, kseg=20).source.rotate_from_angax(-90, axis='y', anchor=0)
# sf2 = magt.spinflipper(isf=2, kseg=20).source.rotate_from_angax(-90, axis='y', anchor=0)
# coll = magpy.Collection([co, octu, comp, sf1, sf2])
# # coll = magpy.Collection([sf1, sf2])

# pl = coll.show(backend='pyvista', plotter=pl, return_fig=True, style_magnetization_show=False)
# da = data()
# da.plotstl(plotter=pl, axes=axes, opacity=0.2)
# pl.show_bounds(color="black")

# # pl.show(interactive_update=True, interactive=True, auto_close=False)
# pl.show(interactive_update=False, interactive=True, auto_close=False)
# # pl.show(interactive_update=True)



#######################

# pl = BackgroundPlotter()

# pl = pv.Plotter()
# pl = BackgroundPlotter()

# co = magt.coil(wx=20, wr=50).source.rotate_from_angax(90, axis='y', anchor=0)
# octu = magt.octupole(nzseg=1, stype="cyl").source.rotate_from_angax(90, axis='y', anchor=0)
# comp = magt.octupole(nzseg=1, stype="cyl", comp=True).source.rotate_from_angax(90, axis='y', anchor=0)
# sf1 = magt.spinflipper(isf=1, kseg=20).source#.rotate_from_angax(-90, axis='y', anchor=0)
# sf2 = magt.spinflipper(isf=2, kseg=20).source#.rotate_from_angax(-90, axis='y', anchor=0)

# # coll = magpy.Collection([co, octu, comp, sf1, sf2])
# # coll = magpy.Collection([sf1, sf2])

# pl = co.show(backend='pyvista', canvas=pl, return_fig=True, style_magnetization_show=False)
# pl = octu.show(backend='pyvista', canvas=pl, return_fig=True, style_magnetization_show=False)
# pl = comp.show(backend='pyvista', canvas=pl, return_fig=True, style_magnetization_show=False)
# pl = sf1.show(backend='pyvista', canvas=pl, return_fig=True, style_magnetization_show=False)
# pl = sf2.show(backend='pyvista', canvas=pl, return_fig=True, style_magnetization_show=False)
# # pl = coll.show(backend='pyvista', canvas=pl, return_fig=True, style_magnetization_show=False)
# magactors = pl.actors.copy()

# da = data()
# pl, stlactors = da.plotstl(plotter=pl, axes=axes, opacity=0.5)
    
# def toggle_stl(flag):
#     for actor in stlactors:
#         actor.SetVisibility(flag)

# def toggle_mag(flag):
#     for actor in magactors:
#         magactors[actor].SetVisibility(flag)


# pl.add_checkbox_button_widget(toggle_stl, position=(10, 10), value=True)
# pl.add_checkbox_button_widget(toggle_mag, position=(100, 10), value=True)
# pl.show_bounds(color="black")

# # pl.enable_terrain_style()
# pl.show()

# pl.update()
# # pl.show(interactive_update=True, interactive=True, auto_close=False)
# # pl.show(interactive_update=False, interactive=True, auto_close=False)
# # pl.show(interactive_update=True)


