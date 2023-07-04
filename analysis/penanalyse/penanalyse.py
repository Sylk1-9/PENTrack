import os as os
import glob as glob
import scipy as scipy
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
import time

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
class SetVisibilityCallback:
    """Helper callback to keep a reference to the actor being modified."""

    def __init__(self, actors):
        self.actors = actors

    def __call__(self, state):
        if isinstance(self.actors, (list, tuple, set, dict)):
            for actor in self.actors:
                actor.SetVisibility(state)
        else:
            self.actors.SetVisibility(state)

            
class data:    
    def __init__(self, dfile="000000000002", dfilext=".h5", keys=('nt'), loadstl=True, plotter=None, ptdir="/home/sly/Work/Physics/Neutrons/tSPECT/PENtrack/pentrack/"):
        self.ptdir = ptdir
        self.keys = keys
        self.dfile = dfile
        self.dfilext = dfilext
        self.df = {}
        self.stlfiles = []
        self.stlactors = []
        self.logsactor = {}
        self.loadata()
        if loadstl:
            self.loadstl()
        # self.axes = None
        self.plotter = plotter
            
    def loadata(self):
        if (("h5" in self.dfilext) or ("hdf5" in self.dfilext))  :
            fh5 = h5py.File(self.ptdir + "out/" + self.dfile + ".h5",'r')
            for key in fh5.keys():
                if (key not in 'config') and (key in keydic):
                    # print("loading " + key + " in dataframe df")
                    
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

            # # Todo :  iterator for large file?
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
        plotter = pvqt.BackgroundPlotter()
        plotter.show_bounds(color="black")
        return plotter

    
    def opacity_stl(self, opacity):
        for actor in self.stlactors:
            actor.GetProperty().SetOpacity(opacity)

            
    def plotstl(self, plotter=None, opacity=1, stlfiles=None):

        self.plotter = self.initplotter() if self.plotter is None else self.plotter
        self.stlactors = []
        stlfiles = self.stlfiles if stlfiles is None else stlfiles

        size = 30
        pos = 2*size
        colors = get_distinct_colors(len(stlfiles))

        
        for stlfile, color in zip(stlfiles, colors):
            mesh = pv.read(stlfile).scale([1000.0, 1000.0, 1000.0], inplace=False)
            # mesh = mesh.rotate_y(-90, point=axes.origin, inplace=False)
            actor = self.plotter.add_mesh(mesh, opacity=opacity, color=color)
            
            self.stlactors.append(actor)
            
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
        
        self.plotter.add_slider_widget(self.opacity_stl, (0.0, 1.0), value=opacity, pointa=(0.04, 0.6), pointb=(0.04, 0.9),  slider_width=0.02, tube_width=0.005, color="black")

        return self.plotter, self.stlactors


    def opacity_mag(self, opacity):
        for actor in self.magactors:
            actor.GetProperty().SetOpacity(opacity)

            
    def plotmag(self, opacity=1):
        
        self.plotter = self.initplotter() if self.plotter is None else self.plotter
        self.stlactors = []
      
        co = magt.coil(wx=20, wr=50).source.rotate_from_angax(90, axis='y', anchor=0)
        octu = magt.octupole(nzseg=1, nring=24, stype="cyl").source.rotate_from_angax(90, axis='y', anchor=0)
        comp = magt.octupole(nzseg=1, nring=5, stype="cyl", comp=True).source.rotate_from_angax(90, axis='y', anchor=0)
        sf1 = magt.spinflipper(isf=1, kseg=20).source#.rotate_from_angax(-90, axis='y', anchor=0)
        sf2 = magt.spinflipper(isf=2, kseg=20).source#.rotate_from_angax(-90, axis='y', anchor=0)

        size = 30
        pos = 2*size
        sources = [co, octu, comp, sf1, sf2]
        prevactors = self.plotter.actors.copy()
        colors = get_distinct_colors(len(sources))
        for source, color in zip(sources, colors):
            source.set_children_styles(color=color)
            source.set_children_styles(opacity=opacity)
             
            oldvactors = self.plotter.actors.copy()
            self.plotter = source.show(backend='pyvista', canvas=self.plotter, return_fig=True, style_magnetization_show=False)
            actor = [v for k,v in self.plotter.actors.items() if k not in oldvactors]

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
        self.plotter.add_slider_widget(self.opacity_mag, (0, 1), value=opacity, pointa=(0.6, 0.03), pointb=(0.9, 0.03), slider_width=0.02, tube_width=0.005, color="black")

        return self.plotter, self.magactors


    def plotlogs(self, state="start", ptype="neutron", color=None):
        '''
        state = "start", "end", or "hit"
        '''
        st = '' if "hit" in state else state
        ltype = 'h' if "hit" in state else 'e'
        if (ptype == 'n' or ptype == "neutron"):
            key = 'n'
        elif (ptype == 'p' or ptype == "proton"):
            key = 'p'
        elif (ptype == 'e' or ptype == "electron"):
            key = 'e'

        if color is None:
            color = "black" if state == "start" else "red"
        
        self.logsactor[key+ltype+state] = self.plotter.add_mesh(1000 * da.df[key+ltype][['x'+st, 'y'+st, 'z'+st]].values, color=color, render_points_as_spheres=True, point_size=6)

        return self.plotter
    

    def animate(self, dur=1, dt=0.005, framerate=50, savegif=False):
        if savegif:
            self.plotter.open_gif("/home/sly/Work/Physics/Neutrons/tSPECT/PENtrack/plots/animations/penanalyse/animation-%s-2.gif" % dfile)

        df_nts = self.df['nt'].sort_values(['particle', 't'])

        points = pv.PolyData(1000 * self.df['ne'].sort_values('particle')[['xstart', 'ystart', 'zstart']].values)
        points['particle'] = self.df['ne'].sort_values('particle')['particle']
        points['t'] = self.df['ne']['tstart']

        cloud = pv.PolyData(points)
        self.plotter.add_points(cloud, color='blue', render_points_as_spheres=True, point_size=6)
        # particles = df_nts['particle'].unique()

        interpolation_functions = {}
        for particle, data in df_nts.groupby('particle'):
            positions = data[['x', 'y', 'z']].values
            times = data['t'].values
            interpolation_functions[particle] = scipy.interpolate.interp1d(times, positions, axis=0, bounds_error=False, fill_value=(data[['x', 'y', 'z']].iloc[0], data[['x', 'y', 'z']].iloc[-1]))


        self.time_stamp = 0
        def update_animation():
            self.time_stamp += dt
            interpolated_positions = []
            for particle, interpolation_func in interpolation_functions.items():
                interpolated_positions.append(interpolation_func(self.time_stamp))

            cloud.points = 1000 * np.array(interpolated_positions)
            self.plotter.update()
            if savegif:
                self.plotter.write_frame()

        
        self.plotter.add_callback(update_animation, int(1000/framerate), int(dur*framerate/dt))
        # self.plotter.add_callback(update_animation, int(0.001/framerate), int(200/(framerate*dt)))
        return self.plotter



dfile = "000000000006"
da = data(dfile)
pl, magactors = da.plotmag(opacity=0.01)
pl, stlactors = da.plotstl(opacity=0.01)
# pl = da.plotends("end"")
pl = da.plotlogs(ptype="n", state="start", color="black")
pl = da.plotlogs(ptype="n", state="start", color="pink")
pl = da.plotlogs(ptype="n", state="hit", color="red")
pl = da.animate(dur=10, framerate=30, dt=0.005, savegif=True)

# gif = False
# if gif:
#     pl.open_gif("/home/sly/Work/Physics/Neutrons/tSPECT/PENtrack/plots/animations/penanalyse/animation-%s-2.gif" % dfile)

# df_nts = da.df['nt'].sort_values(['particle', 't'])

# points = pv.PolyData(1000 * da.df['ne'].sort_values('particle')[['xstart', 'ystart', 'zstart']].values)
# points['particle'] = da.df['ne'].sort_values('particle')['particle']
# points['t'] = da.df['ne']['tstart']

# cloud = pv.PolyData(points)
# pl.add_points(cloud, color='blue', render_points_as_spheres=True, point_size=6)

# particles = df_nts['particle'].unique()


# interpolation_functions = {}
# for particle, data in df_nts.groupby('particle'):
#     positions = data[['x', 'y', 'z']].values
#     times = data['t'].values
    
#     interpolation_functions[particle] = scipy.interpolate.interp1d(times, positions, axis=0, bounds_error=False, fill_value=(data[['x', 'y', 'z']].iloc[0], data[['x', 'y', 'z']].iloc[-1]))


   
# time_stamp = 0.0
# dt = 0.005
# framerate = 5 # s
# def update_animation():
#     global time_stamp
#     time_stamp += dt
#     interpolated_positions = []
#     for particle, interpolation_func in interpolation_functions.items():
#         interpolated_positions.append(interpolation_func(time_stamp))

#     cloud.points = 1000 * np.array(interpolated_positions)
#     pl.update()
#     if gif:
#         pl.write_frame()


# pl.add_callback(update_animation, int(0.001/framerate), int(200/(framerate*dt)))

# # pl.close()

