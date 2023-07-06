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
from datetime import datetime
# import imageio.v3 as iio
import pygifsicle as pygifsicle


def get_distinct_colors(n):
    colors = []
    for i in np.arange(360., 0., -360. / n):
        h = i / 360.
        l = (50 + np.random.rand() * 10) / 100.
        s = (90 + np.random.rand() * 10) / 100.
        colors.append(hls_to_rgb(h, l, s))
    return colors


plt.ion()
plt.show()

lighting = True

keydic = {
        'neutrontrack':'nt', 'neutronhit':'nh', 'neutronend':'ne', 'neutronsnapshot':'ns', 'neutronspin':'nS',
        'protontrack':'pt', 'protonhit':'ph', 'protonend':'pe', 'protonsnapshot':'ps', 'protonspin':'pS',
        'electrontrack':'et', 'electronhit':'eh', 'electronend':'ee', 'electronsnapshot':'es', 'electronspin':'eS',
}

magpy.defaults.display.style.current.arrow.width = 4

# stl_file_paths = glob.glob(workingdir + "in/CAD/*.stl") + glob.glob(workingdir + "in/CAD/store/*.stl")

            
class data:    
    def __init__(self, dfile="000000000002", dfilext=".h5", keys=('nt'), loadstl=True, plotter=None, ptdir="/home/sly/Work/Physics/Neutrons/tSPECT/PENtrack/pentrack/"):
        self.ptdir = ptdir
        self.keys = keys
        self.dfile = dfile
        self.dfilext = dfilext
        self.df = {}
        self.stlfiles = []
        # self.plotter.stlactors = []
        # self.plotter.logsactor = {}
        # self.pactor = None
        self.loadata()
        self.simcount = int(str(self.df['config']["GLOBAL"][3][1]).split(" ")[1].split("\\")[0])
        self.simtime = int(str(self.df['config']["GLOBAL"][4][1]).split(" ")[1].split("\\")[0])
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
        # stl_file_paths = glob.glob(self.ptdir + "in/STLs/*.stl") + glob.glob(self.ptdir + "in/STLs/store/*.stl")
        for geo in geolist:
            strgeo = str(geo[1]).split("\\t")
            if "ignored" not in strgeo[1]:
                self.allstlfiles.append(self.ptdir + "in/" + strgeo[1])

                
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

                
    class ScreenshotCallback:
        """Helper callback to keep a reference to the actor being modified."""

        def __init__(self, plotter, dfile, spath):
            self.plotter = plotter
            self.dfile = dfile
            self.spath = spath
        def __call__(self, state):
            if state:
                current_datetime = datetime.now()
                formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H-%M-%S")
                filename = self.spath + "/screenshot-%s-%s.png" % (self.dfile, formatted_datetime)
                self.plotter.screenshot(filename)
            
    class RecordGifCallback:
        """Helper callback to keep a reference to the actor being modified."""

        def __init__(self, plotter, dfile, spath, fps):
            self.plotter = plotter
            self.dfile = dfile
            self.spath = spath
            self.gifname = None
            self.fps = fps
        def __call__(self, state):
            self.plotter.savegif = state
            if state:
                print("recording gif")
                current_datetime = datetime.now()
                formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H-%M-%S")
                self.gifname = self.spath + "/animation-%s-%s.gif" % (self.dfile, formatted_datetime)
                self.plotter.open_gif(self.gifname, subrectangles=True, fps=self.fps)
            else:
                print("saving and converting gif (please wait)")
                self.plotter.open_gif(self.spath + "/dummy.gif")
                if self.gifname is not None:
                    pygifsicle.optimize(self.gifname) # For overwriting the original one

    
    def initplotter(self):
        plotter = pvqt.BackgroundPlotter(window_size=(2*1024, 2*768), multi_samples=2)
        plotter.show_bounds(color="black")
        return plotter

    
    def opacity_stl(self, opacity):
        for actor in self.plotter.stlactors:
            actor.GetProperty().SetOpacity(opacity)

            
    def plotstl(self, opacity=1, stlfiles=None):

        self.plotter = self.initplotter() if self.plotter is None else self.plotter
        self.plotter.stlactors = []
        stlfiles = self.stlfiles if stlfiles is None else stlfiles

        size = 30
        pos = 2*size
        colors = get_distinct_colors(len(stlfiles))

        for stlfile, color in zip(stlfiles, colors):
            mesh = pv.read(stlfile).scale([1000.0, 1000.0, 1000.0], inplace=False)
            # mesh = mesh.rotate_y(-90, point=axes.origin, inplace=False)
            actor = self.plotter.add_mesh(mesh, opacity=opacity, color=color, lighting=lighting)
            
            self.plotter.stlactors.append(actor)
            
            callback = self.SetVisibilityCallback(actor)
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
        
        self.plotter.add_slider_widget(self.opacity_stl, (0.0, 1.0), value=opacity, pointa=(0.04, 0.6), pointb=(0.04, 0.98),  slider_width=0.03, tube_width=0.006, color="black")

        return self.plotter


    def opacity_mag(self, opacity):
        for actor in self.plotter.magactors:
            actor.GetProperty().SetOpacity(opacity)

            
    def plotmag(self, opacity=1):
        
        self.plotter = self.initplotter() if self.plotter is None else self.plotter
        self.plotter.magactors = []
      
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

            callback = self.SetVisibilityCallback(actor)
            self.plotter.add_checkbox_button_widget(
                callback,
                value=True,
                position=(40, pos),
                size=size,
                border_size=1,
                color_on=color,
                color_off='grey',
                background_color='grey',
            )
            pos = pos + size + (size // 10)

        self.plotter.magactors = [v for k,v in self.plotter.actors.items() if k not in prevactors]
        # self.plotter.add_slider_widget(self.opacity_mag, (0, 1), value=opacity, pointa=(0.6, 0.03), pointb=(0.9, 0.03), slider_width=0.02, tube_width=0.005, color="black")
        self.plotter.add_slider_widget(self.opacity_mag, (0, 1), value=opacity, pointa=(0.1, 0.6), pointb=(0.1, 0.98),  slider_width=0.03, tube_width=0.006, color="black")

        return self.plotter


    def plotlogs(self, state="start", ptype="neutron", pselect=None, color=None):
        '''
        state = "start", "end", or "hit"
        '''
        self.plotter = self.initplotter() if self.plotter is None else self.plotter
        
        if not hasattr(self.plotter, "logsactor"):
            self.plotter.logsactor = {}
             
        st = '' if "hit" in state else state
        ltype = 'h' if "hit" in state else 'e'

        pos = 0
        if (ptype == 'n' or ptype == "neutron"):
            key = 'n'
        elif (ptype == 'p' or ptype == "proton"):
            key = 'p'
            pos = 60
        elif (ptype == 'e' or ptype == "electron"):
            key = 'e'
            pos = 90
            

        if state == "start":
            color = "black" if color is None else color
        elif state == "end":
            color = "red" if color is None  else color
            pos += 30
        else:
            color = "pink" if color is None  else color
            pos += 60

        df_log = da.df[key+ltype]
        
        if pselect is not None:
            df_log = df_log[df_log['particle'].isin(pselect)]
            
        positions = df_log[['x'+st, 'y'+st, 'z'+st]].values
        cloud = pv.PolyData(1000 * positions)

        # cloud["colors"] = np.ones((cloud.n_faces, 3)) * 100
        
        actor = self.plotter.add_mesh(1000 * positions, color=color, render_points_as_spheres=True, point_size=6, lighting=lighting)
        self.plotter.logsactor[key+ltype+state] = actor
        
        callback = self.SetVisibilityCallback(actor)
        self.plotter.add_checkbox_button_widget(
            callback,
            value=True,
            position=(300 + pos, 1),
            size=30,
            border_size=1,
            color_on=color,
            color_off='grey',
            background_color='grey',
        )

        
        return self.plotter


    def speedanime(self, speed):
        self.plotter.dt = self.plotter.idt * np.sign(speed) * speed**2


    
    def animate(self, ti=0, tf=None, dt=0.005, fps=10, spath="/home/sly/Work/Physics/Neutrons/tSPECT/PENtrack/plots/animations/penanalyse/", pselect=None, minp=0, maxp=None, displaypID=False):
        """
        Main function for animation
        """

        pcolor = "darkorange"
        self.plotter = self.initplotter() if self.plotter is None else self.plotter
        self.plotter.savegif = False
        
        df_nts = self.df['nt'] #.sort_values(['particle', 't'])
        df_ne = self.df['ne']
        maxp = self.simcount if maxp is None else maxp
        
        if pselect is not None:
            df_nts = df_nts[df_nts['particle'].isin(pselect)]
            df_ne = df_ne[df_ne['particle'].isin(pselect)]
        # else:
        df_nts = df_nts[(df_nts['particle'] >= minp) & (df_nts['particle'] <= maxp)]
        df_ne = df_ne[(df_ne['particle'] >= minp) & (df_ne['particle'] <= maxp)]
            
        df_nts = df_nts.sort_values(['particle', 't'])
        df_ne = df_ne.sort_values(['particle'])

        if tf is not None:
            df_nts = df_nts[df_nts['t'] < tf]
        else:
            tf = self.simtime
        
        # self.plotter.enable_depth_peeling()
        cloud = pv.PolyData(1000 * df_ne[['xstart', 'ystart', 'zstart']].values)

        rgba = np.array(pv.Color(pcolor).int_rgba)[None] * np.ones((cloud.n_faces, 4))
        rgba[:, 0] = 255//2 + 75*df_ne['polstart'].values
        rgba[:, 1] = 255//2 - 75*df_ne['polstart'].values
        rgba[:, 2] = 120
        rgba[:, -1] = 1
        cloud["rgba"] = rgba

        self.plotter.pactor = self.plotter.add_points(cloud, render_points_as_spheres=True, point_size=8, scalars="rgba", rgba=True, lighting=lighting)

        if displaypID:
            cloud["pID"] = [f"{int(i)}" for i in df_ne['particle']]

            self.plotter.labelpactor = self.plotter.add_point_labels(cloud, "pID", font_size=20, render_points_as_spheres=True, point_size=6)

            callback = self.SetViqsibilityCallback(self.plotter.labelpactor)
            self.plotter.add_checkbox_button_widget(
                callback,
                value=True,
                position=(550, 40),
                size=30,
                border_size=1,
                color_on="black",
                color_off='grey',
                background_color='grey',
            )


        callback = self.RecordGifCallback(self.plotter, self.dfile, spath, fps=5)
        self.plotter.add_checkbox_button_widget(
            callback,
            value=False,
            position=(800, 1),
            size=30,
            border_size=3,
            color_on="maroon",
            color_off='grey',
            background_color='grey',
        )
        
        
        callback = self.ScreenshotCallback(self.plotter, self.dfile, spath)
        self.plotter.add_checkbox_button_widget(
            callback,
            value=False,
            position=(800, 40),
            size=30,
            border_size=1,
            color_on="gray",
            color_off='lime',
            background_color='grey',
        )
        

        interpolation_functions = {}
        for particle, data in df_nts.groupby('particle'):
            positions = data[['x', 'y', 'z', 'polarisation']].values
            times = data['t'].values
            interpolation_functions[particle] = scipy.interpolate.interp1d(times, positions, axis=0, bounds_error=False, fill_value=(data[['x', 'y', 'z', 'polarisation']].iloc[0], data[['x', 'y', 'z', 'polarisation']].iloc[-1]))

        self.plotter.time_stamp = ti
        self.plotter.idt = dt
        self.plotter.timetxt = self.plotter.add_text("t=%g [s]"%self.plotter.time_stamp, position='lower_right', font_size=18)
        
        def update_animation():
            if ((self.plotter.dt < 0) and (self.plotter.time_stamp > ti)) or ((self.plotter.dt > 0) and (self.plotter.time_stamp < tf)):
                self.plotter.time_stamp += self.plotter.dt
                self.plotter.timetxt.SetText(1, "t=%g [s]"%self.plotter.time_stamp)
                # self.plotter.update_text("t=%g [s]"%self.plotter.time_stamp, self.plotter.timetxt)
                
                interpolated_vals = []
                for particle, interpolation_func in interpolation_functions.items():
                    interpolated_vals.append(interpolation_func(self.plotter.time_stamp))

                interpolated_vals = np.array(interpolated_vals)
                cloud.points = 1000 * interpolated_vals[:, :3]

                # rgba = np.ones((len(interpolated_vals), 4))# np.array(pv.Color(pcolor).int_rgba)[None] * np.ones((cloud.n_faces, 4))
                cloud["rgba"][:, 0] = 255//2 + 75*interpolated_vals[:, 3]
                cloud["rgba"][:, 1] = 255//2 - 75*interpolated_vals[:, 3]
                cloud["rgba"][:, 3] = 1.0 * np.array(self.plotter.time_stamp < df_ne['tend']).astype(int)
                
                # self.plotter.render()
                self.plotter.update()

                if self.plotter.savegif:
                    self.plotter.write_frame()
            
        self.plotter.add_slider_widget(self.speedanime, (-4, 4), value=1, pointa=(0.5, 0.03), pointb=(0.8, 0.03), slider_width=0.02, tube_width=0.005, color="white")

        self.plotter.add_callback(update_animation, int(1000/fps), None)
        return self.plotter



# name of datafile
# dfile = "000000000105"
dfile = "000000000010"

# instantiate data object
da = data(dfile)

# plots stl surface and magnetic source
# pl = da.plotmag(opacity=0.01)
pl = da.plotstl(opacity=0.01)

# df_ne = da.df['ne']

# pselect = df_ne[df_ne['xend'] > df_ne['xstart'] + 0.1]['particle']
pselect = None #np.arange(1, da.df['ne'])
# pselect = np.arange(1, 8)

# plots neutrons start, end, and hits points.
pl = da.plotlogs(ptype="n", state="start", pselect=pselect, color="lightgreen")
pl = da.plotlogs(ptype="n", state="end", pselect=pselect, color="deepskyblue")
pl = da.plotlogs(ptype="n", state="hit", pselect=pselect, color="darkorchid")

# # play animation
pl = da.animate(dt=0.01, fps=10, pselect=pselect, minp=0, maxp=None)

