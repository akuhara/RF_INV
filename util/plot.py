import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import struct

class InvRslt:
    def __init__(self, param_file):
        self.param_file = param_file

    #------------------------------------------------------------------
    def read_param_file(self):
        # Parameter list
        param_names = ('outdir', 'nburn', 'niter', 'ncorr', 'nchains', \
                       'ncool', 't_high', 'iseed', 'ntrc', 'rayp', 'a_gus', \
                       'ipha', 'nfft', 'obs_file', 't_start', 't_end', \
                       'conv_mode', 'sdep', 'vel_file', 'vp_mode', \
                       'k_min', 'k_max', 'z_min',  'z_max', 'h_min', 'prior_mode', \
                       'dvs_prior', 'dvp_prior', 'sig_min', 'sig_max', \
                       'dev_z', 'dev_dvs', 'dev_dvp', 'dev_sig', 'nbin_z', \
                       'nbin_vs', 'nbin_vp', 'nbin_sig', 'nbin_amp', \
                       'amp_min', 'amp_max', 'vp_min', 'vp_max', 'vs_min', \
                       'vs_max')
        
        # Read prameter file
        param = {}
        f = open(self.param_file, 'r')
        line = f.readline()
        i = 0
        while line:
            item = line.split()
            item[0] = item[0].replace("'", "")
            item[0] = item[0].replace('"', '')
            # Skip comment line
            if (re.match(r"^#" ,item[0])):
                line = f.readline()
                continue
                
            # Get parameters
            if (param_names[i] == "rayp" or param_names[i] == "a_gus" or \
                param_names[i] == "ipha" or param_names[i] == "obs_file"):
                args = []
                for j in range(0, int(param["ntrc"])):
                    #--- chomp & split for later traces ---
                    if (j > 0): 
                        line = f.readline()
                        item = line.split()
                    #---------------------------------------

                    item[0] = item[0].replace("'", "")
                    item[0] = item[0].replace('"', '')
                    args.append(item[0])

                param[param_names[i]] = args
                        
            elif param_names[i] == "t_start" or \
                 param_names[i] == "k_min" or \
                 param_names[i] == "z_min" or \
                 param_names[i] == "amp_min" or \
                 param_names[i] == "vp_min" or \
                 param_names[i] == "vs_min":
                param[param_names[i]] = item[0]
                i += 1
                param[param_names[i]] = item[1]
                
            elif param_names[i] == "sig_min":
                args0 = []
                args1 = []
                for j in range(0, int(param["ntrc"])):
                    #--- chomp & split for later traces ---
                    if (j > 0): 
                        line = f.readline()
                        item = line.split()
                    #---------------------------------------
                    
                    args0.append(item[0])
                    args1.append(item[1])

                param[param_names[i]] = args0
                i += 1
                param[param_names[i]] = args1
                
            else:
                param[param_names[i]] = item[0]
                
            i += 1
            line = f.readline()
            f.close
            
        self.param = param

    #------------------------------------------------------------------  
    def read_sac(self, trace_id):

        param = self.param
        obs_file = param["obs_file"][trace_id]
        endian = "<"

        # read entire part of SAC file
        f = open(obs_file, 'rb')
        sac_data = f.read()

        # check endian
        header = struct.unpack_from(endian + 'i', sac_data, offset=4*76)
        if (header[0] != 6):
            endian = ">"
            
        # get headers
        header = struct.unpack_from(endian + 'f', sac_data, offset=0)
        delta = header[0]
        header = struct.unpack_from(endian + 'f', sac_data, offset=4*5)
        b = header[0]
        header = struct.unpack_from(endian + 'i', sac_data, offset=4*79)
        npts = header[0]
        
        # get data
        data = struct.unpack_from(endian + 'f' * npts, sac_data, offset=4*158)
        
        # make time index
        t = np.arange(b, b + npts * delta, delta)
        
        return b, delta, t, data
        

    #------------------------------------------------------------------    
    def plot_num_interface(self, fig, ax):
        
        param = self.param
        xlabel = "# of layer interfaces"
        ylabel = "Posteiror probability"
        rslt_file = param["outdir"] + "/" + "num_interface.ppd"

        # read file
        df = pd.read_csv(rslt_file, delim_whitespace=True, \
                         header=None, names=(xlabel, ylabel))
        
        # plot
        df.plot(x=xlabel, y=ylabel, ax=ax, kind="area", legend=None)
        ax.set_ylabel(ylabel)
        
    #------------------------------------------------------------------    
    def plot_sigma(self, fig, ax, trace_id):
        param = self.param
        xlabel = "Standard deviation of data noise"
        ylabel = "Posteiror probability"
        zlabel = "Trace ID"
        rslt_file = param["outdir"] + "/" + "sigma.ppd"
        df = pd.read_csv(rslt_file, delim_whitespace=True, \
                         header=None, names=(xlabel, ylabel, zlabel))
        df[df[zlabel] == trace_id].plot(x=xlabel, y=ylabel, ax=ax, \
                                        kind="area", legend=None)
        ax.set_ylabel(ylabel)
        
    #------------------------------------------------------------------  
    
    def plot_syn_trace(self, fig, ax, trace_id):
        param = self.param
        xlabel = "Time after P (s)"
        ylabel = "Amplitude"
        zlabel = "Posterior probability"
        wlabel = "Trace ID"
        rslt_file = param["outdir"] + "/" + "syn_trace.ppd"
        df = pd.read_csv(rslt_file, delim_whitespace=True, \
                         header=None, names=(xlabel, ylabel, zlabel, wlabel))
        
        amp_min = float(param["amp_min"])
        amp_max = float(param["amp_max"])
        nbin_amp = int(param["nbin_amp"])
        t_start = float(param["t_start"])
        t_end = float(param["t_end"])
        del_amp = (amp_max - amp_min) / nbin_amp
        
        y, x = np.mgrid[slice(amp_min, amp_max + del_amp, del_amp), \
                        slice(t_start - 0.5 * 0.05, t_end + 0.5 * 0.05, 0.05)]
        
        data = df[df[wlabel] == trace_id].pivot(ylabel, xlabel, zlabel)

        mappable = ax.pcolor(x, y, data, cmap='hot_r')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        cbar = fig.colorbar(mappable, ax=ax)
        cbar.ax.set_ylabel(zlabel)

        # plot observation
        b, delta, t, data = self.read_sac(trace_id-1)
        ax.plot(t, data, color="black")
        t_start = float(param["t_start"])
        t_end = float(param["t_end"])
        ax.set_xlim([t_start, t_end])
        
     #------------------------------------------------------------------    

    def plot_v_z(self, fig, ax, vtype):
        param = self.param

        # plot posterior probability
        if vtype == "vs":
            xlabel = "S wave velocity (km/s)"
            rslt_file = param["outdir"] + "/" + "vs_z.ppd"
            v_min = float(param["vs_min"])
            v_max = float(param["vs_max"])
            nbin_v = int(param["nbin_vs"])
        elif vtype == "vp":
            xlabel = "P wave velocity (km/s)"
            rslt_file = param["outdir"] + "/" + "vp_z.ppd"
            v_min = float(param["vp_min"])
            v_max = float(param["vp_max"])
            nbin_v = int(param["nbin_vp"])
            
        ylabel = "Depth (km)"
        zlabel = "Posterior probability"
        del_v = (v_max - v_min) / nbin_v
        z_min = 0.0
        z_max = float(param["z_max"])
        nbin_z = int(param["nbin_z"])
        del_z = (z_max - z_min) / nbin_z
        
        
        df = pd.read_csv(rslt_file, delim_whitespace=True, \
                         header=None, names=(xlabel, ylabel, zlabel))
        y, x = np.mgrid[slice(z_min, z_max + del_z, del_z), \
                        slice(v_min, v_max + del_v, del_v)]
        data = df.pivot(ylabel, xlabel, zlabel)
        mappable = ax.pcolormesh(x, y, data, cmap='hot_r', \
                                 vmin=0.0, vmax=0.2)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_ylim([z_max + del_z, 0])

        if ((vtype == "vs" and int(param["vp_mode"]) == 0) or \
            vtype == "vp"):
            cbar = fig.colorbar(mappable, ax=ax)
            cbar.ax.set_ylabel(zlabel)

        # plot reference velocity 
        xlabel2 = "Depth (km)"
        ylabel2 = "P wave velocity (km/s)"
        zlabel2 = "S wave velocity (km/s)"
        df = pd.read_csv(param["vel_file"], delim_whitespace=True, \
                         header=None, names=(xlabel2, ylabel2, zlabel2))

        # for legend 
        lines = []
        labels = []
        
        # reference velocity 
        if (vtype == "vs"):
            line, = ax.plot(df[zlabel2], df[xlabel2], color="black")
        elif (vtype == "vp"):
            line, = ax.plot(df[ylabel2], df[xlabel2], color="black")
        lines.append(line)
        labels.append("Reference model")
            
        # plot 1-sigma bound

        # get 1-sigma
        if (vtype == "vs"):
            sigma = float(param["dvs_prior"])
            v_min = float(param["vs_min"])
            v_max = float(param["vs_max"])
        elif (vtype == "vp"):
            sigma = float(param["dvp_prior"])
            v_min = float(param["vp_min"])
            v_max = float(param["vp_max"])
        
        # make lower and upper bound for below seafloor depth
        z_min = float(param["z_min"])
        df_up = df + np.array([0, sigma, sigma])
        df_low = df - np.array([0, sigma, sigma])
        df_up = df_up[df_up[xlabel2] >= z_min]
        df_low = df_low[df_low[xlabel2] >= z_min]
        
        #if (vtype == "vs"):
        #    ax.plot(df_up[zlabel2], df_up[xlabel2], color="black", \
        #             linestyle="dashed", label=None)
        #    ax.plot(df_low[zlabel2], df_low[xlabel2], color="black", \
        #             linestyle="dashed", label=None)
        #elif (vtype == "vp"):
        #    ax.plot(df_up[ylabel2], df_up[xlabel2], color="black", \
        #             linestyle="dashed")
        #    ax.plot(df_low[ylabel2], df_low[xlabel2], color="black", \
        #             linestyle="dashed")

        # Meam model
        ylabel3 = "Depth (km)"
        if (vtype == "vs"):
            mean_file = param["outdir"] + "/" + "vs_z.mean"
            xlabel3 = "S wave velocity (km/s)"
        elif (vtype == "vp"):
            mean_file = param["outdir"] + "/" + "vp_z.mean"
            xlabel3 = "P wave velocity (km/s)"
        
        df = pd.read_csv(mean_file, delim_whitespace=True, \
                         header=None, names=(xlabel3, ylabel3))
        line, = ax.plot(df[xlabel3], df[ylabel3], color="blue")
        lines.append(line)
        labels.append("Mean model")
        ax.legend(lines, labels)
        
            
    #------------------------------------------------------------------    

    def plot_likelihood(self, fig, ax):
        param = self.param
        xlabel = "Iteration number"
        ylabel = "Log-likelihood"
        rslt_file = param["outdir"] + "/" + "likelihood"
        x1 = int(param["nburn"]) + 1
        x2 = x1 + int(param["niter"]) - 1
        
        df = pd.read_csv(rslt_file, delim_whitespace=True, \
                         header=None, names=(xlabel, ylabel))
        df.plot(x=xlabel, y=ylabel, ax=ax, kind="line", legend=None)
        #ax.set_xscale("log")
        #ax.set_yscale("log")
        ax.set_ylabel(ylabel)
        ax.set_ylim(400, )
        ax.axvspan(x1, x2, facecolor="orange")
        
    
    #------------------------------------------------------------------    
    def make_figure(self, trace_id):
        
        grid_geom = (5, 2)

        param = self.param
        sns.set()
        sns.set_style('ticks')
        fig = plt.figure(figsize=(9, 12))
        
        # Adjust intervals between suplots
        fig.subplots_adjust(wspace=0.6, hspace=0.6)

        ax = plt.subplot2grid(grid_geom, (0, 0), fig=fig)
        self.plot_num_interface(fig, ax)
        
        if (float(param["sig_max"][trace_id-1]) \
            - float(param["sig_min"][trace_id-1]) > 1.0e-5):
            ax = plt.subplot2grid(grid_geom, (0, 1), fig=fig)
            self.plot_sigma(fig, ax, trace_id)
            
        ax = plt.subplot2grid(grid_geom, (1, 0), fig=fig)
        self.plot_likelihood(fig, ax)

        ax = plt.subplot2grid(grid_geom, (2, 0), colspan=1, fig=fig)
        self.plot_syn_trace(fig, ax, trace_id)
        
        ax = plt.subplot2grid(grid_geom, (3, 0), rowspan=2, fig=fig)
        self.plot_v_z(fig, ax, vtype='vs')
        
        if (int(param["vp_mode"]) == 1):
            ax = plt.subplot2grid(grid_geom, (3, 1), rowspan=2, fig=fig)
            self.plot_v_z(fig, ax, vtype='vp')
            
        png_file = param["outdir"] + "/" + "plot" +  \
                   str(trace_id).zfill(2) + ".ps"
        print(png_file)
        fig.savefig(png_file)
        
  

#=======================================================================

if __name__ == "__main__":
    
    # Get argument
    args = sys.argv
    if len(args) != 2:
        print("USAGE: plot.py [parameter file]")
        sys.exit()
        
    param_file = args[1]
    rslt = InvRslt(param_file)
    rslt.read_param_file()
    rslt.read_sac(0)
    print(rslt.param)
    for itrc in range (1, int(rslt.param["ntrc"]) + 1):
        rslt.make_figure(itrc)

    
    
