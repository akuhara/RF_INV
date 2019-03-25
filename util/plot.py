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
                       'k_min', 'k_max', 'z_min',  'z_max', 'dvs_prior', \
                       'dvp_prior', 'sig_min', 'sig_max', 'dev_z', \
                       'dev_dvs', 'dev_dvp', 'dev_sig', 'nbin_z', \
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
                        
            elif (param_names[i] == "t_start" or param_names[i] == "k_min" or \
                  param_names[i] == "z_min" or param_names[i] == "sig_min" or \
                  param_names[i] == "amp_min" or param_names[i] == "vp_min" or \
                  param_names[i] == "vs_min"):
                param[param_names[i]] = item[0]
                i += 1
                param[param_names[i]] = item[1]
                
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
    def plot_num_interface(self, ax):
        
        param = self.param
        xlabel = "# of layer interfaces"
        ylabel = "Posteiror probability"
        rslt_file = param["outdir"] + "/" + "num_interface.ppd"

        # read file
        df = pd.read_csv(rslt_file, delim_whitespace=True, \
                         header=None, names=(xlabel, ylabel))
        
        # plot
        df.plot(x=xlabel, y=ylabel, ax=ax, kind="area", legend=None)
        plt.ylabel(ylabel)
        
    #------------------------------------------------------------------    
    def plot_sigma(self, ax, trace_id):
        param = self.param
        xlabel = "Standard deviation of data noise"
        ylabel = "Posteiror probability"
        zlabel = "Trace ID"
        rslt_file = param["outdir"] + "/" + "sigma.ppd"
        df = pd.read_csv(rslt_file, delim_whitespace=True, \
                         header=None, names=(xlabel, ylabel, zlabel))
        df[df[zlabel] == trace_id].plot(x=xlabel, y=ylabel, ax=ax, \
                                        kind="area", legend=None)
        plt.ylabel(ylabel)
        
    #------------------------------------------------------------------  
    
    def plot_syn_trace(self, ax, trace_id):
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

        plt.pcolor(x, y, data, cmap='hot_r')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(zlabel)

        # plot observation
        b, delta, t, data = self.read_sac(trace_id-1)
        plt.plot(t, data, color="black")
        t_start = float(param["t_start"])
        t_end = float(param["t_end"])
        plt.xlim([t_start, t_end])
     #------------------------------------------------------------------    

    def plot_v_z(self, ax, vtype):
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
        ax = plt.pcolormesh(x, y, data, cmap='hot_r', vmin=0.0, vmax=0.2)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.ylim([z_max + del_z, 0])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(zlabel)

        # plot reference velocity 
        xlabel2 = "Depth (km)"
        ylabel2 = "P wave velocity (km/s)"
        zlabel2 = "S wave velocity (km/s)"
        df = pd.read_csv(param["vel_file"], delim_whitespace=True, \
                         header=None, names=(xlabel2, ylabel2, zlabel2))
        # reference velocity 
        if (vtype == "vs"):
            plt.plot(df[zlabel2], df[xlabel2], color="black")
        elif (vtype == "vp"):
            plt.plot(df[ylabel2], df[xlabel2], color="black")
            
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
        
        if (vtype == "vs"):
            plt.plot(df_up[zlabel2], df_up[xlabel2], color="black", \
                     linestyle="dashed")
            plt.plot(df_low[zlabel2], df_low[xlabel2], color="black", \
                     linestyle="dashed")
        elif (vtype == "vp"):
            plt.plot(df_up[ylabel2], df_up[xlabel2], color="black", \
                     linestyle="dashed")
            plt.plot(df_low[ylabel2], df_low[xlabel2], color="black", \
                     linestyle="dashed")
            
        
    #------------------------------------------------------------------    
    def make_figure(self, trace_id):
        param = self.param
        sns.set()
        sns.set_style('ticks')
        plt.figure(figsize=(9, 40))
        
        # Adjust intervals between suplots
        plt.subplots_adjust(wspace=0.6, hspace=0.6)

        ax = plt.subplot2grid((4, 2), (0, 0))
        self.plot_num_interface(ax)
        
        if (float(param["sig_max"]) - float(param["sig_min"]) > 1.0e-5):
            ax = plt.subplot2grid((4, 2), (0, 1))
            self.plot_sigma(ax, trace_id)
            
        ax = plt.subplot2grid((4, 2), (1, 0), colspan=2)
        self.plot_syn_trace(ax, trace_id)

        ax = plt.subplot2grid((4, 2), (2, 0), rowspan=2)
        self.plot_v_z(ax, vtype='vs')

        if (param["vp_mode"] == 1):
            ax = plt.subplot2grid((4, 2), (2, 1), rowspan=2)
            self.plot_v_z(ax, vtype='vp')
        
        plt.show()
        
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

    for itrc in range (1, int(rslt.param["ntrc"]) + 1):
        rslt.make_figure(itrc)

    
    
