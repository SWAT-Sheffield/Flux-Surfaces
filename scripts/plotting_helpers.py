# -*- coding: utf-8 -*-
"""
A module of handy plotting functions and shortcuts for Time-Distance diagrams
"""
from __future__ import absolute_import, division, print_function

import os
import glob
import itertools

import numpy as np

from scipy.interpolate import interp1d

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
from mpl_toolkits.axes_grid1 import make_axes_locatable

__all__ = ['glob_files', 'calc_avgs', 'single_plot', 'triple_plot', 'triple_triple_plot',
           'get_single_percentage_flux', 'get_single_velocity', 'get_single_bpert',
           'get_triple', 'get_xy', 'save_fig', 'get_pgf_include', 'make_fig_env',
           'get_fig_env', 'make_sub_figs', 'add_phase_speeds', 'get_phase_speeds']

xxlim = -1
def betaswap(b,pos):
    return "$%3.2f$"%(b)

def glob_files(cfg, tube_r, search):
    """
    Search in the configs data directory for the following pattern:

    {data_dir}+{tube_r}+{search}

    Parameters
    ----------
    tube_r: string
        The tube radius string to search for

    search: string
        The rest of the search string

    Returns
    -------
    files: list
        A sorted files list
    """
    files = glob.glob(os.path.join(cfg.data_dir,tube_r,search))
    files.sort()
    return files

def calc_avgs(cfg, tube_r, periods, amps, runtime=600.):
    """
    Calculate the average values of Fpar, Fperp and Fphi over all time for an
    integer number of periods.

    Parameters
    ----------
    cfg : scripts.sacconfig.SACConfig instance
        The repo config to use

    tube_r : string
         The tube radius to use

    periods : ndarray
        List of all periods

    amps : list of strings
        List of all corresponding amplitudes

    runtime : float
        Total runtime of the simulation (to calculate integer periods)

    Returns
    -------
    AvgsP : ndarray
        A 3 x len(periods) array of average fluxes (Par, Perp, Phi)
    """
    int_periods = np.floor(runtime/periods)*periods
    AvgsP = np.zeros([3,len(periods)])
    for i, period, amp in itertools.izip(range(len(periods)), periods, amps):
        cfg.period = period
        cfg.amp = amp

        times = np.load(os.path.join(cfg.data_dir, 'Times_{}.npy'.format(cfg.get_identifier())))
        max_index = np.argmin(np.abs(int_periods[i] - times))

        Fpar, Fperp, Fphi = map(np.load,glob_files(cfg, tube_r, 'LineFlux*Fp*npy'))
        #Fpar, Fperp, Fphi = map(np.load,glob_files(cfg, tube_r, '*vp*npy'))
        Fpar[np.abs(Fpar)<1e-5], Fperp[np.abs(Fperp)<1e-5], Fphi[np.abs(Fphi)<1e-5] = 0., 0., 0.
        Fpar, Fperp, Fphi = Fpar[:max_index,:], Fperp[:max_index,:], Fphi[:max_index,:]

        Ftot2 = (Fpar**2 + Fperp**2 + Fphi**2)
        Fpar2, Fperp2, Fphi2 = np.array([Fpar, Fperp, Fphi])**2
        FparP, FperpP, FphiP = (Fpar2/Ftot2)*100, (Fperp2/Ftot2)*100, (Fphi2/Ftot2)*100

        FparP = np.mean(np.ma.masked_array(FparP,np.isnan(FparP)))
        FperpP = np.mean(np.ma.masked_array(FperpP,np.isnan(FperpP)))
        FphiP = np.mean(np.ma.masked_array(FphiP,np.isnan(FphiP)))

        AvgsP[:, i] = FparP, FperpP, FphiP

    return AvgsP

def single_plot(data, x, y, axes=None, beta=None, cbar_label='',
                cmap=plt.get_cmap('RdBu'), vmin=None, vmax=None,
                phase_speeds=True, **kwargs):
    """
    Plot a single frame Time-Distance Diagram on physical axes.

    This function uses mpl NonUniformImage to plot a image using x and y arrays,
    it will also optionally over plot in contours beta lines.

    Parameters
    ----------
    data: np.ndarray
        The 2D image to plot
    x: np.ndarray
        The x coordinates
    y: np.ndarray
        The y coordinates
    axes: matplotlib axes instance [*optional*]
        The axes to plot the data on, if None, use plt.gca().
    beta: np.ndarray [*optional*]
        The array to contour over the top, default to none.
    cbar_label: string [*optional*]
        The title label for the colour bar, default to none.
    cmap: A matplotlib colour map instance [*optional*]
        The colourmap to use, default to 'RdBu'
    vmin: float [*optional*]
        The min scaling for the image, default to the image limits.
    vmax: float [*optional*]
        The max scaling for the image, default to the image limits.
    phase_speeds : bool
        Add phase speed lines to the plot

    Returns
    -------
    None
    """
    if axes is None:
        axes = plt.gca()

    x = x[:xxlim]
    data = data[:,:xxlim]

    im = NonUniformImage(axes,interpolation='nearest',
                         extent=[x.min(),x.max(),y.min(),y.max()],rasterized=False, **kwargs)
    im.set_cmap(cmap)
    if vmin is None and vmax is None:
        lim = np.max([np.nanmax(data),
                  np.abs(np.nanmin(data))])
        im.set_clim(vmax=lim,vmin=-lim)
    else:
        im.set_clim(vmax=vmax,vmin=vmin)
    im.set_data(x,y,data)
    im.set_interpolation('nearest')

    axes.images.append(im)
    axes.set_xlim(x.min(),x.max())
    axes.set_ylim(y.min(),y.max())

    cax0 = make_axes_locatable(axes).append_axes("right", size="5%", pad=0.05)
    cbar0 = plt.colorbar(im, cax=cax0, ticks=mpl.ticker.MaxNLocator(7))
    cbar0.set_label(cbar_label)
    cbar0.solids.set_edgecolor("face")
    kwergs = {'levels': [1., 1/3., 1/5., 1/10., 1/20.]}
    kwergs.update(kwargs)

    if beta is not None:
        ct = axes.contour(x,y,beta[:,:xxlim],colors=['k'], **kwergs)
        plt.clabel(ct,fontsize=14,inline_spacing=3, manual=False,
                   fmt=mpl.ticker.FuncFormatter(betaswap))

    axes.set_xlabel("Time [s]")
    axes.set_ylabel("Height [Mm]")

def triple_plot(ax, x, y, par_line, perp_line, phi_line, beta_line=None,
                par_label='', perp_label='', phi_label='', title='', **kwargs):
    """
    Plot the three components of a T-D diagram.

    Parameters
    ----------
    ax: np.ndarray
        An array of axes, 1D.
    x: np.ndarray
        The x coordinates.
    y: np.ndarray
        The y coordinates.
    par_line: np.ndarray
        The first T-D array.
    perp_line: np.ndarray
        The second T-D array.
    phi_line: np.ndarray
        The third T-D array.
    beta_line: np.ndarray [*optional*]
        The beta array to over plot.
    par_label: string [*optional*]
        The colour bar label for the fist array.
    perp_label: string [*optional*]
        The colour bar label for the second array.
    phi_label: string [*optional*]
        The colour bar label for the third array.
    title: string [*optional*]
        The title for above the top image.
    """
    ax[0].set_title(title)

    single_plot(par_line.T[::-1,:], x, y, axes=ax[0],
                cbar_label=par_label, beta=beta_line, **kwargs)

    single_plot(perp_line.T[::-1,:], x, y, axes=ax[1],
                cbar_label=perp_label, beta=beta_line, **kwargs)

    single_plot(phi_line.T[::-1,:], x, y, axes=ax[2],
                cbar_label=phi_label, beta=beta_line, **kwargs)

def triple_triple_plot(title,
                x_r10, y_r10, beta_line_r10, par_line_r10, perp_line_r10, phi_line_r10,
                x_r30, y_r30, beta_line_r30, par_line_r30, perp_line_r30, phi_line_r30,
                x_r60, y_r60, beta_line_r60, par_line_r60, perp_line_r60, phi_line_r60,
                par_label = r'', perp_label = r'', phi_label = r'',**kwargs):
    """
    This function plots all 3 components for all 3 tube radii.

    This function creates and returns a figure, with 9 subplots.

    Parameters
    ----------
    title: string
        The whole figure title
    x_r10: np.ndarray
        x coordinate array for first column.
    y_r10: np.ndarray
        x coordinate array for first column.
    beta_line_r10: np.ndarray
        x coordinate array for first column.
    par_line_r10: np.ndarray
        x coordinate aray for first column.
    perp_line_r10: np.ndarray
        x coordinate array for first column.
    phi_line_r10: np.ndarray
        x coordinate array for first column.
    x_r30: np.ndarray
        x coordinate array for second column.
    y_r30: np.ndarray
        x coordinate array for second column.
    beta_line_r30: np.ndarray
        x coordinate array for second column.
    par_line_r30: np.ndarray
        x coordinate array for second column.
    perp_line_r30: np.ndarray
        x coordinate array for second column.
    phi_line_r30: np.ndarray
        x coordinate array for second column.
    x_r60: np.ndarray
        x coordinate array for third column.
    y_r60: np.ndarray
        x coordinate array for third column.
    beta_line_r60: np.ndarray
        x coordinate array for third column.
    par_line_r60: np.ndarray
        x coordinate array for third column.
    perp_line_r60: np.ndarray
        x coordinate array for third column.
    phi_line_r60: np.ndarray
        x coordinate array for third column.
    par_label: string
        Label for parallel color bar
    perp_label: string
        Label for perp color bar
    phi_label: string
        Label for phi color bar

    Returns
    -------
    figure: matplotlib.figure.Figure
        The figure containing all the plots
    """

    fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(15,6))
    fig.suptitle(title, y=1.05)

    #r10
    triple_plot(ax[:,0], x_r10, y_r10,
                par_line_r10, perp_line_r10, phi_line_r10,
                beta_line=beta_line_r10, par_label=par_label,
                perp_label=perp_label, phi_label=phi_label,
                title="Flux Surface Radius $= 158$ km", **kwargs)

    #r30
    triple_plot(ax[:,1], x_r30, y_r30,
                par_line_r30, perp_line_r30, phi_line_r30,
                beta_line=beta_line_r30, par_label=par_label,
                perp_label=perp_label, phi_label=phi_label,
                title="Flux Surface Radius $= 468$ km", **kwargs)
    #r60
    triple_plot(ax[:,2], x_r60, y_r60,
                par_line_r60, perp_line_r60, phi_line_r60,
                beta_line=beta_line_r60, par_label=par_label,
                perp_label=perp_label, phi_label=phi_label,
                title="Flux Surface Radius $= 936$ km", **kwargs)

    plt.tight_layout()

    return fig, ax

def add_phase_speeds(ax, x, y, va_line, cs_line, x_shift=60, color='k'):
    """
    Plot the for phase speed lines on an axes
    """
    delta_x = np.zeros(y.shape)
    delta_x[1:] = y[1:] - y[:-1]

    delta_t_va = delta_x / va_line[:,0]
    delta_t_cs = delta_x / cs_line[:,0]
    delta_t_vf = delta_x / np.sqrt(cs_line[:,0]**2 + va_line[:,0]**2)
    delta_t_vs = delta_x / np.sqrt(cs_line[:,0]**-2 + va_line[:,0]**-2)**-1

    t_va = np.cumsum(delta_t_va) + x_shift
    t_cs = np.cumsum(delta_t_cs) + x_shift
    t_vf = np.cumsum(delta_t_vf) + x_shift
    t_vs = np.cumsum(delta_t_vs) + x_shift

    ax.plot(t_va, y, label=r"$V_A$", linewidth=2, linestyle=':', color=color)#b
    ax.plot(t_cs, y, label=r"$C_s$", linewidth=2, linestyle='--', color=color)#g
    ax.plot(t_vf, y, label=r"$V_f$", linewidth=2, linestyle='-.', color=color)#r
    ax.plot(t_vs, y, label=r"$V_s$", linewidth=2, linestyle='-', color=color)#c

def get_phase_speeds(cfg, tube_r):
    """
    Read in the data for phase speed plotting.
    """
    x, y = get_xy(cfg, tube_r)
    va_line = np.load(glob_files(cfg, tube_r, '*va.npy')[0]).T * 1e-6
    cs_line = np.load(glob_files(cfg, tube_r, '*cs.npy')[0]).T * 1e-6
    x_shift = cfg.period / 4.

    return {'x':x, 'y':y, 'va_line':va_line,
            'cs_line':cs_line, 'x_shift': x_shift}

def get_xy(cfg, tube_r):
    """
    read in from file and interpolate the x and y coords
    """
    #There is a bug in the gdf files that mean the x array is wrong:
    height_Mm = np.linspace(0.03664122,1.5877863,128)
    f = interp1d(np.linspace(0,128,128),height_Mm)
    all_spoints = np.load(glob_files(cfg, tube_r,'*points*npy')[0])[:,::-1,:]
    all_times = np.load(glob_files(cfg, tube_r,'*times*npy')[0])
    y = f(all_spoints[0,:,2])
    x = all_times[:-1]
    return x,y

def get_single_velocity(cfg, tube_r, search='*_vp*npy', beta=True,
                        triple_triple=False):
    """
    Read in the args for single_plot from files
    """
    kwargs = {}
    vpar, vperp, vphi = map(np.load,glob_files(cfg, tube_r, search))
    if triple_triple:
        tube_r2 = '_'+tube_r
    else:
        tube_r2 =''

    [kwargs['par_line%s'%tube_r2],
     kwargs['perp_line%s'%tube_r2],
     kwargs['phi_line%s'%tube_r2]] = vpar, vperp, vphi


    kwargs['x%s'%tube_r2],kwargs['y%s'%tube_r2] = get_xy(cfg, tube_r)
    if beta:
        kwargs['beta_line%s'%tube_r2] = np.load(glob_files(cfg, tube_r,'*beta*npy')[0]).T
    else:
        kwargs['beta_line%s'%tube_r2] = None
    return kwargs

def get_single_bpert(cfg, tube_r, search='*bpert*npy', beta=True,
                     triple_triple=False):
    """
    Read in the args for single_plot from files
    """
    kwargs = {}
    bppar, bpperp, bpphi = map(np.load,glob_files(cfg, tube_r, search))
    if triple_triple:
        tube_r2 = '_'+tube_r
    else:
        tube_r2 =''

    [kwargs['par_line%s'%tube_r2],
     kwargs['perp_line%s'%tube_r2],
     kwargs['phi_line%s'%tube_r2]] = bppar/1e11, bpperp/1e11, bpphi/1e11


    kwargs['x%s'%tube_r2],kwargs['y%s'%tube_r2] = get_xy(cfg, tube_r)
    if beta:
        kwargs['beta_line%s'%tube_r2] = np.load(glob_files(cfg, tube_r,'*beta*npy')[0]).T
    else:
        kwargs['beta_line%s'%tube_r2] = None
    return kwargs

def get_single_percentage_flux(cfg, tube_r, search='LineFlux*Fp*npy', beta=True,
                               triple_triple=False):
    """
    Read in the args for single_plot from files
    """
    kwargs = {}
    Fpar, Fperp, Fphi = map(np.load,glob_files(cfg, tube_r, search))
    Ftot = np.sqrt(Fpar**2 + Fperp**2 + Fphi**2)
    Fpar[Ftot<1e-5], Fperp[Ftot<1e-5], Fphi[Ftot<1e-5] = 0., 0., 0.

    if triple_triple:
        tube_r2 = '_'+tube_r
    else:
        tube_r2 =''

    [kwargs['par_line%s'%tube_r2],
     kwargs['perp_line%s'%tube_r2],
     kwargs['phi_line%s'%tube_r2]] = Fpar/Ftot, Fperp/Ftot, Fphi/Ftot


    kwargs['x%s'%tube_r2],kwargs['y%s'%tube_r2] = get_xy(cfg, tube_r)
    if beta:
        kwargs['beta_line%s'%tube_r2] = np.load(glob_files(cfg, tube_r,'*beta*npy')[0]).T
    else:
        kwargs['beta_line%s'%tube_r2] = None

    return kwargs

def get_triple(cfg, single='velocity',**kwergs):
    """
    Read in the args for a triple triple plot
    """
    if single == 'velocity':
        get_single = get_single_velocity
    elif single == 'percentage_flux':
        get_single = get_single_percentage_flux
    elif single == 'bpert':
        get_single = get_single_bpert
    kwargs = {}
    kwargs.update(get_single(cfg, 'r10', triple_triple=True, **kwergs))
    kwargs.update(get_single(cfg, 'r30', triple_triple=True, **kwergs))
    kwargs.update(get_single(cfg, 'r60', triple_triple=True, **kwergs))
    return kwargs


#==============================================================================
# Latex helpers
#==============================================================================
# Some LaTeX code snippiets
fig_str = r"""
\begin{{figure*}}
    {myfig}
    \caption{{{caption}}}
    \label{{{label}}}
\end{{figure*}}
"""

sub_fig = r"""\begin{{subfigure}}[{placement}]{{{width}}}
        {myfig}
        \caption{{{caption}}}
        \label{{{label}}}
    \end{{subfigure}}"""

def save_fig(cfg, fig_count=0, fig=None, fname=None):
    """
    Save a figure into the figures directory
    """

    if fname is None:
        fname = 'paperfig{}'.format(fig_count)
    if fig is None:
        fig = plt.gcf()

    fname = os.path.join(cfg.fig_dir, fname)
    fig.savefig(fname)

    return fname

def get_pgf_include(fname, width=r"\columnwidth"):
    return r"\pgfimage[width={1}]{{{0}}}".format(fname, width)

def get_fig_env(figure_str, fig_count=0, **kwargs):
    """
    Print a whole figure environment
    """
    global fig_str

    default_kwargs = {'placement':'H', 'caption':'caption',
                      'label':'fig_{}'.format(fig_count)}

    default_kwargs.update(kwargs)

    return fig_str.format(myfig=figure_str, **default_kwargs)

def get_subfig_env(subfig_str, caption, label, subplacement='b',
                   subwidth=r'\textwidth'):

    return sub_fig.format(placement=subplacement, width=subwidth,
                          caption=caption, label=label, myfig=subfig_str)

def make_fig_env(cfg, **kwargs):
    fname = kwargs.pop('fname', None)
    fig = kwargs.pop('fig', None)
    width = kwargs.pop('width', r"\columnwidth")

    fname = save_fig(cfg, fig=fig, fname=fname)
    fig_str = get_pgf_include(fname, width=width)

    print(get_fig_env(fig_str, **kwargs))

    return fname

def make_sub_figs(cfg, figs, captions, fig_count=0, labels=None, 
                  subplacement='b', subwidth=r'\textwidth', subbefore='\n',
                  subafter='\n', **kwargs):

    if len(figs) != len(captions) != len(labels):
        raise ValueError("Please specify the same number of labels as figures")

    sub_figs_str = r""
    fnames = []

    if labels is None:
        labels = ['fig_{fig_count}']*len(figs)

    for fig, caption, label in itertools.izip(figs, captions, labels):
        fname = save_fig(cfg, fig=fig)
        fnames.append(fname)
        pgf = get_pgf_include(fname, width=subwidth)
        label = label.format(fig_count=fig_count)
        sub_figs_str += subbefore + get_subfig_env(pgf, caption, label,
                                                   subplacement=subplacement,
                                                  subwidth=subwidth) + subafter

    print(get_fig_env(sub_figs_str, **kwargs))
    return fnames
