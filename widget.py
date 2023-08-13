#!/usr/bin/env python
# coding: utf-8

import optool
import matplotlib.pyplot as plt
import numpy as np
from functools import lru_cache
from matplotlib.widgets import Slider
import matplotlib.pyplot as plt
import numpy as np

@lru_cache(maxsize=10)
def optool_call(s):
    """
    Call optool with caching of results for performance.

    Parameters
    ----------
    s : str
        Command line string for optool.

    Returns
    -------
    optool.particle : object
        Optool particle object constructed from parameters in s.
    """
    p = optool.particle(s,None,True)
    return p

def compute_model(amin, amax, exp, por, material, wavelength_micron, sizebins):
    """
    Compute the phase function for a given set of parameters and materials.

    Materials with abundances lower than 1e-2 are ignored.

    Example material dict:
    material = {
        "pyr-mg70": 0.87,
        "c": 0.13,
    }

    Parameters
    ----------
    amin : float
        Minimum grain size.
    amax : float
        Maximum grain size.
    exp : float
        Exponent of the power-law grain size distribution.
    por : float
        Porosity.
    material : dict
        Dictionary with material names and their corresponding abundances.
    wavelength_micron : float
        Wavelength in micrometers.
    sizebins : int
        Number of size bins for the grain size distribution.

    Returns
    -------
    angles_model : ndarray
        Array of angles for which the model is computed.
    pol_flux_model : ndarray
        Computed polarized flux for each angle.
    """
    materials_str = ""
    for m, x in material.items():
        # Ignore materials with abundance < 0.01
        if float(x) > 1e-2:
            materials_str += f"{m} {x} "
    print(materials_str)
    cmd = f"optool {materials_str} -p {por} -a {amin} {amax} {exp} {sizebins} -l {wavelength_micron} -s"
    p = optool_call(cmd)
    f12 = np.abs(p.f12[0][0])
    angles_model = np.arange(0, 180)
    pol_flux_model = f12/f12[len(f12)//2]
    return angles_model, pol_flux_model

def compute_model_nk(amin, amax, exp, por, n, k, wavelength_micron, sizebins):
    """
    Compute the phase function for a given set of parameters and materials.

    Materials with abundances lower than 1e-2 are ignored.

    Parameters
    ----------
    amin : float
        Minimum grain size.
    amax : float
        Maximum grain size.
    exp : float
        Exponent of the power-law grain size distribution.
    por : float
        Porosity.
    material : dict
        Dictionary with material names and their corresponding abundances.
    wavelength_micron : float
        Wavelength in micrometers.
    sizebins : int
        Number of size bins for the grain size distribution.

    Returns
    -------
    angles_model : ndarray
        Array of angles for which the model is computed.
    pol_flux_model : ndarray
        Computed polarized flux for each angle.
    """
    cmd = f"optool -p {por} -a {amin} {amax} {exp} {sizebins} -l {wavelength_micron} -c {n}:{k}:1.0 -s"
    p = optool_call(cmd)
    f12 = np.abs(p.f12[0][0])
    angles_model = np.arange(0, 180)
    pol_flux_model = f12/f12[len(f12)//2]
    return angles_model, pol_flux_model

def polarized_flux_widget(angles_obs, 
                          pol_flux_obs, 
                          materials, 
                          wavelength_micron,     
                          init_amin=0.1,
                          init_amax=1000,
                          init_exp=3,
                          init_por=0.25):
    """
    Interactively display a plot of observed and modeled polarized flux, allowing user
    to adjust parameters via sliders.

    Parameters
    ----------
    angles_obs : ndarray
        Observed angles.
    pol_flux_obs : ndarray
        Observed polarized flux for each angle.
    materials : dict
        Dictionary with material names and their initial relative abundances.
    wavelength_micron : float
        Wavelength in micrometers for the plot.
    init_amin : float, optional
        Initial minimum grain size (default is 0.1).
    init_amax : float, optional
        Initial maximum grain size (default is 1000).
    init_exp : float, optional
        Initial exponent of the grain size distribution (default is 3).
    init_por : float, optional
        Initial porosity (default is 0.25).
    """

    # Create a figure with axes for each parameter slider and the plot.
    axes_names = [
        "plot", ".", ".", ".", "slamin", "slamax", "slexp", "slpor"
    ] + [f"sl{m}" for m in materials]

    height_ratios = [1] + [0.03] * (len(axes_names) - 1)

    fig, axd = plt.subplot_mosaic([[a] for a in axes_names],
                                  height_ratios=height_ratios,
                                  figsize=(8, 6+0.25*(len(axes_names)-2)))

    angles_model, pol_flux_model = compute_model(
        init_amin, init_amax, init_exp, init_por, materials, wavelength_micron=wavelength_micron, sizebins=30)
    
    line, = axd["plot"].plot(angles_model, pol_flux_model, label="model")
    axd["plot"].plot(angles_obs, pol_flux_obs, label="observed")
    axd["plot"].legend()
    axd["plot"].set_xlabel("Angle (deg)")
    axd["plot"].set_ylabel("Normalized flux")
    axd["plot"].set_ylim(0, 1.2*np.max(pol_flux_obs))

    # Make a horizontal sliders to control model parameters.
    sliders = {}

    amin_slider = Slider(
        ax=axd["slamin"],
        label='amin [micron]',
        valmin=0.1,
        valmax=1,
        valinit=init_amin,
    )
    sliders["amin"] = amin_slider

    amax_slider = Slider(
        ax=axd["slamax"],
        label='amax [micron]',
        valmin=50,
        valmax=2000,
        valinit=init_amax,
    )
    sliders["amax"] = amax_slider

    exp_slider = Slider(
        ax=axd["slexp"],
        label='exp',
        valmin=0.1,
        valmax=5,
        valinit=init_exp,
    )
    sliders["exp"] = exp_slider

    por_slider = Slider(
        ax=axd["slpor"],
        label='por',
        valmin=0,
        valmax=1,
        valinit=init_por,
    )
    sliders["por"] = por_slider

    for m, c in materials.items():
        sl = Slider(
            ax=axd[f"sl{m}"],
            label=m,
            valmin=0,
            valmax=1,
            valinit=c,
        )
        sliders[m] = sl

    # The function to be called anytime a slider's value changes
    def update(val):
        amin = amin_slider.val
        amax = amax_slider.val
        exp = exp_slider.val
        por = por_slider.val
        mats = {m: sliders[f"{m}"].val for m in materials.keys()}
        angles_model, pol_flux_model = compute_model(
            amin, amax, exp, por, mats, wavelength_micron=wavelength_micron, sizebins=30)
        line.set_ydata(pol_flux_model)
        fig.canvas.draw_idle()

    # Update the plot when mouse button is released instead of on update.
    # Otherwise the widget quickly becomes unresponsive due to exessive optool calls.
    fig.canvas.mpl_connect("button_release_event", update)

    plt.show()

# Specify the materials used in the model 
materials = {
    "pyr-mg70": 0.87,
    "c": 0.13,
    "ol-mg40": 0.0,
    "c-z": 0.0,
    "c-p": 0.0,
    "c-gra": 0.0,
    "c-nano": 0.0,
    "c-org": 0.0,
    "sio2": 0.0,
    "cor-c": 0.0,
    "fe-c": 0.0,
    "fes": 0.0,
    "sic": 0.0,
    "h2o-w": 0.0,
    "h2o-a": 0.0,
    "co2-w": 0.0,
    "nh3-m": 0.0,
    "co-a": 0.0,
    "co2-a": 0.0,
    "co2-c": 0.0,
    "ch4-a": 0.0,
    "ch4-c": 0.0,
    "ch3oh-a": 0.0,
    "ch3oh-c": 0.0
}


def fit_nk_widget(angles_obs, 
                  flux_obs,
                  wavelength_micron,
                  init_n = 1.0,  init_k = 0.1,
                  init_amin=0.1, init_amax=1000,
                  init_exp=3.5,
                  init_por=0.25):
    """
    Interactively display a plot of observed and modeled polarized flux, allowing user
    to adjust parameters via sliders.

    Parameters
    ----------
    angles_obs : ndarray
        Observed angles.
    flux_obs : ndarray
        Observed polarized flux for each angle.
    wavelength_micron : float
        Wavelength in micrometers for the plot.
    init_n : float
        Initial value for real part of refractive index, default is 1.0
    init_k : float
        Initial value for imag part of refractive index, default is 0.1
    init_amin : float, optional
        Initial minimum grain size (default is 0.1).
    init_amax : float, optional
        Initial maximum grain size (default is 1000).
    init_exp : float, optional
        Initial exponent of the grain size distribution (default is 3.5).
    init_por : float, optional
        Initial porosity (default is 0.25).
    """

    # Create a figure with axes for each parameter slider and the plot.
    axes_names = [
        "plot", ".", ".", ".", "slamin", "slamax", "slexp",
        "slpor", "sln", "slk" ]

    height_ratios = [1] + [0.03] * (len(axes_names) - 1)

    fig, axd = plt.subplot_mosaic([[a] for a in axes_names],
                                  height_ratios=height_ratios,
                                  figsize=(8, 6+0.25*(len(axes_names)-2)))

    angles_model, pol_flux_model = compute_model_nk(
        init_amin, init_amax, init_exp, init_por, init_n, init_k,
        wavelength_micron=wavelength_micron, sizebins=30)
    
    line, = axd["plot"].plot(angles_model, pol_flux_model, label="model")
    axd["plot"].plot(angles_obs, flux_obs, label="observed")
    axd["plot"].legend()
    axd["plot"].set_xlabel("Angle (deg)")
    axd["plot"].set_ylabel("Normalized flux")
    axd["plot"].set_ylim(0, 1.2*np.max(flux_obs))

    # Make horizontal sliders to control model parameters.
    sliders = {}

    amin_slider = Slider(
        ax=axd["slamin"],
        label='amin [micron]',
        valmin=0.1,
        valmax=1,
        valinit=init_amin,
    )
    sliders["amin"] = amin_slider

    amax_slider = Slider(
        ax=axd["slamax"],
        label='amax [micron]',
        valmin=50,
        valmax=2000,
        valinit=init_amax,
    )
    sliders["amax"] = amax_slider

    exp_slider = Slider(
        ax=axd["slexp"],
        label='exp',
        valmin=0.1,
        valmax=5,
        valinit=init_exp,
    )
    sliders["exp"] = exp_slider

    por_slider = Slider(
        ax=axd["slpor"],
        label='por',
        valmin=0,
        valmax=1,
        valinit=init_por,
    )
    sliders["por"] = por_slider

    n_slider = Slider(
        ax=axd["sln"],
        label='n',
        valmin=0,
        valmax=10,
        valinit=init_n,
    )
    sliders["n"] = n_slider

    k_slider = Slider(
        ax=axd["slk"],
        label='k',
        valmin=0.1,
        valmax=10,
        valinit=init_k,
    )
    sliders["k"] = k_slider

    # The function to be called anytime a slider's value changes
    def update(val):
        amin = amin_slider.val
        amax = amax_slider.val
        exp = exp_slider.val
        por = por_slider.val
        n   = n_slider.val
        k   = k_slider.val
        angles_model, pol_flux_model = compute_model_nk(
            amin, amax, exp, por, n, k, wavelength_micron=wavelength_micron, sizebins=30)
        line.set_ydata(pol_flux_model)
        fig.canvas.draw_idle()

    # Update the plot when mouse button is released instead of on update.
    # Otherwise the widget quickly becomes unresponsive due to exessive optool calls.
    fig.canvas.mpl_connect("button_release_event", update)

    plt.show()

# Specify the materials used in the model 
materials = {
    "pyr-mg70": 0.87,
    "c": 0.13,
    "ol-mg40": 0.0,
    "c-z": 0.0,
    "c-p": 0.0,
    "c-gra": 0.0,
    "c-nano": 0.0,
    "c-org": 0.0,
    "sio2": 0.0,
    "cor-c": 0.0,
    "fe-c": 0.0,
    "fes": 0.0,
    "sic": 0.0,
    "h2o-w": 0.0,
    "h2o-a": 0.0,
    "co2-w": 0.0,
    "nh3-m": 0.0,
    "co-a": 0.0,
    "co2-a": 0.0,
    "co2-c": 0.0,
    "ch4-a": 0.0,
    "ch4-c": 0.0,
    "ch3oh-a": 0.0,
    "ch3oh-c": 0.0
}




#fit_nk_widget(np.linspace(0.5,179.5,180), np.linspace(0.,1.,180),
#                      wavelength_micron=1.25)
#polarized_flux_widget(angles_obs, pol_flux_obs,
#                      materials, wavelength_micron=1.25)
polarized_flux_widget(np.linspace(0.5,179.5,180), np.linspace(0.,1.,180),
                      materials, wavelength_micron=1.25)
