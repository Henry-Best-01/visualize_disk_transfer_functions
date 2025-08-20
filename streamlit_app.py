import numpy as np
from amoeba.Classes.accretion_disk import AccretionDisk
from amoeba.Util.util import (
    create_maps,
    calculate_gravitational_radius,
    accretion_disk_temperature
)
import matplotlib.pyplot as plt
import streamlit as st
from astropy.io import fits
import astropy.constants as const
import astropy.units as u


path_to_raytraces = "data/"

st.title("The reverberating Accretion Disk")
st.write("This is a toy GUI which uses Amoeba as described in Best et al. 2025 to compute the accretion disk's response function using the lamp-post model.")
st.write("The accretion disk model stems from the Shakura-Sunyaev thin-disk model which includes general relatvisitc (GR) corrections.")
st.write("Additions to this model are the lamp-post heating term as outlined in Cackett et al. 2007 and the Disk-wind term from Yong et al. 2019.")
st.write("GR ray tracing is computed with Sim5 as described in Bursa 2017.")
st.write("References:")
col_stack = st.columns(5)
col_stack[0].link_button("Best et al., 2025", "https://ui.adsabs.harvard.edu/abs/2025MNRAS.539.1269B/abstract")
col_stack[1].link_button("Shakura + Sunyaev, 1973", "https://ui.adsabs.harvard.edu/abs/1973A%26A....24..337S/abstract")
col_stack[2].link_button("Cackett et al., 2007", "https://ui.adsabs.harvard.edu/abs/2007MNRAS.380..669C/abstract")
col_stack[3].link_button("Yong et al., 2019", "https://ui.adsabs.harvard.edu/abs/2017PASA...34...42Y/abstract")
col_stack[4].link_button("Bursa, 2017", "https://ui.adsabs.harvard.edu/abs/2017bhns.work....7B/abstract")


## The definition of the disk can't be in a fragment, since changing it 
# means we need to recalculate everything later. 
def define_disk():
    left_col, right_col = st.columns(2)
    mexp = left_col.slider("mass exponent", min_value=7.0, max_value=10.0, step=0.1, value=8.0)
    inclination = left_col.slider("inclination angle [deg.]", min_value=2, max_value=89, step=3, value=20)
    edd_ratio = left_col.slider("Eddington ratio", min_value=0.01, max_value=0.3, step=0.01, value=0.10)
    wind_beta = right_col.slider("wind strength", min_value=0.0, max_value=0.8, step=0.01, value=0.0)
    redshift = right_col.slider("redshift", min_value=0.0, max_value=9.0, step=0.1, value=1.0)
    apply_gr = left_col.toggle("apply GR")
    use_nt = right_col.toggle("use Novikov-Thorne")
    if not apply_gr:
        disk_flaring = st.slider("disk flaring", min_value=0, max_value=45, step=1, value=0)
    else:
        st.write("disk flaring is disabled due to its incompatibility with the GR ray traced maps")
        disk_flaring = 0
    param_bundle = [mexp, inclination, edd_ratio, wind_beta, redshift]

    # grab the GR-ray trace
    fname = path_to_raytraces+"RayTrace"+str(int(inclination)).zfill(2)+".fits"
    with fits.open(fname) as f:
        r_map = f[0].data
        phi_map = f[1].data
        g_map = f[2].data
        header = f[0].header

    if use_nt:
        temp_prof = "NT"
    else:
        temp_prof = "SS"

    # do some amoeba construction
    acc_disk_dict = create_maps(
        mexp,
        redshift_source=redshift+0.0001,
        number_grav_radii=header['numgrs'],
        inclination_angle=inclination,
        resolution=np.size(r_map, 0),
        eddington_ratio=edd_ratio,
        visc_temp_prof=temp_prof,
        temp_beta=wind_beta
    )

    # adjust maps to include GR
    if apply_gr:
        acc_disk_dict['radii_array'] = r_map
        acc_disk_dict['phi_array'] = phi_map
        acc_disk_dict['g_array'] = g_map
        grav_rad = calculate_gravitational_radius(10**mexp)
        t_map = accretion_disk_temperature(r_map * grav_rad, 6.0 * grav_rad, 10**mexp, edd_ratio, beta=wind_beta, visc_temp_prof=temp_prof)
        acc_disk_dict['temp_array'] = t_map
    else:
        acc_disk_dict['height_array'] = r_map * np.tan(disk_flaring * np.pi / 180)
    disk = AccretionDisk(**acc_disk_dict)

    return disk, param_bundle

disk, param_bundle = define_disk()
mexp, inclination, edd_ratio, wind_beta, redshift = param_bundle


# This can be in a fragment, but it needs up also update the figures
@st.fragment
def define_transfer_function(disk):
    left_col, right_col = st.columns(2)
    flare_offset = left_col.slider(r"source offset [$r_{\rm{g}}$]", min_value=0, max_value=1000, step=5, value=0)
    angle_offset = right_col.slider(r"phi offset [deg.]", min_value=0, max_value=360, step=1, value=0)
    wavelength = st.slider("observer frame wavelength range [nm]", min_value=150, max_value=2000, step=5, value=400)
    response_function = disk.construct_accretion_disk_transfer_function(
        wavelength, 
        axis_offset_in_gravitational_radii=flare_offset, 
        angle_offset_in_degrees=angle_offset
    )
    mean_response = np.sum(response_function * np.linspace(
        0,
        len(response_function) - 1,
        len(response_function)
    )) / np.sum(response_function)
    centroid_mask = response_function > 0.5 * np.max(response_function)
    centroid_mean = np.sum(response_function * np.linspace(
        0,
        len(response_function) - 1,
        len(response_function)
    ) * centroid_mask) / np.sum(response_function * centroid_mask)
    response_map, time_lags = disk.construct_accretion_disk_transfer_function(
        wavelength, 
        axis_offset_in_gravitational_radii=flare_offset, 
        angle_offset_in_degrees=angle_offset,
        return_response_array_and_lags=True
    )

    flux_array = disk.calculate_surface_intensity_map(wavelength)

    lags_in_days = (time_lags * disk.rg / (const.c.to(u.m/u.day))).value * (1 + redshift)
    tau_ax_in_days = np.linspace(
        0,
        np.max(lags_in_days),
        len(response_function)
    )
    spacing = 0.1
    if np.max(lags_in_days) > 1:
        spacing = 1
    if np.max(lags_in_days) > 10:
        spacing = 10
    if np.max(lags_in_days) > 100:
        spacing = 100

    axis_range = st.slider(r"axis range [$r_{\rm{g}}$]", min_value=10, max_value=1000, value=1000)

    lag_contours = np.linspace(0, 2000, 11)
    lag_contours_in_days = np.linspace(
        0, 
        np.max(lags_in_days)//spacing * spacing,
        int(np.max(lags_in_days)//spacing + 1)
    )

    X, Y = flux_array.get_plotting_axes()

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.contourf(X, Y, np.log10(response_map), 20)
    ax.contour(X, Y, time_lags, lag_contours, colors='white', linewidths=0.5)
    ax.contour(X, Y, lags_in_days, lag_contours_in_days, colors='white', alpha=0.7, linewidths=0.5,linestyles='dashed')

    ax.set_title(f"dashed contours placed every {spacing} day(s) \n solid contours placed every 200" + r" $r_{\rm{g}}$")
    ax.set_xlabel(r"x [$r_{\rm{g}}$]")
    ax.set_ylabel(r"y [$r_{\rm{g}}$]")
    ax.set_xlim(-axis_range, axis_range)
    ax.set_ylim(-axis_range, axis_range)
    ax.set_aspect(1)
    st.write(fig)

    fig2, ax2 = plt.subplots(figsize=(8, 3))

    ax2.plot(response_function)
    ax2.set_xlabel(r"$\tau$ [$r_{\rm{g}}$]")
    ax2.set_ylabel(r"$\Psi$ [arb.]")

    
    x_lim_slider = st.slider(
        r"max time lag [$r_{\rm{g}}$]",
        min_value=0,
        max_value=2000,
        step=1,
        value=200,
    )
    ax2_twin = ax2.twiny()
    ax2_twin.plot(tau_ax_in_days, response_function)
    value_in_days = x_lim_slider * disk.rg / (const.c.to(u.m/u.day)).value * (1 + redshift)
    
    ax2_twin.set_xlabel(r"$\tau$ [day]")

    ax2.plot([mean_response, mean_response], [0, response_function[int(mean_response)]], '--', color='grey', alpha=0.7, label="total mean")
    ax2.plot([centroid_mean, centroid_mean], [0, response_function[int(centroid_mean)]], '-.', color='black', alpha=0.7, label="centroid mean")
    
    ax2.set_xlim(-1, x_lim_slider)
    ax2_twin.set_xlim(
        -disk.rg / (const.c.to(u.m/u.day)).value * (1 + redshift),
        value_in_days
    )
    ax2.legend()
    st.write(fig2)

define_transfer_function(disk)





