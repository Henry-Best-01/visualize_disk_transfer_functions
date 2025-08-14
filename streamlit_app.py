import numpy as np
from amoeba.Classes.accretion_disk import AccretionDisk
from amoeba.Classes.blr_streamline import Streamline
from amoeba.Classes.blr import BroadLineRegion
from amoeba.Util.util import (
    create_maps,
    calculate_gravitational_radius,
    accretion_disk_temperature
)
import matplotlib.pyplot as plt
import streamlit as st
from astropy.io import fits
from pandas import DataFrame as df


path_to_raytraces = "data/"

st.title("The Accretion Disk")
st.write("This is a toy GUI which uses Amoeba as described in Best et al. 2025, designed to explore how the flux distribution of the accretion disk and broad line region (BLR) is related to various parameters.")
st.write("The accretion disk model stems from the Shakura-Sunyaev thin-disk model which includes general relatvisitc (GR) corrections.")
st.write("Additions to this model are the lamp-post heating term as outlined in Cackett et al. 2007 and the Disk-wind term from Yong et al. 2019.")
st.write("GR ray tracing is computed with Sim5 as described in Bursa 2017.")
col_stack = st.columns(5)
col_stack[0].link_button("Best et al., 2025", "https://ui.adsabs.harvard.edu/abs/2025MNRAS.539.1269B/abstract")
col_stack[1].link_button("Shakura + Sunyaev, 1973", "https://ui.adsabs.harvard.edu/abs/1973A%26A....24..337S/abstract")
col_stack[2].link_button("Cackett et al., 2007", "https://ui.adsabs.harvard.edu/abs/2007MNRAS.380..669C/abstract")
col_stack[3].link_button("Yong et al., 2019", "https://ui.adsabs.harvard.edu/abs/2017PASA...34...42Y/abstract")
col_stack[4].link_button("Bursa, 2017", "https://ui.adsabs.harvard.edu/abs/2017bhns.work....7B/abstract")


left_col, right_col = st.columns(2)

mexp = left_col.slider("mass exponent", min_value=6.0, max_value=10.0, step=0.1, value=8.0)
redshift = right_col.slider("redshift", min_value=0.0, max_value=9.0, step=0.1, value=1.0)
inclination = left_col.slider("inclination angle [deg.]", min_value=2, max_value=89, step=3, value=20)
edd_ratio = left_col.slider("Eddington ratio", min_value=0.01, max_value=0.3, step=0.01, value=0.10)
wind_beta = right_col.slider("wind strength", min_value=0.0, max_value=0.8, step=0.01, value=0.0)

axis_range = right_col.slider(r"axis range [$r_{\rm{g}}$]", min_value=10, max_value=1000, step=10, value=100)
apply_gr = st.toggle("apply GR")
wavelength = st.slider("observer frame wavelength range [nm]", min_value=150, max_value=2000, step=5, value=(400, 600))

# grab the GR-ray trace
fname = path_to_raytraces+"RayTrace"+str(int(inclination)).zfill(2)+".fits"
with fits.open(fname) as f:
    r_map = f[0].data
    phi_map = f[1].data
    g_map = f[2].data
    header = f[0].header
    
# work on conversion to mags
min_wavelength = np.min(wavelength) * 10**-9
max_wavelength = np.max(wavelength) * 10**-9
min_frequency = 3e8 / max_wavelength
max_frequency = 3e8 / min_wavelength
delta_freq = abs(max_frequency - min_frequency)
delta_lam = abs(max_wavelength - min_wavelength)


# do some amoeba construction
acc_disk_dict = create_maps(
    mexp,
    redshift_source=redshift+0.0001,
    number_grav_radii=header['numgrs'],
    inclination_angle=inclination,
    resolution=np.size(r_map, 0),
    eddington_ratio=edd_ratio,
    visc_temp_prof="NT",
    temp_beta=wind_beta
)

# adjust maps to include GR
if apply_gr:
    acc_disk_dict['radii_array'] = r_map
    acc_disk_dict['phi_array'] = phi_map
    acc_disk_dict['g_array'] = g_map
    grav_rad = calculate_gravitational_radius(10**mexp)
    t_map = accretion_disk_temperature(r_map * grav_rad, 6.0 * grav_rad, 10**mexp, edd_ratio, beta=wind_beta, visc_temp_prof="NT")
    acc_disk_dict['temp_array'] = t_map


disk = AccretionDisk(**acc_disk_dict)

emission_array = disk.calculate_surface_intensity_map(np.mean(wavelength))
flux_array = emission_array.flux_array
total_flux = emission_array.total_flux
flux_exp = round(np.log10(total_flux), 0)
X, Y = emission_array.get_plotting_axes()


lum_dist = disk.lum_dist
approx_ab_mag = round(-2.5 * np.log10(total_flux / (4 * np.pi * lum_dist**2) * delta_lam / delta_freq * 1000) - 48.6, 1)

title_string = "emitted flux density: "+str(total_flux)[:4]+"e"+str(flux_exp)[:-2]+" W/m"+r", AB mag $\approx$"+str(approx_ab_mag)

fig, ax = plt.subplots(figsize=(8, 6))
contours = ax.contourf(X, Y, (flux_array), 50)
cbar = plt.colorbar(contours, ax=ax, label=r'flux density [W/m$^{2}$/m]')
ax.set_xlabel(r"x [$r_{\rm{g}}$]")
ax.set_ylabel(r"y [$r_{\rm{g}}$]")
ax.set_title(title_string)
ax.set_xlim(-axis_range, axis_range)
ax.set_ylim(-axis_range, axis_range)
ax.set_aspect(1)


st.write(fig)

fig2, ax2 = plt.subplots(figsize=(8, 3))
ax2.plot(X[0, :], flux_array[500, :])
ax2.set_xlim(-axis_range, axis_range)
ax2.set_xlabel(r"x [$r_{\rm{g}}$]")
ax2.set_ylabel(r"flux density across center [W/m$^{2}$/m]")

st.write(fig2)

@st.fragment()
def save_data(current_data, current_metadata, my_key):
    file_name = st.text_input("Type output file name below", value="output.csv", key=my_key)
    metadata_file_name = file_name[:file_name.find(".")]+"_metadata.txt"
    output_data = df(current_data).to_csv().encode("utf-8")
    st.download_button("Download "+my_key+" flux distribution", output_data, file_name, key=my_key+" button")
    st.download_button("Download "+my_key+" metadata", current_metadata, metadata_file_name, key=my_key+" metadata button")

current_metadata = "Refer to this file for further information \n"+f"axis_range: (-1000, 1000) \n"+f"mass_exponent: {mexp} \n"+f"redshift: {redshift} \n"+f"inclination: {inclination} \n"+f"eddington_ratio: {edd_ratio} \n"+f"with GR: {apply_gr} \n"+f"wavelength range: {wavelength}"               

save_data(flux_array, current_metadata, my_key="accretion disk")


st.title("The BLR")

left_col, right_col = st.columns(2)

blr_radii = left_col.slider(r"wind radius [$r_{\rm{g}}$]", min_value=10, max_value=2000, step=10, value=(100, 500))
blr_angles = right_col.slider("launch angles [deg.]", min_value=1, max_value=80, step=1, value=(20, 45))
blr_characteristic_distance_1 = left_col.slider(r"inner characteristic distance [$r_{\rm{g}}$]", min_value=10, max_value=2000, step=10, value=200)
blr_characteristic_distance_2 = right_col.slider(r"outer characteristic distance [$r_{\rm{g}}$]", min_value=10, max_value=2000, step=10, value=200)
asympt_vel_1 = left_col.slider("inner asymptotic velocity [c]", min_value=0.001, max_value=0.4, step=0.001, value=0.05)
asympt_vel_2 = right_col.slider("outer asymptotic velocity [c]", min_value=0.001, max_value=0.4, step=0.001, value=0.05)

blr_opt_radius = left_col.slider(r"optimal blr emission radius [$r_{\rm{g}}$]", min_value=10, max_value=4000, step=10, value=300)
blr_opt_width = right_col.slider(r"optimal blr emission width [$r_{\rm{g}}$]", min_value=10, max_value=1000, step=10, value=50)
emitted_wavelength = st.slider("rest frame emission line [nm]", min_value=50, max_value=2000, step=1, value=210)


# Do the Amoeba production of a BLR object
stream_1 = Streamline(blr_radii[0], blr_angles[0], 1000, blr_characteristic_distance_1+0.001, asympt_vel_1, height_step=20)
stream_2 = Streamline(blr_radii[1], blr_angles[1], 1000, blr_characteristic_distance_2+0.001, asympt_vel_2, height_step=20)

my_blr = BroadLineRegion(mexp, 1000, emitted_wavelength, redshift, radial_step=20, height_step=20)
my_blr.add_streamline_bounded_region(stream_1, stream_2)

# define some emission efficiency array based on the optimal emission radius
R, Z = my_blr.get_density_axis()
r_spherical = np.sqrt(R**2 + Z**2)
efficiency = np.exp(-(r_spherical-blr_opt_radius)**2/(2 * blr_opt_width**2))

blr_projection = my_blr.project_blr_intensity_over_velocity_range(inclination, observed_wavelength_range_in_nm=wavelength, emission_efficiency_array=efficiency)

X, Y = blr_projection.get_plotting_axes()
blr_flux = blr_projection.flux_array

fig3, ax3 = plt.subplots(figsize=(8, 6))
contours3 = ax3.contourf(X, Y, blr_flux, 50)
cbar3 = plt.colorbar(contours3, ax=ax3, label="relative blr emission")
ax3.set_xlabel(r"x [$r_{\rm{g}}$]")
ax3.set_ylabel(r"y [$r_{\rm{g}}$]")
ax3.set_title(f"BLR emitting at rest frame {emitted_wavelength} nm \n observed in range {wavelength[0]}-{wavelength[1]} nm at z={redshift}")
ax3.set_aspect(1)

st.write(fig3)

st.write("Note that this figure may appear blank if the emission line does not project into your range of wavelengths.")

current_metadata = "Refer to this file for further information \n"+f"x_axis_range: ({np.min(X)}, {np.max(X)}) \n"+f"y_axis_range: ({np.min(Y)}, {np.max(Y)}) \n"+f"mass_exponent: {mexp} \n"+f"redshift: {redshift} \n"+f"inclination: {inclination} \n"+f"wavelength range: {wavelength} \n"+f"wind launching radii: {blr_radii} \n"+f"wind launching angles: {blr_angles} \n"+f"inner characteristic distance: {blr_characteristic_distance_1} \n"+f"outer characteristic distance: {blr_characteristic_distance_2} \n"+f"inner asymptotic velocity: {asympt_vel_1} \n"+f"outer asymptotic velocity: {asympt_vel_2} \n"+f"optimal BLR emission radius: {blr_opt_radius} \n"+f"BLR emission width: {blr_opt_width} \n"+f"rest frame wavelength: {emitted_wavelength}"

save_data(blr_flux, current_metadata, my_key="broad line region")

st.write("Here are the relevant R-Z projections of the density and poloidal velocities")

fig4, ax4 = plt.subplots(1, 2, figsize=(10, 4))

densities = my_blr.density_grid
velocities = np.sqrt(my_blr.z_velocity_grid**2 + my_blr.r_velocity_grid**2)

density_contours = ax4[0].contourf(R, Z, np.log10(densities*efficiency/np.nanmax(densities*efficiency)), 50)
velocity_contours = ax4[1].contourf(R, Z, velocities, 50)

density_cbar = plt.colorbar(density_contours, ax=ax4[0], label="log relative particle density * efficiency [arb.]")
velocity_cbar = plt.colorbar(velocity_contours, ax=ax4[1], label="outflowing velocities [c]")

plt.subplots_adjust(wspace=0.4)

for axis in ax4:
    axis.set_xlabel(r"R [$r_{\rm{g}}$]")
    axis.set_ylabel(r"Z [$r_{\rm{g}}$]")
    axis.set_aspect(1)

st.write(fig4)

st.title("Accretion disk params")
st.write("The sliders provided represent physical parameters of the active galactic nucleus.")
st.write(r"The mass exponent is related to the physical mass of the central supermassive black hole as mass_exp = log$_{10}(M_{\bullet}/M_{\odot})$. This also defines the size scale--the gravitational radius $r_{\rm{g}} = GM_{\bullet}c^{-2}$.")  
st.write("The inclination angle represents the observer's position with respect to the axis of symmetry. Zero degrees is face-on, while ninety degrees is edge-on.")
st.write("The Eddington ratio is related to how much matter is drawn into the supermassive black hole in terms of the Eddington limit. The Eddington limit is the amount of matter required to balance the gravitational force with radiation pressure in a system with zero angular momentum. Thin disks are not effective above ~0.3 Eddington ratio, where other accretion modes are required to affectively accrete material.")
st.write("The redshift of the system represents the physical distance. This impacts everything from the wavelengths observed to its relative intensity. Within this model, the Flat Lambda-CDM Universe is assumed with H0 = 70 km/s/Mpc, Omega_M = 0.3, Omega_0 = 0.7.")

st.write("The wind strength parameter represents how much material from the accretion disk is 'blown off' the accretion disk as it moves inwards to the black hole. Zero represents no material is removed, while 0.8 represents 80 percent of the material is removed every 6 $r_{\rm{g}}$.")
st.write(r"Axis range simply controls the zoom of the plots in the accretion disk section in units of $r_{\rm{g}}$.")
st.write("The 'apply GR' toggle allows you to switch between natively calculated accretion disks and those computed with GR corrections. Most prominent effects are the light bending over the black hole and the relativistic Doppler boosting of the approaching side (left) for high inclinations. The GR effects introduce redshift and blueshift to the accretion disk depending on line-of-sight velocities.")
st.write("The observer frame wavelength range represents the wavelengths an observer sees. Amoeba is designed for optical, ultra-violet, and infrared emission. Redshifting of the filter range is applied and an estimate of the AB-magnitude are related to this range.")

st.title("Broad line region params")
st.write("Wind radius represents the base of where the BLR wind is launched from. Amoeba assumes a biconical outflow for some region of the accretion disk. Note that not all regions contained within this wind radius need to be at the proper ionization state to emit the simulated emission line (see optimal blr emission radius).")
st.write("inner/outer characteristic distance represent the poloidal distance at which the velocity accelerates to the asymptotic velocity. After one characteristic distance, the wind will be outflowing poloidally at 0.5 * asymptotic velocity. ")
st.write("innner/outer asymptotic velocity represents the maximum outflowing velocity as the poloidal distance extends towards infinity.")
st.write("optimal BLR emission radius is a proxy for the ideal emission region of the BLR. Closer than this region, it is assumed that the simulated emission species will be completely ionized. Further out, it is not likely for these species to emit their emission line efficiently.")
st.write("optimal BLR emission radius represents how wide this optimal emission region is. This current model is a radial Gaussian, and in practice the emission efficiency should be calculated with a photoionization code.")
st.write("launch angles represent the inner and outer angles of the BLR biconical geometry, in degrees. Zero degrees represents normal to the accretion disk, while ninety degrees represents traveling along the surface of the accretion disk.")
st.write("rest frame emission line represents the wavelength of the BLR if it wasn't moving and was here next to the observer. This may be tied to realistic emission lines, and the emission will then be transformed into the observer frame's reference. This effectively allows a calculation of the line shape as you change the filter, but is not included in this app because it will have to compute all wavelengths every time something is changed.")




