
import os.path
import numpy as np
import pickle

from .functions import read_data, topo_phase_removal, calculate_covariance_matrix, capon_beam_forming_inversion, \
    plot_tomo_slices, plot_profiles

# ****************************************************************************
# Main Function for SAR Tomography Processing
# Nesrin Salepci
# special thanks to Marivi Tello Alonso, Victor Cazcarra Bes and Felix Cremer
# January 2017 update: January 2018 Jena
# revision by John Truckenbrodt, March/April 2018
# ****************************************************************************
#
# ----------------------------------------------------------------------------
# ****************************************************************************
# flags for the functions execute =1, not-execute = 0
func_read_data = 1
func_topo_phase_removal = 1
func_calculate_covariance_matrix = 1
func_capon_beam_forming_inversion = 1
func_plot_tomo_slices = 1
func_plot_profiles = 1

# *****************************************************************************
# --------------------------- Define the Input Parameters ---------------------
# *****************************************************************************


inpath = "I:\\SAR-EDU_Tomography_Module\\Input_Files\\"
outpath = "I:\\SAR-EDU_Tomography_Module\\Results\\"

if not outpath or not inpath:
    print("-------------------------------------------------------------------")
    print("!!!!!")
    print("Define the input and output folders and their paths: path/folder")
    print("-------------------------------------------------------------------")

master = "SLC_0_20151013_L_hv"

slaves = ["SLC_1_20151013_L_hv",
          "SLC_2_20151013_L_hv",
          "SLC_3_20151013_L_hv",
          "SLC_4_20151013_L_hv",
          "SLC_5_20151013_L_hv"]

phases = ["Pha_1_20151013_L_hv",
          "Pha_2_20151013_L_hv",
          "Pha_3_20151013_L_hv",
          "Pha_4_20151013_L_hv",
          "Pha_5_20151013_L_hv"]

kzs = ["Kz_1_20151013_L_hv",
       "Kz_2_20151013_L_hv",
       "Kz_3_20151013_L_hv",
       "Kz_4_20151013_L_hv",
       "Kz_5_20151013_L_hv"]

# define the boxcar smoothing dimension
multi_look = 10

# define the max height for the inversion
height = 70

# define the selected pixels for reflectivity profile plotting:
# define the method of plotting 0: selected pixel, 1: around a central pixel
central_pixel = 0

if central_pixel == 1:
    # select a central pixel (range: midr and azimuth: mida) and plot the reflectivity profiles of the pixel
    # around the central pixel in a range of inrange (defined by the user)
    midr = 304
    mida = 17
    inrange = 5
    list_pixels = []
    for i_rg in range(0, inrange):
        for j_az in range(0, inrange):
            pixel = [midr + i_rg, mida + j_az]
            list_pixels.append(pixel)

if central_pixel == 0:
    # define pixel coordinates for plotting reflectivity profiles, [range pixel, azimuth pixel]
    list_pixels = [[189, 253], [190, 254], [195, 260], [304, 17], [306, 20], [350, 30]]

# ******************************************************************************
# end of input parameters
# ******************************************************************************


if not os.path.exists(outpath):
    os.makedirs(outpath)

SLC_list = list()
SLC_list.append(os.path.join(inpath, master))
print(SLC_list)
slaves_path = [os.path.join(inpath, slave)for slave in slaves]
print(slaves_path)
SLC_list.extend(slaves_path)
print(SLC_list)
nTrack = len(SLC_list)
numSlaves = len(slaves)

phase_list = list()
phase_path = [os.path.join(inpath, phase)for phase in phases]
phase_list.extend(phase_path)
print(phase_list)

kz_list = list()
kz_path = [os.path.join(inpath, kz)for kz in kzs]
kz_list.extend(kz_path)
print(kz_list)

j_complex = complex(0, 1)
height_vector = np.matrix(np.arange(-height, height+1, 1))
# *****************************************************************************
# ----------------------------------------------------------------------------
# *****************************************************************************
# =============================================================================
#     ========= Call the Functions =========
# =============================================================================
if func_read_data == 1:
    print('------------------------------------------------------------')
    print('read the SAR data ')
    print('------------------------------------------------------------')
    print('SLC_stack')
    filout = 'SLC_stack'
    SLC_stack = read_data(SLC_list, filout)

    print('phase_stack')
    filout = 'phase_stack'
    phase_stack = read_data(phase_list, filout)

    print('kz_stack')
    filout = 'kz_stack'
    kz_stack = read_data(kz_list, filout)

if func_topo_phase_removal == 1:
    print('------------------------------------------------------------')
    print('removal of the flat-earth phase and the topographical phase ')
    print('------------------------------------------------------------')
    # restoring the input variables
    SLC_stack = pickle.load(open(os.path.join(outpath + 'SLC_stack'), "rb"))
    phase_stack = pickle.load(open(os.path.join(outpath + 'phase_stack'), "rb"))

    Normalized_Stack = topo_phase_removal(SLC_stack, phase_stack)

if func_calculate_covariance_matrix == 1:
    print('------------------------------------------------------------')
    print('estimation of coherence and generation of covariance matrix')
    print('------------------------------------------------------------')
    # restoring the input variables
    Normalized_Stack = pickle.load(open(os.path.join(outpath + 'normalized_stack'), "rb"))

    Covariance_Matrix = calculate_covariance_matrix(Normalized_Stack, multi_look)

if func_capon_beam_forming_inversion == 1:
    print('------------------------------------------------------------')
    print('CAPON beam forming inversion')
    print('------------------------------------------------------------')
    # restoring the input variables
    Covariance_Matrix = pickle.load(open(os.path.join(outpath + 'cov_matrix'), "rb"))
    kz_stack = pickle.load(open(os.path.join(outpath + 'kz_stack'), "rb"))

    Capon_bf = capon_beam_forming_inversion(Covariance_Matrix, kz_stack, height_vector)

if func_plot_tomo_slices == 1:
    print('------------------------------------------------------------')
    print('Plot Tomoslices')
    print('------------------------------------------------------------')
    # restoring the input variable
    Capon_bf = pickle.load(open(os.path.join(outpath + 'capon_bf'), "rb"))

    plot_tomo_slices(np.absolute(Capon_bf))

if func_plot_profiles == 1:
    print('------------------------------------------------------------')
    print('Plot reflectivity profile')
    print('------------------------------------------------------------')
    # restoring te variable
    Capon_bf = pickle.load(open(os.path.join(outpath + 'capon_bf'), "rb"))

    plot_profiles(np.absolute(Capon_bf), list_pixels)
# *****************************************************************************
# end of script
# *****************************************************************************
