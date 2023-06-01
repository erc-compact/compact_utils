#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse, Circle
from matplotlib import patheffects as PathEffects
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
import sys, argparse, copy, glob, random, os
import matplotlib.patches as mpatches
from mosaic.plot import get_fwhm_of_telescope

def parse_argument():

    for i, arg in enumerate(sys.argv):
        if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg
    parser = argparse.ArgumentParser()

    parser.add_argument('--tiling_plot', nargs=1, metavar="file", help='filename for the tiling plot')
    parser.add_argument('--tiling_files', nargs=1, metavar="file", help='Comma separated list of files for the boresight, shape and coordinates of the beams')
    parser.add_argument('--beam_size_scaling', nargs=1, metavar="scaling", help='scaling factor for the size of the beam')
    parser.add_argument("--flip", action="store_true", help='flip the orientation of the beam')
    parser.add_argument("--index", action="store_true", help='wether to add index to beams')
    parser.add_argument("--inner", nargs=1, metavar="number", help='highlight the most [number] inner beams')
    parser.add_argument("--extra_source", nargs=1, metavar="file", help='extra point sources to plot')
    parser.add_argument("--core", help="Core radius of GC in deg", type=float)
    parser.add_argument("--half_light", help="Half light radius of GC in deg", type=float)
    parser.add_argument("--half_mass", help="Half mass radius of GC in deg", type=float)


    args = parser.parse_args()

    return args

args = parse_argument()

if args.beam_size_scaling is not None:
    scaling = float(args.beam_size_scaling[0])
else:
    scaling = 1.0

fileName = args.tiling_plot[0]

if args.inner is not None:
    inner = int(args.inner[0])
else:
    inner = 0

index = args.index

step = 1/10000000000.
wcs_properties = wcs.WCS(naxis=2)
wcs_properties.wcs.crpix = [0, 0]
wcs_properties.wcs.cdelt = [-step, step]
wcs_properties.wcs.ctype = ["RA---TAN", "DEC--TAN"]
resolution = step
core_radius = args.core
half_light_radius = args.half_light
half_mass_radius = args.half_mass

fig = plt.figure(figsize=(3000./96, 3000./96), dpi=96)
axis = fig.add_subplot(111,aspect='equal', projection=wcs_properties)
zoom_out_axis = fig.add_axes([0.745, 0.06, 0.2, 0.2], projection=wcs_properties)

color_palette = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']
gbt_dish_diameter = 100
gbt_freq = 1.4e9
meerkat_dish_diameter = 13.5
gbt_beam_size = get_fwhm_of_telescope(gbt_dish_diameter, gbt_freq)
gbt_beam_size_radius = gbt_beam_size/2
meerkat_beam_size = get_fwhm_of_telescope(meerkat_dish_diameter, gbt_freq)
meerkat_beam_size_radius = meerkat_beam_size/2
gbt_beam = Circle(xy=(0,0), radius = gbt_beam_size_radius/step, linestyle='dotted', linewidth=4, fill=False, label='GBT L-BAND Primary beam', color='g')
copy_gbt_beam = Circle(xy=(0,0), radius = gbt_beam_size_radius/step, linestyle='dotted', linewidth=4, fill=False, label='GBT L-BAND Primary beam', color='g')
merkat_beam = Circle(xy=(0,0), radius = meerkat_beam_size_radius/step, linestyle='dotted', linewidth=4, fill=False, label='MeerKAT L-BAND Primary beam', color='r')
copy_merkat_beam = Circle(xy=(0,0), radius = meerkat_beam_size_radius/step, linestyle='dotted', linewidth=4, fill=False, label='MeerKAT L-BAND Primary beam', color='r')

axis.add_artist(gbt_beam)
#axis.add_artist(merkat_beam)
circle_core = Circle(xy=(0,0), radius = core_radius/step, linestyle='--', linewidth=4, fill=False, label='Core')
circle_half_light = Circle(xy=(0,0), radius = half_light_radius/step, linestyle='-.', linewidth=4, fill=False, label='Half light')
circle_half_mass = Circle(xy=(0,0), radius = half_mass_radius/step, linestyle=':', linewidth=4, fill=False, label='Half mass')

copy_circle_core = Circle(xy=(0,0), radius = core_radius/step, linestyle='--', linewidth=0.5, fill=False, label='Core')
copy_circle_half_light = Circle(xy=(0,0), radius = half_light_radius/step, linestyle='-.', linewidth=0.5, fill=False, label='Half light')
copy_circle_half_mass = Circle(xy=(0,0), radius = half_mass_radius/step, linestyle=':', linewidth=0.5, fill=False, label='Half mass')

axis.add_artist(circle_core)
axis.add_artist(circle_half_light)
axis.add_artist(circle_half_mass)
beam_color_legend = []

file_paths = []

for file_names in args.tiling_files:
    file_names_list = file_names.split(",")
    for file_name in file_names_list:
        if not os.path.isfile(file_name):
            print("File does not exist: " + file_name, "Exiting...")
            sys.exit()
        file_paths.extend(glob.glob(file_name))

for i, file_path in enumerate(file_paths):
    
    # Process each file as needed
    df = pd.read_csv(file_path, nrows=10, header=None)
    new_tiling_metadata = dict(zip(df[0].str.split(': ').str[0], df[0].str.split(': ').str[1]))
    boresight_new_tiling = str(new_tiling_metadata['Boresight'])
    beamshape_x_new_tiling = float(new_tiling_metadata['Beam_x_deg'])
    beamshape_y_new_tiling = float(new_tiling_metadata['Beam_y_deg'])
    beamshape_angle_new_tiling = float(new_tiling_metadata['Beam_angle'])
    date = str(new_tiling_metadata['DATE'])
    project = str(new_tiling_metadata['PROJECT'])
    cluster = str(new_tiling_metadata['CLUSTER'])
    freq_band = str(new_tiling_metadata['BAND'])
    overlap = float(new_tiling_metadata['Overlap'])
   

    beam_coords = pd.read_csv(file_path, skiprows=10, delim_whitespace=True)
    BEAM_RA = beam_coords['RA'].astype(str)
    BEAM_DEC = beam_coords['DEC'].astype(str)
    BEAM_NAMES = beam_coords['Beam'].astype(str)

    equatorialCoordinates = SkyCoord(BEAM_RA, BEAM_DEC, frame='fk5', unit=(u.hourangle, u.deg))

    equatorialCoordinates = np.array([equatorialCoordinates.ra.astype(float), equatorialCoordinates.dec.astype(float)]).T

    if args.flip:
        beam_angle = 180 - beamshape_angle_new_tiling
        axis1, axis2, angle = (beamshape_x_new_tiling * scaling, beamshape_y_new_tiling * scaling, beam_angle)  
    else:
        beam_angle = beamshape_angle_new_tiling
        axis1, axis2, angle = (beamshape_x_new_tiling * scaling, beamshape_y_new_tiling * scaling, beam_angle)
    boresight_ra, boresight_dec = boresight_new_tiling.split()
    equatorialBoresight = SkyCoord(boresight_ra, boresight_dec, frame='fk5', unit=(u.hourangle, u.deg))
    boresight = (equatorialBoresight.ra.deg , equatorialBoresight.dec.deg)
    wcs_properties.wcs.crval = boresight
    center = boresight
    inner_idx = []

    scaled_pixel_coordinats = wcs_properties.wcs_world2pix(equatorialCoordinates, 0)
    beam_coordinate = np.array(scaled_pixel_coordinats)

    if inner > 0:
        index_sort = np.argsort(np.sum(np.square(beam_coordinate), axis=1))
        beam_coordinate = beam_coordinate.take(index_sort, axis=0)
        indice = indice.take(index_sort, axis=0)

    axis.plot(beam_coordinate[0][0], beam_coordinate[0][0], marker='+', markersize=15, color='black')

    r = random.random()
    g = random.random()
    b = random.random()
    edge_color = (r, g, b)  # Tuple of RGB value

    for idx in range(len(beam_coordinate)):
        coord = beam_coordinate[idx]
        if index == True:
            num = indice[idx].split('cfbf')[-1]
            axis.text(coord[0], coord[1], int(num), size=6, ha='center', va='center', color="gray", alpha=0.5)
        ellipse = Ellipse(xy=coord,
                width=2.*axis1/resolution,height=2.*axis2/resolution, angle=angle)
        ellipse.fill = False
        if inner > 0 and idx < inner:
            if index == True:
                inner_idx.append(int(num))
            ellipse.set(edgecolor = '#0066ff')
        else:
            # Select a contrasting color from the color palette
            edge_color = color_palette[i % len(color_palette)]
            ellipse.set(edgecolor = edge_color)
        

        axis.add_artist(ellipse)

    # Create a legend patch for the ellipse color
    #label = project + "_" + cluster + "_" + freq_band + "_" + date.replace("-", "/")
    label = project + "_" + cluster + "_" + freq_band + "_" + "overlap_" + str(overlap)

    #label = cluster + "_" + freq_band + "_" + "overlap_" + str(overlap)

    
    legend_patch = mpatches.Patch(color=edge_color, label=label)
    beam_color_legend.append(legend_patch)

     
    #margin = 1.5 * max(np.sqrt(np.sum(np.square(beam_coordinate), axis=1)))
    margin = 1.5 * max(np.sqrt(np.sum(np.square(beam_coordinate), axis=1)))
    axis.set_xlim(center[0]-margin, center[0]+margin)
    axis.set_ylim(center[1]-margin, center[1]+margin)


    # Plot the zoomed out view
    for idx in range(len(beam_coordinate)):
        coord = beam_coordinate[idx]
        ellipse = Ellipse(xy=coord,
                width=2.*axis1/resolution,height=2.*axis2/resolution, angle=angle)
        ellipse.fill = False
        ellipse.set_alpha(0.4)
        if inner > 0 and idx < inner:
            ellipse.set(edgecolor = '#0066ff')
        else:
            ellipse.set(edgecolor = edge_color)

        zoom_out_axis.add_artist(ellipse)
    zoom_out_axis.add_artist(copy_circle_core)
    zoom_out_axis.add_artist(copy_circle_half_light)
    zoom_out_axis.add_artist(copy_circle_half_mass)
    zoom_out_axis.add_artist(copy_gbt_beam)
    #zoom_out_axis.add_artist(copy_merkat_beam)
    margin2 = 3.0 * margin
    zoom_out_axis.set_xlim(center[0]-margin2, center[0]+margin2)
    zoom_out_axis.set_ylim(center[1]-margin2, center[1]+margin2)
    zoom_out_axis.patch.set_alpha(0.8)
    ra_zoom_out = zoom_out_axis.coords[0]
    dec_zoom_out = zoom_out_axis.coords[1]
    ra_zoom_out.tick_params(axis='x', which='both', bottom=False, top=False,
            labelbottom=False, right=False, left=False, labelleft=False)
    dec_zoom_out.tick_params(axis='y', which='both', bottom=False, top=False,
            labelbottom=False, right=False, left=False, labelleft=False)
    for spine in zoom_out_axis.spines.values():
        spine.set_edgecolor('grey')

legend_handles=[circle_core, circle_half_light, circle_half_mass, gbt_beam]
for idx in range(len(beam_color_legend)):
    legend_handles.append(beam_color_legend[idx])
axis.legend(handles = legend_handles, prop={'size': 34})


if args.extra_source is not None:
    extra_coords = np.genfromtxt(args.extra_source[0], dtype=None)
    if len(extra_coords.shape) == 1:
        extra_coords = extra_coords.reshape(1, -1)
    extra_equatorial_coordinates = SkyCoord(extra_coords[:,1].astype(str),
            extra_coords[:,2].astype(str),
        frame='fk5', unit=(u.hourangle, u.deg))
    extra_equatorial_coordinates = np.array([extra_equatorial_coordinates.ra.astype(float),
            extra_equatorial_coordinates.dec.astype(float)]).T
    scaled_extra_pixel_coordinats = np.array(wcs_properties.wcs_world2pix(
                extra_equatorial_coordinates, 0))
    source_ransom = axis.scatter(scaled_extra_pixel_coordinats[0,0], scaled_extra_pixel_coordinats[0,1],s=90, color='b', label="Ransom et al. 2004")
    source_trapum = axis.scatter(scaled_extra_pixel_coordinats[1,0], scaled_extra_pixel_coordinats[1,1],s=90, color='k', label="TRAPUM (this work)")

    for ext_coord_idx in range(len(extra_coords[0:5])):
            # font_weight = 'bold' if ext_coord_idx > 1 else 'regular'
            font_weight = 'bold'
            axis.annotate(str(extra_coords[ext_coord_idx][0])[2:-1],
            xy=scaled_extra_pixel_coordinats[ext_coord_idx, :],
            textcoords='offset pixels', xytext=(5,-3),
            size=80 , fontweight=font_weight, color='black', path_effects=[PathEffects.withStroke(linewidth=4, foreground="black")])

    for ext_coord_idx in range(len(extra_coords[5:18])):
        font_weight = 'bold' if ext_coord_idx > 5 else 'regular'
        
        axis.annotate(str(extra_coords[ext_coord_idx + 5][0])[2:-1],
                xy=scaled_extra_pixel_coordinats[ext_coord_idx + 5, :],
                textcoords='offset pixels', xytext=(5,-3),
                size=40 , fontweight=font_weight)

ra = axis.coords[0]
dec = axis.coords[1]
ra.set_ticklabel(size=20)
dec.set_ticklabel(size=20, rotation="vertical")
dec.set_ticks_position('l')
ra.set_ticks_position('b')
ra.set_axislabel("Right ascension (RA)", size=40)
dec.set_axislabel("Declination (DEC)", size=40)
#plt.tight_layout()
plt.subplots_adjust(left=0.01, bottom=0.05, right=0.999, top=0.95)
plt.suptitle('M30 LBAND Tiling: TRAPUM Vs COMPACT', size=55)
#plt.suptitle('M30 LBAND Tiling: COMPACT 90 vs 95 % overlap', size=45)

plt.savefig(fileName, dpi=96)    
#plt.show()
    
sys.exit()

























