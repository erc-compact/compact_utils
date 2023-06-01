#!/usr/bin/env python

import numpy as np
import sys, datetime, os
import argparse
import logging
import seaborn as sns
import matplotlib
# matplotlib.use('Agg')

import matplotlib.pyplot as plt
loggerFormat = '%(asctime)-15s %(filename)s  %(message)s'
logging.basicConfig(format = loggerFormat, level=logging.WARNING)
logger = logging.getLogger()

import katpoint
from mosaic.beamforming import PsfSim, generate_nbeams_tiling
from mosaic.coordinate import createTilingRegion, readPolygonRegion, convert_sexagesimal_to_degree, convert_equatorial_coordinate_to_pixel
from mosaic import __version__
from astropy.coordinates import SkyCoord
from astropy import units as u
from mosaic.plot import get_fwhm_of_telescope
import matplotlib.gridspec as gridspec
import pandas as pd
import matplotlib.patches as mpatches

# Generate color palette
colors = sns.color_palette("muted", 6)

# Define color mapping
overlap_colors = {
    "0.9": colors[0],
    "0.95": colors[1],
    "0.99": colors[2],
    "0.8": colors[3],
    "0.6": colors[4],
    "0.5": colors[5]
}




def plot_nbeams_versus_antennas(data, output_dir):
    data['nantennas_requested'] = data['nantennas_requested'].astype(str)
    cluster = data['cluster'].iloc[0]
    freq_band = data['band'].iloc[0]
    fig = plt.figure(figsize=(24, 14))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1,1]) 
    ax1 = plt.subplot(gs[0, 0])  # first row, first column
    ax2 = plt.subplot(gs[0, 1])  # first row, second column
    ax3 = plt.subplot(gs[1, :])  #
    gs.update(wspace=0.35, hspace=0.35)
    overlap_order = ["0.9", "0.95", "0.99", "0.8", "0.6", "0.5"]
    hue_order_top_row = ["0.9", "0.95", "0.99"]
    hue_order_bottom_row = ["0.9", "0.8", "0.6", "0.5"]

    plot_data_half_light = data.loc[(data['region'] == 'half_light') & ((data['overlap'] == 0.9) | (data['overlap'] == 0.95) | (data['overlap'] == 0.99))]
    plot_data_half_light['overlap'] = plot_data_half_light['overlap'].astype(str)
    plot_data_core = data.loc[(data['region'] == 'core') & ((data['overlap'] == 0.9) | (data['overlap'] == 0.95) | (data['overlap'] == 0.99))]
    plot_data_core['overlap'] = plot_data_core['overlap'].astype(str)
    plot_data_twice_half_light = data.loc[(data['region'] == 'twice_half_light_radius') & ((data['overlap'] == 0.9) | (data['overlap'] == 0.8) | (data['overlap'] == 0.6) | (data['overlap'] == 0.5))]
    plot_data_twice_half_light['overlap'] = plot_data_twice_half_light['overlap'].astype(str)
    sns.barplot(x = 'nantennas_requested', y = 'beams_required', hue = 'overlap', data = plot_data_half_light, ax = ax1, hue_order=hue_order_top_row, palette=overlap_colors)
    sns.barplot(x = 'nantennas_requested', y = 'beams_required', hue = 'overlap', data = plot_data_core, ax = ax2, hue_order=hue_order_top_row, palette=overlap_colors)
    sns.barplot(x = 'nantennas_requested', y = 'beams_required', hue = 'overlap', data = plot_data_twice_half_light, ax = ax3, hue_order=hue_order_bottom_row, palette=overlap_colors)

    # Increase this value to increase space between the bars and their labels
    vertical_space = 5

    for p in ax1.patches:
        height = p.get_height()
        if not np.isnan(height):
            ax1.text(p.get_x() + p.get_width() / 2., height + vertical_space, f'{int(height)}', ha='center', fontsize=18)

    vertical_space = 0.2

    for p in ax2.patches:
        height = p.get_height()
        if not np.isnan(height):
            ax2.text(p.get_x() + p.get_width() / 2., height + vertical_space, f'{int(height)}', ha='center', fontsize=18)

    vertical_space = 20

    for p in ax3.patches:
        height = p.get_height()
        if not np.isnan(height):
            ax3.text(p.get_x() + p.get_width() / 2., height + vertical_space, f'{int(height)}', ha='center', fontsize=20)

    ax1.set_title('Region: Max.(Half-light/Half-Mass)', fontsize=26, fontweight='bold', y=1.05)
    ax2.set_title('Region: Core', fontsize=26, fontweight='bold', y=1.05)
    ax3.set_title('Region: Twice Max.(Half-light/Half-Mass)', fontsize=26, fontweight='bold')

    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax3.set_yscale('log')
    ax1.set_ylabel('Beams Required', fontsize=26)
    ax2.set_ylabel('Beams Required', fontsize=26)
    ax3.set_ylabel('Beams Required', fontsize=26)
    ax1.set_xlabel('Number of Antennas', fontsize=26)
    ax2.set_xlabel('Number of Antennas', fontsize=26)
    ax3.set_xlabel('Number of Antennas', fontsize=26)

    #Increase tick label size
    ax1.tick_params(axis='both', which='major', labelsize=22)
    ax2.tick_params(axis='both', which='major', labelsize=22)
    ax3.tick_params(axis='both', which='major', labelsize=22)
    plt.suptitle(f"COMPACT: BEAMS vs NANTENNAS - {cluster}: {freq_band}", fontsize=30, fontweight='bold')
    
    overlap_legend_1 = ["0.9", "0.95", "0.99"]
    overlap_legend_2 = ["0.9", "0.95", "0.99"]
    overlap_legend_3 = ["0.9", "0.8", "0.6", "0.5"]

    ax1.legend([mpatches.Patch(color=overlap_colors[key]) for key in overlap_legend_1],
           overlap_legend_1, title='Overlap fraction', fontsize=16, title_fontsize=18)

    ax2.legend([mpatches.Patch(color=overlap_colors[key]) for key in overlap_legend_2],
           overlap_legend_2, title='Overlap fraction', fontsize=16, title_fontsize=18)

    ax3.legend([mpatches.Patch(color=overlap_colors[key]) for key in overlap_legend_3],
           overlap_legend_3, title='Overlap fraction', fontsize=22, title_fontsize=18)
    
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir + '/' + f"beams_vs_nantennas_{cluster}_{freq_band}.png", dpi=300)




def makeKatPointAntenna(antennaString):
    antennaKat = []

    for antenna in antennaString:
        antkat = katpoint.Antenna(antenna)
        antennaKat.append(antkat)

    return antennaKat

def meerkat_band(central_freq):
    if central_freq > 2100e6:
        return 'SBAND'
    elif central_freq > 1200e6:
        return 'LBAND'
    elif central_freq > 800e6:
        return 'UHF'

def grid_search_required_beams(antennaCoords, sourceCoord, observeTime, frequencies,
        duration, overlap, subarray, update_interval, overlay_source, overlay_source_name,
        size, resolution, tilingMethod, tilingShape, tilingParameter, tilingParameterCoordinateType,
        weights, interpolation, output, core_radius, half_light_radius, half_mass_radius):
    """
    This function is used to find the number of beams required 
    to cover the region of interest using a binary search algorithm.
    """
  
    if subarray != []:
        antennaKat = makeKatPointAntenna(
                [antennaCoords[ant] for ant in subarray])
    else:
        antennaKat = makeKatPointAntenna(antennaCoords)

    
    # boresight = sourceCoord
    boresight = katpoint.Target('boresight, radec, {}, {}'.format(
                sourceCoord[0], sourceCoord[1]))
    reference = makeKatPointAntenna(
            ["ref, -30:42:39.8, 21:26:38.0, 1035.0",])[0]
    psf = PsfSim(antennaKat, frequencies[0])

    newBeamShape = psf.get_beam_shape(sourceCoord, observeTime, size, resolution, weights)
    
    meerkat_dish = 13.5
    meerkat_fwhm = get_fwhm_of_telescope(meerkat_dish, frequencies[0])
    meerkat_radius = meerkat_fwhm/2

    gbt_dish = 100
    gbt_fwhm = get_fwhm_of_telescope(gbt_dish, frequencies[0])
    gbt_radius = gbt_fwhm/2
    # Giving a higher upper threshold to be safe and for a more denser packing.
    upper_threshold = 0.4
    lower_threshold = 0.02
    
    # Initialize the beam number bounds
    min_beamNum, max_beamNum = 1, 100000 # replace 10000 with a suitable upper limit
    half_light_radius = max(half_light_radius, half_mass_radius)
    beamNums = {"core": None, "half_light": None, "twice_half_light_radius": None}
    twice_half_light_radius = 2 * half_light_radius

    for region, radius in [("core", core_radius), ("half_light", half_light_radius), ("twice_half_light_radius", twice_half_light_radius)]:
        # Simple binary search implementation to quickly find the number of beams within a 5 percent threshold of the desired radius
        while min_beamNum <= max_beamNum:
            beamNum = (min_beamNum + max_beamNum) // 2
            
            tiling = generate_nbeams_tiling(newBeamShape, beamNum, overlap, tilingMethod, tilingShape, core_radius, half_light_radius, half_mass_radius, frequencies,
                                            tilingParameter, tilingParameterCoordinateType)
            
            beam_x_deg, beam_y_deg, beam_angle = tiling.meta["axis"][:3]
            
            equatorial_coordinates = tiling.get_equatorial_coordinates()
            equatorial_coordinates = SkyCoord(ra=equatorial_coordinates[:,0]*u.deg, dec=equatorial_coordinates[:,1]*u.deg)
            ra = equatorial_coordinates.ra.to_string(unit=u.hour, sep=':', pad=True, precision=2)
            dec = equatorial_coordinates.dec.to_string(unit=u.deg, sep=':', pad=True, precision=2)
            band = meerkat_band(frequencies[0])
            cluster = output["pointing_name"]
            boresight_ra = sourceCoord[0]
            boresight_dec = sourceCoord[1]
            boresight_position_string = f"{boresight_ra}{boresight_dec}"
            boresight = SkyCoord(ra=boresight_ra, dec=boresight_dec, unit=(u.hourangle, u.deg))
            #Calculate beam furthest from boresight in degrees
            distance_to_each_beam = np.sqrt((equatorial_coordinates.ra-boresight.ra)**2 + (equatorial_coordinates.dec-boresight.dec)**2)
            max_beam_distance = np.max(distance_to_each_beam)
            max_beam_distance = max_beam_distance.to(u.deg).value

            # Special Case when one beam covers the entire region
            semi_major_axis = max(beam_x_deg, beam_y_deg)
            
            if semi_major_axis > radius:
                beamNums[region] = 1
                break
            

            elif (1 - lower_threshold) * radius <= max_beam_distance <= (1 + upper_threshold) * radius:
                beamNums[region] = beamNum
                break
        
            elif max_beam_distance < (1 - lower_threshold) * radius:
               
                min_beamNum = beamNum + 1
            
            else: 
               
                max_beamNum = beamNum - 1
               

        #Corner case when min_beamNum >= max_beamNum. This happens due to our large threshold and oscillations. Just take the max value
        if beamNums[region] is None:
            print(f"Beam Number = {beamNum}, region = {region}, too few beams")
            beamNums[region] = max_beamNum
      
        min_beamNum, max_beamNum = 1, 100000

    for region, beamNum in beamNums.items():
        
        tiling = generate_nbeams_tiling(newBeamShape, beamNum, overlap, tilingMethod, tilingShape, core_radius, half_light_radius, half_mass_radius, frequencies,
                                        tilingParameter, tilingParameterCoordinateType)
        beam_x_deg, beam_y_deg, beam_angle = tiling.meta["axis"][:3]
        generated_beams_mosiac = len(tiling.coordinates)
        print(f"Generated {generated_beams_mosiac} beams for {region} region")

        if subarray != []:
            nantennas_requested = len(subarray)
        else:
            nantennas_requested = 64

        if not os.path.isfile("beams_required_versus_antennas_grid.csv"):
            with open ("beams_required_versus_antennas_grid.csv", 'w') as f:
                f.write("cluster,band,central_frequency,boresight_position,beam_x_deg,beam_y_deg,beam_angle,date,nantennas_requested,overlap,region,beams_required\n")

        with open ("beams_required_versus_antennas_grid.csv", 'a') as f:
            f.write(f"{cluster},{band},{frequencies[0]},{boresight_position_string},{beam_x_deg},{beam_y_deg},{beam_angle},{observeTime},{nantennas_requested},{overlap},{region},{generated_beams_mosiac}\n")

        if "tiling_plot" in output:
            overlap_string = str(overlap).replace(".", "_")
            plot_filename = output["tiling_plot"][0].replace(".png", f"_{band}_{region}_{nantennas_requested}_{overlap_string}.png")  # modify filename
            if generated_beams_mosiac < 10000:
                tiling.plot_tiling(plot_filename, HD=True, edge=True, extra_coordinates = overlay_source, extra_coordinates_text = overlay_source_name)


    
def parseOptions(parser):
    # enable interpolation when ploting
    parser.add_argument('--no_interpolation', action='store_true', help=argparse.SUPPRESS)
    # plot the beam images
    parser.add_argument('--psf_plot', nargs='+', metavar="file [paras]", help='filename for the psf plot')
    parser.add_argument('--psf_fits', nargs='+', metavar="file [paras]", help='name for the psf fits file')
    parser.add_argument('--tiling_plot', nargs='+', metavar="file [paras]", help='filename for the tiling plot')
    parser.add_argument('--tiling_coordinate', nargs='+', metavar="file [paras]", help='filename for the tiling coordinate')
    parser.add_argument('--tiling_region', nargs=1, metavar="file", help='filename for the tiling region')
    parser.add_argument('--ants', nargs=1, metavar="file", help='antenna coodinates files')
    parser.add_argument('--resolution', nargs=1, metavar="asec", help='resolution in arcsecond')
    parser.add_argument('--size', nargs=1, metavar="num", help='width in pixels')
    parser.add_argument('--freqrange', nargs=3, metavar=('s', 't', 'i'), help='freqencies range as start stop interval')
    parser.add_argument('--freq', nargs='+',  help='multiple freqencies')
    parser.add_argument('--frame', nargs=1, metavar="RADEC/AziAlt", help='source coordinate frame')
    parser.add_argument('--source', nargs=2, metavar=('RA/Azi', 'DEC/Alt'), help='source position in RADEC or AziAlt')
    parser.add_argument('--datetime', nargs=2, metavar=('date', 'time'), help='observation time in format: 03/10/2015 15:23:10.000001')
    parser.add_argument('--overlap', nargs=1, metavar="ratio", help='overlap point between beams')
    parser.add_argument('--overlay_source', nargs=1, metavar="file", help='extra overlay sources')
    parser.add_argument('--subarray', nargs='+', metavar="num", help='list of antennas, saperated by comma')
    parser.add_argument('--subarray_file', nargs=1, metavar="file", help='File containing list of antennas, in each line')
    parser.add_argument('--tiling_method', nargs=1, metavar="method",
            help='tiling method, such as \"variable_overlap\" or \"variable_size\".')
    parser.add_argument('--tiling_shape', nargs=1, metavar="shape",
            help='shape of the tiling boundary, such as \"circle\", \"hexagon\", \"ellipse\", \"polygon\", \"annulus\".')
    parser.add_argument('--tiling_parameter', nargs='+', metavar='parameters', help='parameters for the tiling')
    parser.add_argument('--tiling_parameter_file', nargs=1, metavar='parameter_file', help='parameter_file for the tiling')
    parser.add_argument('--tiling_parameter_coordinate_type', nargs=1, metavar='coordinate_type',
            help='type of the coordinate of the tiling parameter, such as \"equatorial\", default is image(pixel) coordinate')
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument("--version", help="show the version of this package", action="store_true")
    parser.add_argument("--weight", action="store_true",
            help='apply weights to individual antenna, attach weight after the item in --subarray, e.g., 0:0.5, 1:0.7, 2:0.5 ')

    parser.add_argument("--core", help="Core radius of GC in deg", type=float)
    parser.add_argument("--half_light", help="Half light radius of GC in deg", type=float)
    parser.add_argument("--half_mass", help="Half mass radius of GC in deg", type=float)
    parser.add_argument("--pointing_name", help="Name of pointing usually GC name", type=str)
    parser.add_argument("--overlap_range", nargs='+', help="Single or csv list of overlap fraction between beams", type=str)
    parser.add_argument('--subarray_files', nargs='+', metavar="file", help='Files containing list of antennas, in each line')
    parser.add_argument('--plot_nbeams_antennas', help='Plot Nbeams vs antenna plot. Default: False', action='store_true')


    args = parser.parse_args()
   

   
    if args.verbose:
        logger.setLevel(logging.INFO)

    if args.version is True:
        print(__version__)
        exit(0)
    else:
        logger.info("Mosaic " + __version__)

    if args.weight is True:
        weights = []
    else:
        weights = None
    

    frequencies = [1.4e9,]
    plotting_data = []
    paras = None
    resolution = None #in arcsecond
    size = 400
    duration = 0
    overlap = 0.5
    update_interval = 3600
    overlay_source = None
    output = {}
    beam_antenna_output_dir = 'BEAMS_NANTENNA_PLOTS'
    if args.psf_plot is not None:
        output["psf_plot"] = args.psf_plot
    if args.psf_fits is not None:
        output["psf_fits"] = args.psf_fits
    if args.tiling_plot is not None:
        output["tiling_plot"] = args.tiling_plot
    if args.tiling_coordinate is not None:
        output["tiling_coordinate"] = args.tiling_coordinate
    if args.tiling_region is not None:
        output["tiling_region"] = args.tiling_region[0]
    if args.pointing_name is not None:
        output["pointing_name"] = args.pointing_name
    if args.plot_nbeams_antennas:
        output["plot_nbeams_antennas"] = 'True'
        output["plot_nbeams_antennas_dir"] = beam_antenna_output_dir
 
    if args.ants is not None:
        with open(args.ants[0], 'r') as antFile:
            antennaCoords = antFile.readlines()
    else:
        parser.error("no antennas file, try --ants file")

    if args.datetime is not None:
        try:
            observeTime=datetime.datetime.strptime(args.datetime[0] + " "
                    + args.datetime[1], '%Y.%m.%d %H:%M:%S.%f')
        except: 
            observeTime=datetime.datetime.strptime(args.datetime[0] + " "
                    + args.datetime[1], '%Y-%m-%d %H:%M:%S.%f')
    else:
        parser.error("no time specified, try --datetime date time")
    
    #try:
    if args.subarray is not None:
        subarray = []
        arrayString = "".join(args.subarray)
        ant_weights = arrayString.split(",")
        for ant_weight in ant_weights:
            ant_weight_pair = ant_weight.split(':')
            if weights is not None:
                subarray.append(int(ant_weight_pair[0]))
                if len(ant_weight_pair) > 1:
                    complex_weight = complex(ant_weight_pair[1])
                    if complex_weight.imag == 0:
                        weights.append(float(ant_weight_pair[1]))
                    else:
                        weights.append(complex_weight)
                else:
                    weights.append(1.0)
            else:
                subarray.append(int(ant_weight_pair[0]))
    else:
        subarray = []


    if args.source is not None:
        sourceCoord = args.source
    else:
        parser.error("no source specifed, try --source RA DEC")

    if args.overlap is not None:
        overlap = float(args.overlap[0])
    else:
        overlap = 0.5

    if args.resolution is not None:
        resolution = float(args.resolution[0])
    if args.size is not None:
        size = int(args.size[0])
    if args.freqrange is None and args.freq is None:
        parser.error("no freqencies or frequency range specified")
    elif args.freq is not None:
        frequencies = [float(freq) for freq in args.freq]
    elif args.freqrange is not None:
        frequencies = np.arange(float(args.freqrange[0]),
                float(args.freqrange[1]), float(args.freqrange[2]))

    if args.tiling_shape is not None:
        tilingShape = args.tiling_shape[0]
    else:
        tilingShape = "circle"

    if args.tiling_parameter_coordinate_type is not None:
        tilingParameterCoordinateType = args.tiling_parameter_coordinate_type[0]
    else:
        tilingParameterCoordinateType  = 'equatorial'

    if args.tiling_method is not None:
        tilingMethod = args.tiling_method[0]
        if args.tiling_method[0] == "variable_overlap":
            if args.tiling_parameter is None and tilingShape != "polygon":
                parser.error("no parameter for \"variable_overlap\" method!")
                exit(-1)
            if tilingShape == "circle":
                tilingParameter = float(args.tiling_parameter[0])
            elif tilingShape == "hexagon":
                tilingParameter = [float(args.tiling_parameter[0]),
                                    float(args.tiling_parameter[1])]
            elif tilingShape == "ellipse":
                tilingParameter = [float(args.tiling_parameter[0]),
                                    float(args.tiling_parameter[1]),
                                    float(args.tiling_parameter[2])]
            elif tilingShape == "polygon":
                if args.tiling_parameter is not None:
                    parameterString = "".join(args.tiling_parameter).split(",")
                    tilingParameter = [float(coord) for coord in parameterString]
                    tilingParameter = np.array(tilingParameter).reshape(-1,2).tolist()
                elif args.tiling_parameter_file is not None:
                    tilingParameter = readPolygonRegion(
                            args.tiling_parameter_file[0]).tolist()
                else:
                    parser.error("no parameter for polygon tiling!")
                    exit(-1)
        else:
            tilingParameter = None
    else:
        tilingMethod = "variable_size"
        tilingParameter = None

    if args.no_interpolation:
        interpolation = False
    else:
        interpolation = True

    if args.overlay_source is not None:
        overlay_coords = np.genfromtxt(args.overlay_source[0], dtype=None)
        if len(overlay_coords.shape) == 1:
            overlay_coords = overlay_coords.reshape(1, -1)
        bore_sight = convert_sexagesimal_to_degree([sourceCoord,])[0]
        overlay_source_degree = convert_sexagesimal_to_degree(overlay_coords[:, 1:])
        overlay_coordinates = convert_equatorial_coordinate_to_pixel(
                overlay_source_degree, bore_sight)

        overlay_source = overlay_coordinates
        overlay_source_name = overlay_coords[:, 0].astype('str')

    else:
        overlay_source = []
        overlay_source_name = []
   
    if args.core:
        core_radius = args.core
    else:
        parser.error("no core radius specified, try --core radius")

    if args.half_light:
        half_light_radius = args.half_light
    else:
        parser.error("no half light radius specified, try --half_light radius")

    if args.half_mass:
        half_mass_radius = args.half_mass
    else:
        parser.error("no half mass radius specified, try --half_light radius")

    if args.overlap_range is not None:
        # Check if a single value or multiple values are provided
        if len(args.overlap_range) == 1:
            # Single string value possibly with commas
            overlap_range = [float(value) for value in args.overlap_range[0].split(',')]
        else:
            # Multiple strings values
            overlap_range = [float(value) for value in args.overlap_range]

    
    if args.subarray_files is not None:
        subarray_files = [str(file) for file in ','.join(args.subarray_files).split(',')]
    else:
        subarray_files = []
    band = meerkat_band(frequencies[0])
    
    for i in range(len(overlap_range)):
        for subarray_file in subarray_files:
            subarray = []
            print("Reading subarray file: %s" % subarray_file)
            try:
                with open(subarray_file, 'r') as f:
                    for line in f:
                        ant_weight = line.strip()
                        ant_weight_pair = ant_weight.split(':')
                        if weights is not None:
                            subarray.append(int(ant_weight_pair[0]))
                            if len(ant_weight_pair) > 1:
                                complex_weight = complex(ant_weight_pair[1])
                                if complex_weight.imag == 0:
                                    weights.append(float(ant_weight_pair[1]))
                                else:
                                    weights.append(complex_weight)
                            else:
                                weights.append(1.0)
                        else:
                            subarray.append(int(ant_weight_pair[0]))
                subarray.sort()
                paras  = {"antennaCoords": antennaCoords,
                "sourceCoord": sourceCoord,
                "observeTime": observeTime,
                "frequencies":frequencies,
                "duration":duration,
                "overlap":overlap_range[i],
                "subarray":subarray,
                "update_interval":update_interval,
                "size":size,
                "resolution":resolution,
                "tilingMethod":tilingMethod,
                "tilingShape":tilingShape,
                "tilingParameter":tilingParameter,
                "tilingParameterCoordinateType":tilingParameterCoordinateType,
                "overlay_source":overlay_source,
                "overlay_source_name":overlay_source_name,
                "weights":weights,
                "interpolation":interpolation,
                "output":output,
                "core_radius":core_radius,
                "half_light_radius":half_light_radius,
                "half_mass_radius":half_mass_radius}
                grid_search_required_beams(**paras)

            except FileNotFoundError:
                print(f"File {subarray_file} not found.")
    
    if "plot_nbeams_antennas" in output:
        df = pd.read_csv('beams_required_versus_antennas_grid.csv')
        # Select the cluster and right freq band
        df = df.loc[(df['cluster'] == output["pointing_name"]) & (df['band'] == band)]
        beam_antenna_output_dir = output["plot_nbeams_antennas_dir"]
        if not os.path.exists(beam_antenna_output_dir):
            os.makedirs(beam_antenna_output_dir)
        
        # Plot the data
        plot_nbeams_versus_antennas(df, beam_antenna_output_dir)

def captureNegetiveNumber():
    for i, arg in enumerate(sys.argv):
          if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

def main():
    captureNegetiveNumber()
    parser = argparse.ArgumentParser()
    parseOptions(parser)

if __name__ == "__main__":
    main()
