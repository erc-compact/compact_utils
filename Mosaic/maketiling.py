#!/usr/bin/env python

import numpy as np
import sys, datetime, os
import argparse
import logging

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

def createBeamMatrix(antennaCoords, sourceCoord, observeTime, frequencies, beamNum,
        duration, overlap, subarray, update_interval, overlay_source, overlay_source_name,
        size, resolution, tilingMethod, tilingShape, tilingParameter, tilingParameterCoordinateType,
        weights, interpolation, output, core_radius, half_light_radius, half_mass_radius, append_data_core, append_data_half_light, append_data_gbt):
  
    if subarray != []:
        antennaKat = makeKatPointAntenna(
                [antennaCoords[ant] for ant in subarray])
    else:
        antennaKat = makeKatPointAntenna(antennaCoords)

    boresight = katpoint.Target('boresight, radec, {}, {}'.format(
                sourceCoord[0], sourceCoord[1]))
    reference = makeKatPointAntenna(
            ["ref, -30:42:39.8, 21:26:38.0, 1035.0",])[0]
    psf = PsfSim(antennaKat, frequencies[0])
    newBeamShape = psf.get_beam_shape(sourceCoord, observeTime, size, resolution, weights)
    if "psf_plot" in output:
        newBeamShape.plot_psf(output["psf_plot"][0], overlap = overlap,
                shape_overlay=True, interpolation=interpolation)
    if "psf_fits" in output:
        newBeamShape.psf.write_fits(output["psf_fits"][0])
    if ("tiling_plot" in output) or ("tiling_coordinate" in output) or ("tiling_region" in output) or ("tiling_position_file" in output):
       
        tiling = generate_nbeams_tiling(newBeamShape, beamNum, overlap, tilingMethod, tilingShape, core_radius, half_light_radius, half_mass_radius, frequencies,
        tilingParameter, tilingParameterCoordinateType)
        beam_x_deg, beam_y_deg, beam_angle = tiling.meta["axis"][:3]
        #date = datetime.datetime.now().strftime("%d-%m-%Y")
        generated_beams_mosiac = len(tiling.coordinates)

        project = "COMPACT"
        if "tiling_plot" in output:
            tiling.plot_tiling(output["tiling_plot"][0], HD=True, edge=True,
                    extra_coordinates = overlay_source,
                    extra_coordinates_text = overlay_source_name)
        if "tiling_coordinate" in output:
            equatorial_coordinates = tiling.get_equatorial_coordinates()
            np.savetxt(output["tiling_coordinate"][0], equatorial_coordinates)
        if "tiling_region" in output:
            equatorial_coordinates = tiling.get_equatorial_coordinates()
            actualShape = tiling.meta["axis"]
            createTilingRegion(equatorial_coordinates, actualShape, output["tiling_region"])

        if "tiling_position_file" in output:
            print("Writing tiling position file to ", output["tiling_position_file"])
            equatorial_coordinates = tiling.get_equatorial_coordinates()
            equatorial_coordinates = SkyCoord(ra=equatorial_coordinates[:,0]*u.deg, dec=equatorial_coordinates[:,1]*u.deg)
            ra = equatorial_coordinates.ra.to_string(unit=u.hour, sep=':', pad=True, precision=2)
            dec = equatorial_coordinates.dec.to_string(unit=u.deg, sep=':', pad=True, precision=2)
            band = meerkat_band(frequencies[0])
            cluster = output["pointing_name"]
            boresight_ra = sourceCoord[0]
            boresight_dec = sourceCoord[1]
           
            with open(output["tiling_position_file"], 'w') as f:
                # Write the cluster, band, boresight, beam information, date, and project
                f.write(f"CLUSTER: {cluster}\n")
                f.write(f"BAND: {band}\n")
                f.write(f"Central_freq: {frequencies[0]}\n")
                f.write(f"Boresight: {boresight_ra} {boresight_dec}\n")
                f.write(f"Beam_x_deg: {beam_x_deg}\n")
                f.write(f"Beam_y_deg: {beam_y_deg}\n")
                f.write(f"Beam_angle: {beam_angle}\n")
                f.write(f"Overlap: {overlap}\n")
                f.write(f"DATE: {observeTime}\n")
                f.write(f"PROJECT: {project}\n\n")
                
                #  Write the header of the data
                f.write("Beam\tRA\tDEC\n")
                for i in range(len(ra)):
                    beam_name = 'cfbf{:04d}'.format(i+1)
                    f.write(f"{beam_name}\t{ra[i]}\t{dec[i]}\n")
        if subarray != []:
            nantennas_requested = len(subarray)
        else:
            nantennas_requested = 64
        boresight_ra, boresight_dec = sourceCoord
        boresight_position_string = f"{boresight_ra}{boresight_dec}"
        
        if not os.path.isfile("beams_required_versus_antennas.csv"):
            with open ("beams_required_versus_antennas.csv", 'w') as f:
                f.write("cluster,band,central_frequency,boresight_position,beam_x_deg,beam_y_deg,beam_angle,date,nantennas_requested,overlap,region,beams_required\n")
        if append_data_core:
            region = 'core'
            with open ("beams_required_versus_antennas.csv", 'a') as f:
                f.write(f"{cluster},{band},{frequencies[0]},{boresight_position_string},{beam_x_deg},{beam_y_deg},{beam_angle},{observeTime},{nantennas_requested},{overlap},{region},{generated_beams_mosiac}\n")
        if append_data_half_light:
            region = 'half_light'
            with open ("beams_required_versus_antennas.csv", 'a') as f:
                f.write(f"{cluster},{band},{frequencies[0]},{boresight_position_string},{beam_x_deg},{beam_y_deg},{beam_angle},{observeTime},{nantennas_requested},{overlap},{region},{generated_beams_mosiac}\n")
        if append_data_gbt:
            region = 'gbt'
            with open ("beams_required_versus_antennas.csv", 'a') as f:
                f.write(f"{cluster},{band},{frequencies[0]},{boresight_position_string},{beam_x_deg},{beam_y_deg},{beam_angle},{observeTime},{nantennas_requested},{overlap},{region},{generated_beams_mosiac}\n")


def parseOptions(parser):
    # enable interpolation when ploting
    parser.add_argument('--no_interpolation', action='store_true', help=argparse.SUPPRESS)
    # plot the beam images
    parser.add_argument('--psf_plot', nargs='+', metavar="file [paras]", help='filename for the psf plot')
    parser.add_argument('--psf_fits', nargs='+', metavar="file [paras]", help='name for the psf fits file')
    parser.add_argument('--tiling_plot', nargs='+', metavar="file [paras]", help='filename for the tiling plot')
    parser.add_argument('--tiling_coordinate', nargs='+', metavar="file [paras]", help='filename for the tiling coordinate')
    parser.add_argument('--tiling_position_file', nargs='+', metavar="file [paras]", help='filename for the tiling position for plotting')
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
    parser.add_argument('--beamnum', nargs=1, metavar="num", help='beam number of tiling')
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
    parser.add_argument("--append_data_core", help="Append n-antennas, overlap, generated beams required to cover core region for plotting later", action='store_true')
    parser.add_argument("--append_data_half_light", help="Append n-antennas, overlap, generated beams required to cover larger of half-mass/half light region for plotting later", action='store_true')
    parser.add_argument("--append_data_gbt", help="Append n-antennas, overlap, generated beams required to cover larger of GBT primary beam for plotting later", action='store_true')



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
    
    append_data_core = False
    append_data_half_light = False
    append_data_gbt = False
    if args.append_data_core:
        append_data_core = True
    if args.append_data_half_light:
        append_data_half_light = True
    if args.append_data_gbt:
        append_data_gbt = True
    

    frequencies = [1.4e9,]
    plotting_data = []
    paras = None
    resolution = None #in arcsecond
    size = 400
    duration = 0
    overlap = 0.5
    beamnum = 400
    update_interval = 3600
    overlay_source = None
    output = {}
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
    if args.tiling_position_file is not None:
        output["tiling_position_file"] = args.tiling_position_file[0]
    if args.pointing_name is not None:
        output["pointing_name"] = args.pointing_name
 
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
    
    try:
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
    except:
        if args.subarray_file is not None:
            with open(args.subarray_file[0], 'r') as subarrayFile:
                subarray = subarrayFile.readlines()

        print("No subarray specified, will read from file")

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
        
        print("Trying to read from file...")
        if args.subarray_file is not None:
            subarray_file = args.subarray_file[0]
            subarray = []
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
            except FileNotFoundError:
                print(f"File {args.subarray_file} not found.")
            
        else:
            print("No subarray file or list provided.")
            subarray = []
   
    if args.source is not None:
        sourceCoord = args.source
    else:
        parser.error("no source specifed, try --source RA DEC")

    if args.beamnum is not None:
        beamnum = int(args.beamnum[0])

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
        parser.error("no half mass radius specified, try --half_mass radius")

    paras  = {"antennaCoords": antennaCoords,
        "sourceCoord": sourceCoord,
        "observeTime": observeTime,
        "frequencies":frequencies,
        "beamNum":beamnum,
        "duration":duration,
        "overlap":overlap,
        "subarray":subarray,
        "update_interval":update_interval,
        "overlay_source":overlay_source,
        "overlay_source_name":overlay_source_name,
        "size":size,
        "resolution":resolution,
        "tilingMethod":tilingMethod,
        "tilingShape":tilingShape,
        "tilingParameter":tilingParameter,
        "tilingParameterCoordinateType":tilingParameterCoordinateType,
        "weights":weights,
        "interpolation":interpolation,
        "output":output,
        "core_radius":core_radius,
        "half_light_radius":half_light_radius,
        "half_mass_radius":half_mass_radius,
        "append_data_core":append_data_core,
        "append_data_half_light":append_data_half_light,
        "append_data_gbt":append_data_gbt}

    createBeamMatrix(**paras)

def captureNegetiveNumber():
    for i, arg in enumerate(sys.argv):
          if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

def main():
    captureNegetiveNumber()
    parser = argparse.ArgumentParser()
    parseOptions(parser)

if __name__ == "__main__":
    main()
