import numpy as np
import sys
from pathlib import Path
import argparse
import yaml, os
import argparse

def read_yaml(file_path):
    with open(file_path, 'r') as file:
        data = yaml.safe_load(file)
    return data

def parse_args():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-r', '--skycleaver_out_root', type=str, required=True , help='Root directory of skycleaver output')
    argparser.add_argument('-o', '--output_csv', type=str, required=True, help='Output csv file')
    argparser.add_argument('-v', '--skycleaver_naming_version', type=int, default=1, help='Version of the skycleaver naming convention')
    argparser.add_argument('-s', '--source_name', type=str, default='J0514-4002A', help='Name of the source')
    argparser.add_argument('-u', '--utc_start_dir', type=str, default='2024-05-19-15:50:23', help='UTC start time of observation')
    argparser.add_argument('-c', '--coherent_dms', type=str, default='', help='Coherent DMs, default ALL')
    argparser.add_argument('-i', '--incoherent_dms', type=str, default='', help='Incoherent DMs, default ALL')
    argparser.add_argument('--stokes', type=str, default='', help='Stokes parameter, default ALL')
    argparser.add_argument('-t', '--targets_file', type=str, help='Targets file')
    argparser.add_argument('-y', '--yaml_file', type=str, help='YAML file')
    argparser.add_argument('-b', '--beam_sets', type=str, help='only choose beamsets with this prefix', default='')
    argparser.add_argument("--beamformed_dir", type=str, default='01_BEAMFORMED', help='Directory containing beamformed files')
    argparser.add_argument("--filterbank_dir", type=str, default='02_FILTERBANKS', help='Directory containing filterbank files')
    return argparser.parse_args()

def main():

    args = parse_args()
    source_name = args.source_name
    utc_start_dir = args.utc_start_dir
    root = Path(args.skycleaver_out_root) 
    utc_start = utc_start_dir if '_' not in utc_start_dir else utc_start_dir.split('_')[0]
    print(f'Root directory: {root}')
    beamformed_root = root / args.beamformed_dir / source_name / utc_start_dir

    if args.targets_file == None:
        print(f"Targets file not provided, will take it from {beamformed_root} if there is only one file")
        targets_files = list(beamformed_root.glob('**/*.targets'))
        if len(targets_files) == 1:
            targets_file = targets_files[0]
        elif len(targets_files) > 1:
            print(f"Multiple targets files found in {beamformed_root}, explicitly provide the targets file using --targets_file")
            sys.exit()
        else:
            print(f"No targets file found in {beamformed_root}, explicitly provide the targets file using --targets_file")
            sys.exit()
    else:
        targets_file = beamformed_root / args.targets_file    
    if args.yaml_file == None:
        print(f"YAML file not provided, will take it from {beamformed_root} if there is only one file")
        yaml_files = list(beamformed_root.glob('**/*.yaml'))
        if len(yaml_files) == 1:
            yaml_file = yaml_files[0]
        elif len(yaml_files) > 1:
            print(f"Multiple yaml files found in {beamformed_root}, explicitly provide the yaml file using --yaml_file")
            sys.exit()
        else:
            print(f"No yaml file found in {beamformed_root}, explicitly provide the yaml file using --yaml_file")
            sys.exit()
    else:
        yaml_file = beamformed_root / args.yaml_file



    print(f'Targets file: {targets_file}')
    print(f'YAML file: {yaml_file}')



    if not root.exists():
        print(f'Root directory does not exist: {root}')
        sys.exit()
    if not targets_file.exists():
        print(f'Targets file does not exist: {targets_file}')
        sys.exit()
    if not yaml_file.exists():
        print(f'YAML file does not exist: {yaml_file}')
        sys.exit()
    
    #load file containing name,ra,dec,x,y,angle,beam_set_id,overlap,nantennas,type as numpy array
    targets = np.loadtxt(targets_file, delimiter=',', dtype={'names': ('beam_name', 'ra', 'dec', 'x', 'y', 'angle', 'beam_set_id', 'overlap', 'nantennas', 'type'), 
                                                             'formats': ('U50', 'U50', 'U50', 'f8', 'f8', 'f8', 'i4', 'f8', 'i4', 'U50')}, skiprows=1)
    
    data = read_yaml(yaml_file)
    
    beamsets = data.get('beam_sets', [])

    filterbank_dir = root /   args.filterbank_dir /source_name / utc_start_dir
    print(f'Filterbank directory: {filterbank_dir}')
    

    fil_files = list(filterbank_dir.glob(f'*/**/*.fil'))
   


    print(f'Found {len(fil_files)} fil files')
    if len(fil_files) == 0:
        print(f'No fil files found in {root}')
        sys.exit()



    lines = []
    out_header="target,beam_name,utc_start,filenames,filstr,coherent_dm,subband_dm,hardware_name,beam_type,bf_reference_freq,bf_ra,bf_dec,bf_utc,bf_sub_array,bf_tiling_method,bf_tiling_shape,bf_overlap,beam_shape_x,beam_shape_y,beam_shape_angle"
    lines.append(out_header)


    cdms = args.coherent_dms.split(',') if args.coherent_dms != '' else [f.name.split('_')[3] for f in fil_files]
    stokes = args.stokes.split(',') if args.stokes != '' else [f.name.split('_')[-2] for f in fil_files]
    cdms =np.unique(cdms)
    stokes = np.unique(stokes)
    print(f'CDMs: {cdms}')
    print(f'Stokes: {stokes}')
    for target in targets: 
        beam_name = target['beam_name']
        if args.beam_sets != '' and not beam_name.startswith(args.beam_sets):
            continue
        print(f'Processing beam: {beam_name}')
        antennas = None
        reference_frequency = None
        for beamset in beamsets:
            if beamset['name'] == beam_name.split('_')[0]:
                antennas = " ".join(beamset['antenna_set']).replace('m','')
                tilings = beamset.get('tilings', [])
                if len(tilings) > 0:
                    reference_frequency = tilings[0].get('reference_frequency')
                break
        for stok in stokes:        
            for cdm in cdms:
                cdm_filterbanks = [f for f in fil_files if f"cdm_{cdm}_" in f.name]
                idms = args.incoherent_dms.split(',') if args.incoherent_dms != '' else [f.name.split('_')[5] for f in cdm_filterbanks]   
                idms = list(set(idms))      
                print(f'IDMs: {idms}')
                for idm in idms:
                    shortlisted = [f for f in cdm_filterbanks if f"cdm_{cdm}_idm_{idm}_{beam_name}" in f.name]
                    shortlisted = sorted(shortlisted, key=lambda x: x.name.replace(".fil", "").split('_')[-1])
                    fil_names = ' '.join([f.absolute().as_posix() for f in shortlisted])
                    filstr = '_'.join(shortlisted[0].name.split('_')[1:-1])
                    lines.append(f"{source_name},{beam_name},{utc_start},{fil_names},{filstr},{cdm},{idm},'contra','STOKES_{stok}',816000000.0,{target['ra']},{target['dec']},{utc_start},{antennas},variable_size,circle,{target['overlap']},{target['x']},{target['y']},{target['angle']}")
            
    with open(args.output_csv, 'w') as f:
        for line in lines:
            f.write(f"{line}\n")











    





if __name__ == '__main__':
    main()
   