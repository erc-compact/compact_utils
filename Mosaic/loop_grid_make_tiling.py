import pandas as pd
import sys, subprocess

df = pd.read_csv('COMPACT_GC_CATALOG_FINAL.csv', skiprows=31, na_values='X')
df['Half_light_radius_Baumgardt_arcmin'] = df['Half_light_radius_Baumgardt_arcmin'].astype(float)
df['Half_light_radius_harris_arcmin'] = df['Half_light_radius_harris_arcmin'].astype(float)
central_freq = [816e6, 1284e6, 2406.25e6]
antenna_list = ['subset_40_antennas.txt,subset_52_antennas.txt,subset_56_antennas.txt,subset_60_antennas.txt,full_array.txt']
tiling_shape = 'hexagon'
antenna_file = 'antenna.csv'
overlap_range = ['0.5,0.6,0.8,0.9,0.95,0.99']
tiling_method = 'variable_size'
for index, row in df.iterrows():
    source_name = row['GC']
    boresight = row['RA_hms'] + ' ' + row['DEC_dms']
    mosaic_utc = row['Mosaic_utc']
    max_core_radius_arcmin = max(row['Core_radius_Baumgardt_arcmin'], row['Core_radius_harris_arcmin'])
    max_core_radius_deg = max_core_radius_arcmin / 60
    max_half_light_radius_arcmin = max(row['Half_light_radius_Baumgardt_arcmin'], row['Half_light_radius_harris_arcmin'])
    max_half_light_radius_deg = max_half_light_radius_arcmin / 60
    max_half_mass_radius_deg = row['Half_mass_radius_Baumgardt_arcmin']/60
 
    for freq in central_freq:
        run_command = f"python grid_maketiling.py --freq {freq} --source {boresight} --datetime {mosaic_utc} --subarray_files {antenna_list[0]} --tiling_method {tiling_method} --tiling_shape {tiling_shape} --ants {antenna_file} --core {max_core_radius_deg} --half_light {max_half_light_radius_deg} --half_mass {max_half_mass_radius_deg} --overlap_range {overlap_range[0]} --pointing_name {source_name} --plot_nbeams_antennas"
        #subprocess.check_output(run_command, shell=True)
        process = subprocess.Popen(run_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        # Print output in real-time
        for line in iter(process.stdout.readline, b''):
            print(line.decode().strip())

        process.stdout.close()
        process.wait()  # Wait for the process to finish

        