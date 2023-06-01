import pandas as pd
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy import units as u
import sys, math
# read data from file
df = pd.read_csv('meerkat_antenna_positions.csv', delim_whitespace=True, skiprows=1, names=['Antenna', 'East', 'North', 'Up'], dtype='str')

# Clean antenna names
df['Antenna'] = df['Antenna'].apply(lambda x: x.strip("m"))
df['North'] = df['North'].astype(float)
df['East'] = df['East'].astype(float)
df['Up'] = df['Up'].astype(float)
# mkt_location = {
#         "elevation": 1086.599484882955,
#         #"elevation": 0,
#         "name": "MeerKAT",
#         "longitude_unit": "degree",
#         "latitude_unit": "degree",
#         "latitude": -30.711055553291878,
#         "elevation_unit": "meter",
#         "longitude": 21.443888889697842,
#         "source": "MEERKAT, used in timing mode.\n\n    The origin of this data is unknown but as of 2021 June 8 it agrees exactly with\n    the values used by TEMPO and TEMPO2.\nvia PINT",
# 	"timezone": "Africa/Johannesburg",
# 	"aliases": ["MeerKAT"]}

# MeerKAT phase centre
ref_location = EarthLocation.from_geodetic(mkt_location['longitude'] * u.deg, \
    mkt_location['latitude'] * u.deg, mkt_location['elevation'] * u.m)




# def convert_to_abs_location(row):
#     # Convert offsets to angular measures
    
#     # Compute the absolute latitude, longitude, and height
#     abs_lat =  ref_location.y  + row['North'] * u.m
#     abs_lon = ref_location.x  + row['East'] * u.m
#     abs_height = ref_location.z + row['Up'] * u.m
#     return EarthLocation(abs_lat, abs_lon, abs_height)


# Convert ENU to EarthLocation for each antenna
df['location'] = df.apply(convert_to_abs_location, axis=1)


meerkat_phase_center = [0 * u.m, 0 * u.m, 0 * u.m]

distance = []
for index, row in df.iterrows():
    distance_from_center = np.sqrt((meerkat_phase_center[0] - row['East'] * u.m)**2 + (meerkat_phase_center[1] - row['North'] * u.m)**2 + (meerkat_phase_center[2] - row['Up'] * u.m)**2)
    distance.append(distance_from_center)

# Calculate distances to the phase center
df['distance_from_center_meter'] = distance

# Sort the DataFrame by distance
df.sort_values(by='distance_from_center_meter', inplace=True)
df['distance_from_center_meter'] = df.apply(lambda x: x['distance_from_center_meter'].value, axis=1)
# Save full array
df['Antenna'].to_csv('full_array.txt', index=False, header=False)

#df1 = df.loc[df['distance_from_center_meter'] < 1000] 

# Save the subsets
# for antennas in [32, 40, 16, 8, 4]:
#     subset = df['Antenna'][:antennas]
#     subset.to_csv(f'subset_{antennas}.txt', index=False, header=False)










#df['distance'] = np.sqrt(((df['coord'] - ref_coord)**2).sum(axis=1))
#print(df)
# Sort the DataFrame by distance
df.sort_values(by='distance', inplace=True)

# Save full array
df['Antenna'].to_csv('full_array.txt', index=False, header=False)

# Save the subsets
for antennas in [32, 40, 16, 8, 4]:
    subset = df['Antenna'][:antennas]
    subset.to_csv(f'subset_{antennas}.txt', index=False, header=False)
