import numpy as np
import pandas as pd
from astro_coordinates_converter import hours_to_decimal, decimal_to_hours, equatorial_to_galactic
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u
from astroplan import Observer, FixedTarget
from pytz import timezone
from astropy.time import Time
from astropy.table import Table
import sys
from core_half_mass_pc_to_arcmin import calculate_theta

distance_pc_sun_earth = 4.84e-6
df = pd.read_csv('COMPACT_GC_CLUSTERS.txt', skiprows=13, delim_whitespace=True)
RA_hour_angle = df['RA(HA)'].values
DEC_degrees = df['DEC(DEG)'].values
# Convert coordinates to standard units
RA_deg, DEC_deg = hours_to_decimal(RA_hour_angle, DEC_degrees)
coordinates_hms_dms = decimal_to_hours(RA_deg, DEC_deg)
RA_hms, DEC_dms = zip(*[val.split() for val in coordinates_hms_dms])
l, b = equatorial_to_galactic(RA_deg, DEC_deg, unit='deg')

df['RA_deg'] = RA_deg
df['DEC_deg'] = DEC_deg
df['RA_hms'] = RA_hms
df['DEC_dms'] = DEC_dms
df['Gl'] = l
df['Gb'] = b

core_radius_pc = df['CORE_RADIUS_Baumgardt'].values
half_light_radius_pc = df['HALF_LIGHT_RADIUS_Baumgardt'].values
half_mass_radius_pc = df['HALF_MASS_RADIUS_Baumgardt'].values
distance_pc = (df['R_Sun_Baumgardt'].values * 1000) + distance_pc_sun_earth

core_radius_arcmin = calculate_theta(core_radius_pc, distance_pc) 
half_light_radius_arcmin = calculate_theta(half_light_radius_pc, distance_pc)
half_mass_radius_arcmin = calculate_theta(half_mass_radius_pc, distance_pc)

df['Core_radius_Baumgardt_arcmin'] = core_radius_arcmin
df['Half_light_radius_Baumgardt_arcmin'] = half_light_radius_arcmin
df['Half_mass_radius_Baumgardt_arcmin'] = half_mass_radius_arcmin

# Calculate thetas for the Baumgardt catalog.



mkt_location = {
        "elevation": 1086.599484882955,
        "name": "MeerKAT",
        "longitude_unit": "degree",
        "latitude_unit": "degree",
        "latitude": -30.711055553291878,
        "elevation_unit": "meter",
        "longitude": 21.443888889697842,
        "source": "MEERKAT, used in timing mode.\n\n    The origin of this data is unknown but as of 2021 June 8 it agrees exactly with\n    the values used by TEMPO and TEMPO2.\nvia PINT",
	"timezone": "Africa/Johannesburg",
	"aliases": ["MeerKAT"]}

# MeerKAT Location
location = EarthLocation.from_geodetic(mkt_location['longitude'] * u.deg, \
    mkt_location['latitude'] * u.deg, mkt_location['elevation'] * u.m)


meerkat = Observer(name='MeerKAT',
               location=location,
               timezone=timezone('Africa/Johannesburg'),
               description="MeerKAT telescope in SA")


fixed_time = Time("2023-05-16 00:00")

targets_df = df[['GC', 'RA_hms', 'DEC_dms']]
target_table = Table.from_pandas(targets_df)
rise_set_error_margin = 10*u.minute
targets = [FixedTarget(coord=SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg)), name=name)
            for name, ra, dec in target_table]

rise_time_utc_20_deg_list = []
set_time_utc_20_deg_list = []
rise_time_utc_50_deg_list = []
set_time_utc_50_deg_list = []
rise_time_sidereal_20_deg_list = []
set_time_sidereal_20_deg_list = []
rise_time_sidereal_50_deg_list = []
set_time_sidereal_50_deg_list = []
mosaic_utc_list = []

for i in range(len(targets)):
    rise_time_utc_20_deg = meerkat.target_rise_time(fixed_time, targets[i], n_grid_points=10, horizon=20 * u.deg).iso
    set_time_utc_20_deg = meerkat.target_set_time(fixed_time, targets[i], n_grid_points=10, horizon=20 * u.deg).iso

    rise_time_utc_50_deg = meerkat.target_rise_time(fixed_time, targets[i], n_grid_points=10, horizon=50 * u.deg).iso
    set_time_utc_50_deg = meerkat.target_set_time(fixed_time, targets[i], n_grid_points=10, horizon=50 * u.deg).iso

    # Rise set time Sidereal at 20 deg.
    rise_time_sidereal_20_deg = meerkat.local_sidereal_time(rise_time_utc_20_deg).to_string(sep=':')
    set_time_sidereal_20_deg = meerkat.local_sidereal_time(set_time_utc_20_deg).to_string(sep=':')
       
    # Rise set time Sidereal at 50 deg.
    rise_time_sidereal_50_deg = meerkat.local_sidereal_time(rise_time_utc_50_deg).to_string(sep=':')
    set_time_sidereal_50_deg = meerkat.local_sidereal_time(set_time_utc_50_deg).to_string(sep=':')
    mosaic_utc = Time(rise_time_utc_20_deg) + 20 * u.minute
    rise_time_utc_20_deg_list.append(rise_time_utc_20_deg)
    set_time_utc_20_deg_list.append(set_time_utc_20_deg)
    rise_time_utc_50_deg_list.append(rise_time_utc_50_deg)
    set_time_utc_50_deg_list.append(set_time_utc_50_deg)
    rise_time_sidereal_20_deg_list.append(rise_time_sidereal_20_deg)
    set_time_sidereal_20_deg_list.append(set_time_sidereal_20_deg)
    rise_time_sidereal_50_deg_list.append(rise_time_sidereal_50_deg)
    set_time_sidereal_50_deg_list.append(set_time_sidereal_50_deg)
    mosaic_utc_list.append(mosaic_utc.iso)
    


df['Rise_time_utc_20_deg'] = rise_time_utc_20_deg_list
df['Set_time_utc_20_deg'] = set_time_utc_20_deg_list
df['Rise_time_utc_50_deg'] = rise_time_utc_50_deg_list
df['Set_time_utc_50_deg'] = set_time_utc_50_deg_list
df['Rise_time_sidereal_20_deg'] = rise_time_sidereal_20_deg_list
df['Set_time_sidereal_20_deg'] = set_time_sidereal_20_deg_list
df['Rise_time_sidereal_50_deg'] = rise_time_sidereal_50_deg_list
df['Set_time_sidereal_50_deg'] = set_time_sidereal_50_deg_list
df['Mosaic_utc'] = mosaic_utc_list
del df['RA(HA)']
del df['DEC(DEG)']

# Rearrange Pandas columns to make it easier to read.
df = df[['GC', 'RA_hms', 'DEC_dms', 'RA_deg', 'DEC_deg', 'Gl', 'Gb', 'Gamma', 'gamma', \
        'Core_radius_harris_arcmin', 'Core_radius_Baumgardt_arcmin', 'Half_light_radius_harris_arcmin', \
        'Half_light_radius_Baumgardt_arcmin', 'Half_mass_radius_Baumgardt_arcmin', 'CORE_RADIUS_Baumgardt', \
         'HALF_LIGHT_RADIUS_Baumgardt', 'HALF_MASS_RADIUS_Baumgardt', 'R_Sun_harris', 'R_Sun_Baumgardt', \
         'Rise_time_utc_20_deg', 'Set_time_utc_20_deg', 'Rise_time_utc_50_deg', 'Set_time_utc_50_deg', \
         'Rise_time_sidereal_20_deg', 'Set_time_sidereal_20_deg', 'Rise_time_sidereal_50_deg', \
         'Set_time_sidereal_50_deg', 'Mosaic_utc', 'CORE_COLLAPSED']]

print(df.columns)
# Save the dataframe to a csv file.
df.to_csv('COMPACT_GC_CATALOG_FINAL.csv', index=False)