import pandas as pd
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy import units as u
import sys, math
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist, squareform
import seaborn as sns
import matplotlib.pyplot as plt
# read data from file
df = pd.read_csv('meerkat_antenna_positions.csv', delim_whitespace=True, skiprows=1, names=['Antenna', 'East', 'North', 'Up'], dtype='str')

# Clean antenna names
df['Antenna'] = df['Antenna'].apply(lambda x: x.strip("m"))
df['North'] = df['North'].astype(float)
df['East'] = df['East'].astype(float)
df['Up'] = df['Up'].astype(float)


df.set_index('Antenna', inplace=True)

# Step 2: Create a copy of the original DataFrame
df_original = df.copy()

# Step 3: Prepare a list to store data for the plot
plot_data = []

# Step 4: Calculate the distance matrix and maximum distance
#dist_matrix = squareform(pdist(df.values))
dist_matrix = distance_matrix(df.values, df.values)
# Add max.baseline for full array.
plot_data.append(('None', dist_matrix.max(), len(df)))
full_baseline = dist_matrix.max()

# Step 5: Remove the antenna that is farthest from the center and update the plot data
while len(df) > 1:
    df['dist_from_center'] = np.sqrt(df['East']**2 + df['North']**2 + df['Up']**2)
 
    farthest_antenna = df['dist_from_center'].idxmax()
    
    df.drop(farthest_antenna, inplace=True)
    new_dist_matrix = squareform(pdist(df.values))
    
    # Add data for the plot
    plot_data.append((farthest_antenna, new_dist_matrix.max(), len(df)))
    
    # Update the distance matrix
    dist_matrix = new_dist_matrix

# Step 6: Create a DataFrame for the plot data and plot it
plot_df = pd.DataFrame(plot_data, columns=['DroppedAntenna', 'MaxBaseline', 'RetainedAntennas'])


def max_baseline_to_percentage_drop_off(x):
    return (x/full_baseline) * 100


def percentage_drop_off_to_max_baseline(x):
    return (x/100) * full_baseline


#print(plot_df.loc[plot_df['MaxBaseline'] <=6500])

fig, ax1 = plt.subplots(figsize=(10,6))

color = 'tab:red'
ax1.set_xlabel('Dropped Antenna')
ax1.set_ylabel('Maximum Baseline (m)', color=color)
ax1.scatter(plot_df['DroppedAntenna'], plot_df['MaxBaseline'], color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.tick_params(axis='x', labelcolor=color, rotation=90)
ax1.axhline(6031.596988 , color='black', linewidth=0.2, linestyle='solid')
ax1.axhline(4422.287242 , color='black', linewidth=0.2, linestyle='solid')
ax1.axhline(3918.131479 , color='black', linewidth=0.2, linestyle='solid')
ax1.axhline(2343.512346 , color='black', linewidth=0.2, linestyle='solid')
ax1.axhline(1146.109969 , color='black', linewidth=0.2, linestyle='solid')

ax1.axvline(5, color='black', linewidth=0.5, linestyle='--')
ax1.axvline(7, color='black', linewidth=0.5, linestyle='--')
ax1.axvline(9, color='black', linewidth=0.5, linestyle='--')
ax1.axvline(13, color='black', linewidth=0.5, linestyle='--')
ax1.axvline(21, color='black', linewidth=0.5, linestyle='--')


ax2 = ax1.twiny()  # instantiate a second axes that shares the same y-axis

color = 'tab:blue'
ax2.set_xlabel('Number of Retained Antennas', color=color)  # we already handled the x-label with ax1
ax2.plot(plot_df['RetainedAntennas'].astype(str), np.zeros_like(plot_df['RetainedAntennas']), color=color)  # we already handled the y-label with ax1
ax2.tick_params(axis='x', labelcolor=color, rotation=90)
ax2.set_xticks(['64', '59', '58', '57', '55', '51', '43', '41', '40', '32', '16', '8'])  # set x-ticks manually
plt.ylim(0.1, 8000)
plt.suptitle('MeerKAT Max. Baseline vs No. of Antennas', fontsize=16)
#sns.despine()

#ax3 = ax1.twinx()  # instantiate a second y-axis

color = 'tab:blue'
secay = ax1.secondary_yaxis('right', functions=(max_baseline_to_percentage_drop_off, percentage_drop_off_to_max_baseline), color=color)
secay.set_ylabel('Percentage Max. Baseline Drop-off', color=color)


#percentage_dropoff = 100 * (1 - plot_df['MaxBaseline'] / full_baseline)
# ax3.plot(plot_df['DroppedAntenna'], percentage_dropoff, color=color)
# ax3.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('meerkat_baseline_drop_off.png', dpi=600)
#plt.show()



