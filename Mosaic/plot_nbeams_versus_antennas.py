import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import sys
import matplotlib.gridspec as gridspec


data = pd.read_csv('beams_required_versus_antennas_grid.csv')  

data['nantennas_requested'] = data['nantennas_requested'].astype(str)
data = data.loc[data['band'] == 'UHF']

fig = plt.figure(figsize=(22, 14))
gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1,1]) 
ax1 = plt.subplot(gs[0, 0])  # first row, first column
ax2 = plt.subplot(gs[0, 1])  # first row, second column
ax3 = plt.subplot(gs[1, :])  #
gs.update(wspace=0.35, hspace=0.35)


col_order = ["half_light", "gbt", "core"]
overlap_order = ["0.9", "0.95", "0.99"]
overlap_order_gbt = ["0.9", "0.8", "0.6", "0.5"]
#fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(8, 14))
plot_data_half_light = data.loc[(data['region'] == 'half_light') & ((data['overlap'] == 0.9) | (data['overlap'] == 0.95) | (data['overlap'] == 0.99))]
plot_data_half_light['overlap'] = plot_data_half_light['overlap'].astype(str)

plot_data_core = data.loc[(data['region'] == 'core') & ((data['overlap'] == 0.9) | (data['overlap'] == 0.95) | (data['overlap'] == 0.99))]
plot_data_core['overlap'] = plot_data_core['overlap'].astype(str)

plot_data_gbt = data.loc[(data['region'] == 'gbt') & ((data['overlap'] == 0.9) | (data['overlap'] == 0.8) | (data['overlap'] == 0.6) | (data['overlap'] == 0.5))]
plot_data_gbt['overlap'] = plot_data_gbt['overlap'].astype(str)

sns.barplot(x = 'nantennas_requested', y = 'beams_required', hue = 'overlap', data = plot_data_half_light, ax = ax1)
sns.barplot(x = 'nantennas_requested', y = 'beams_required', hue = 'overlap', data = plot_data_core, ax = ax2, hue_order=overlap_order)
sns.barplot(x = 'nantennas_requested', y = 'beams_required', hue = 'overlap', data = plot_data_gbt, ax = ax3)

# Increase this value to increase space between the bars and their labels
vertical_space = 5

for p in ax1.patches:
    ax1.text(p.get_x() + p.get_width() / 2., p.get_height() + vertical_space, f'{int(p.get_height())}', ha='center', fontsize=18)

vertical_space = 0.2

for p in ax2.patches:
    ax2.text(p.get_x() + p.get_width() / 2., p.get_height() + vertical_space, f'{int(p.get_height())}', ha='center', fontsize=18)

vertical_space = 20

for p in ax3.patches:
    ax3.text(p.get_x() + p.get_width() / 2., p.get_height() + vertical_space, f'{int(p.get_height())}', ha='center', fontsize=20)

ax1.set_title('Region: Max.(Half-light/Half-Mass)', fontsize=26, fontweight='bold', y=1.05)
ax2.set_title('Region: Core', fontsize=26, fontweight='bold', y=1.05)
ax3.set_title('Region: GBT L-BAND Primary Beam', fontsize=26, fontweight='bold')

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
plt.suptitle('COMPACT: BEAMS vs NANTENNAS - M30: UHF', fontsize=30, fontweight='bold')
ax1.legend(fontsize=16, title='Overlap fraction', title_fontsize=18, loc=[0.001, 0.7])
ax2.legend(fontsize=16, title='Overlap fraction', title_fontsize=18)
ax3.legend(fontsize=22, title='Overlap fraction', title_fontsize=18)
sns.despine()
plt.tight_layout()
plt.savefig('Nbeams_versus_Antennas_M30_UHF.png', dpi=300)
#plt.show()
print(plot_data_gbt.loc[(plot_data_gbt['nantennas_requested'] == '64') & (plot_data_gbt['overlap'] == '0.99')]['beams_required'])