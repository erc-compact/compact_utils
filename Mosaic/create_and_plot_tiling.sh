#!/bin/bash
# Compact overlap 0.9 tiling
python maketiling.py --freq 1.284e9 --source 21:40:22.12 -23:10:47.5 --datetime 2023.05.15 23:43:00.0 --subarray 000,001,002,003,004,005,006,007,008,009,010,011,012,013,014,015,016,017,018,019,020,021,022,023,024,025,026,027,028,029,030,031,032,033,034,035,036,037,038,039,040,041,042,043,044,045,046,047,048,049,050,051,052,053,054,055,056,057,058,059,060,061,062,063 --verbose --tiling_method variable_size --tiling_shape hexagon --tiling_plot tiling3.png --ants antenna.csv --tiling_coordinate coordinate.csv --overlay_source extra_source.txt --core 0.001 --half_light 0.01766666666 --half_mass 0.02891881893 --overlap 0.9 --beamnum 380 --tiling_position_file compact_m30_overlap_90.txt --pointing_name M30
#Compact overlap 0.95 tiling
python maketiling.py --freq 1.284e9 --source 21:40:22.12 -23:10:47.5 --datetime 2023.05.15 23:43:00.0 --subarray 000,001,002,003,004,005,006,007,008,009,010,011,012,013,014,015,016,017,018,019,020,021,022,023,024,025,026,027,028,029,030,031,032,033,034,035,036,037,038,039,040,041,042,043,044,045,046,047,048,049,050,051,052,053,054,055,056,057,058,059,060,061,062,063 --verbose --tiling_method variable_size --tiling_shape hexagon --tiling_plot tiling3.png --ants antenna.csv --tiling_coordinate coordinate.csv --overlay_source extra_source.txt --core 0.001 --half_light 0.01766666666 --half_mass 0.02891881893 --overlap 0.95 --beamnum 1400 --tiling_position_file compact_m30_overlap_95.txt --pointing_name M30 

python plot_multiple_tiling.py --tiling_files compact_m30_overlap_90.txt,compact_m30_overlap_95.txt  --tiling_plot M30_LBAND_COMPACT_OVERLAP_90_95 --core 0.001 --half_light 0.01766666666 --half_mass 0.02891881893 --extra_source extra_source.txt



# Plot compact and trapum tiling together
python plot_multiple_tiling.py --tiling_files TRAPUM_M30_LBAND_DEC_29_overlap_0_7.txt,compact_m30_overlap_90.txt  --tiling_plot M30B_TRAPUM_COMPACT_BEAM_TILING --core 0.001 --half_light 0.01766666666 --half_mass 0.02891881893 --extra_source extra_source.txt

