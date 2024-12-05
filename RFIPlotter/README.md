# RFI Plotter

This script provides a quick look at .fpa files output by skyweaver. Example usage:

python rfi_file_plotter.py --dir /b/PROCESSING_NGC1851_15MIN/J0514-4002A/2024-05-19-15:50:23/0/*/ --total-antennas 64 --used-antennas 62

The resultant plot is the result of thresholding the kurtosis statistic with a value of 0.08 and then integrating over the different combinations of dimensions. In effect this gives the masked fraction for the given threshold.
