import argparse
import glob
import matplotlib.pyplot as plt
import numpy as np

def read_file(fname, total_antennas, used_antennas):
    dtype = [("mean", "float32"), 
             ("std", "float32"), 
             ("skew", "float32"), 
             ("kurtosis", "float32")]
    ar = np.fromfile(fname, dtype=dtype, offset=4096)
    ar = ar.reshape(-1, 64, 2, total_antennas)[:, :, :, :used_antennas]
    return ar

def get_files(dirname):
    return sorted(glob.glob(f"{dirname}/*.fpa"), key=lambda x: int(x.split("_")[-2]))

def make_views(ar, threshold=0.08):
    # ar is in TFPA order
    ar_t  = (abs(ar['kurtosis']) > threshold)
    views = {"t": ar_t.mean(axis=(1,2,3)),
             "tf": ar_t.mean(axis=(2,3)),
             "tp": ar_t.mean(axis=(1,3)),
             "ta": ar_t.mean(axis=(1,2)),
             "f": ar_t.mean(axis=(0,2,3)),
             "fp": ar_t.mean(axis=(0,3)),
             "fa": ar_t.mean(axis=(0,2)),
             "p": ar_t.mean(axis=(0,1,3)),
             "pa": ar_t.mean(axis=(0,1)),
             "a": ar_t.mean(axis=(0,1,2))
             }
    return views

def merge_views(v0, v1):
    views = {"t": np.mean((v0["t"], v1["t"]), axis=0),
             "tf": np.concatenate((v0["tf"], v1["tf"]), axis=1),
             "tp": np.mean((v0["tp"], v1["tp"]), axis=0),
             "ta": np.mean((v0["ta"], v1["ta"]), axis=0),
             "f": np.concatenate((v0["f"], v1["f"])),
             "fp": np.concatenate((v0["fp"], v1["fp"]), axis=0),
             "fa": np.concatenate((v0["fa"], v1["fa"]), axis=0),
             "p": np.mean((v0["p"], v1["p"]), axis=0),
             "pa": np.mean((v0["pa"], v1["pa"]), axis=0),
             "a": np.mean((v0["a"], v1["a"]), axis=0)
             }
    return views

def plot_views(views):
    imshow_args = {"aspect": "auto", "interpolation": "nearest"}
    
    ax_ta = plt.subplot2grid((4,4), (3,0))
    ax_ta.imshow(views["ta"].T, **imshow_args)
    ax_ta.set_xlabel("Time")
    ax_ta.set_ylabel("antenna")
    
    ax = plt.subplot2grid((4,4), (0,0), sharex=ax_ta)
    ax.plot(views["t"])
    ax.xaxis.set_visible(False)
    ax.set_ylabel("masked fraction")
    
    ax = plt.subplot2grid((4,4), (1,0), sharex=ax_ta)
    ax.imshow(views["tf"].T, **imshow_args)
    ax.xaxis.set_visible(False)
    ax.set_ylabel("channel")
    
    ax = plt.subplot2grid((4,4), (2,0), sharex=ax_ta)
    ax.imshow(views["tp"].T, **imshow_args)
    ax.xaxis.set_visible(False)
    ax.set_ylabel("polarisation")
    
    ax = plt.subplot2grid((4,4), (1,1))
    ax.plot(views["f"])
    ax.xaxis.set_visible(False)
    
    ax = plt.subplot2grid((4,4), (2,1))
    ax.imshow(views["fp"].T, **imshow_args)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    
    ax = plt.subplot2grid((4,4), (3,1))
    ax.imshow(views["fa"].T, **imshow_args)
    ax.set_xlabel("channel")
    ax.yaxis.set_visible(False)
    
    ax = plt.subplot2grid((4,4), (2,2))
    ax.plot(views["p"])
    ax.xaxis.set_visible(False)
    
    ax = plt.subplot2grid((4,4), (3,2))
    ax.imshow(views["pa"].T, **imshow_args)
    ax.yaxis.set_visible(False)
    ax.set_xlabel("polarisation")
    
    ax = plt.subplot2grid((4,4), (3,3))
    ax.plot(views["a"])
    ax.set_xlabel("antenna")
    
def run(files, total_antennas, used_antennas):
    print(files)
    views = make_views(read_file(files[0], total_antennas, used_antennas))
    for ii, fname in enumerate(files[1:]):
        print(f"Processing file {ii} of {len(files) - 1}")
        views = merge_views(views, make_views(read_file(fname, total_antennas, used_antennas)))
    plt.figure()
    plot_views(views)
    plt.show()
    return views
    
    
def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Calculate remaining antennas.")
    
    # Add arguments for total and used antennas
    parser.add_argument(
        '--total-antennas', 
        type=int, 
        default=64, 
        help='The total number of antennas available.'
    )
    parser.add_argument(
        '--used-antennas', 
        type=int, 
        required=True, 
        help='The number of valid antennas use.'
    )
    parser.add_argument(
        '--dir', 
        type=str, 
        required=True, 
        help='The directory holding the fpa files'
    )
    
    # Parse the arguments
    args = parser.parse_args()
    run(get_files(args.dir), args.total_antennas, args.used_antennas)

if __name__ == "__main__":
    main()

