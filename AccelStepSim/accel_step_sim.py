import numpy as np
import matplotlib.pyplot as plt
from sigpyproc.timeseries import TimeSeries
from sigpyproc.header import Header

def make_signal(profile, period, nsamples, tsamp, amp=1):
    """Make a noiseless timeseries with a repeating profile

    Args:
        profile (iterable): An array of arbitrary size containing the desired pulse shape
        period (float): The repitition period of the pulse in seconds
        nsamples (int): The total number of samples to produce
        tsamp (float): The sampling interval in seconds
        amp (int, optional): A scale factor to apply to the signal. Defaults to 1.

    Returns:
        np.ndarray: An array containing the signal
    """
    profile = np.asarray(profile) * amp
    prof_phase = np.arange(profile.size) / profile.size
    phases = np.modf((np.arange(nsamples) * tsamp) / period)[0]
    return np.interp(phases, prof_phase, profile)

def plot_spec(t, acc):
    """Plot the spectrum and power distribution per fold

    Args:
        t (_type_): A sigpyproc timeseries
        acc (_type_): An acceleration to resample to

    Returns:
        tuple(axis, axis): The plotting axes
    """
    s = t.resample(acc).rfft().form_spec(interpolate=True)
    freq = s.bin2freq(1) * np.arange(s.data.size)
    s.data[:] -= s.data.mean()
    s.data[:] /= np.sqrt(((s.data.std()**2)/4))
    s.data[:] += 2
    s.data[:1000] = 0
    peaks = [s.data.max()]
    folds = s.harmonic_fold(5)
    ax = plt.subplot(121)
    plt.plot(freq, s.data, label="f1")
    for ii, fold in enumerate(folds):
        n = 2**(ii+1)
        normed = (fold.data - n)/ np.sqrt(n)
        peaks.append(normed.max())
        plt.plot(freq/n, normed, label=f"f{n}")
    plt.legend()
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power")
    ax1 = plt.subplot(122, sharey=ax)
    plt.plot(peaks)
    plt.xlabel("Fold number")
    return ax, ax1

def get_peak(t, acc):
    """Get the peak power for the folds
    """
    s = t.resample(acc).rfft().form_spec(interpolate=True)
    s.data[:] -= s.data.mean()
    s.data[:] /= np.sqrt(((s.data.std()**2)/4))
    s.data[:] += 2
    s.data[:1000] = 0
    peaks = [s.data.max()]
    folds = s.harmonic_fold(5)
    for ii, fold in enumerate(folds):
        n = 2**(ii+1)
        normed = (fold.data - n)/ np.sqrt(n)
        peaks.append(normed.max())
    return peaks

def make_sensitivity_map(t, amax, steps, daccel=[]):
    """Make a map of sensitivity loss as a function of harmonic and acceleration

    Args:
        t (_type_): A sigpyproc timeseries
        amax (_type_): The maximum acceleration
        steps (_type_): The step size
        daccel (list, optional): a list of thresholds to overplot with form (thresh, name, colour). Defaults to [].
    """
    ar = []
    acc = np.linspace(0, amax, steps)
    for a in acc:
        print(a)
        ar.append(get_peak(t, a))
    ar = np.array(ar)
    ar /= ar[0].max()
    plt.imshow(ar, aspect='auto', extent=[-0.5,5.5, amax, 0])
    plt.xlabel("Fold number")
    plt.ylabel("Acceleration offset (m/s/s)")
    c = plt.colorbar()
    c.set_label("Fractional sensitivity")
    for dacc in daccel:
        plt.hlines(dacc[0], -0.5, 5.5, label=dacc[1], color=dacc[2])
    plt.legend(loc="lower left")
    
if __name__ == "__main__":
    # Below is an example of how this can be used
    nsamples = 2**23
    tsamp = 76e-6
    c = 299792458.0
    tobs = tsamp * nsamples
    period = 0.023 # seconds
    
    # make a 1% duty cycle square wave
    delta_1percent = np.zeros(100)
    delta_1percent[50] = 1
    
    ralph = 64 * c * tsamp / tobs**2
    z1_nyquist = 2 * c / (tobs**2 * 1 / (2 * tsamp))
    peasoup = lambda tol: 2.0 * tsamp * 24.0 * c/tobs**2 * np.sqrt((tol*tol)-1.0)
    
    # some noise
    base = np.random.normal(0, 1, nsamples)
    
    # nonsense header needed for sigpyproc
    header = Header(filename="dummy.tim", data_type=2, nchans=1, foff=1.0, fch1=1500, nbits=32, tsamp=tsamp, tstart=0, nsamples=nsamples)
    
    
    t = TimeSeries(base + make_signal(delta_1percent, period, nsamples, tsamp, amp=0.4), header)
    plt.figure()
    make_sensitivity_map(t, 10, 100, daccel=[
        (ralph, "ralph", "m"), 
        (peasoup(1.1), "peasoup(1.1)", "k"), 
        (z1_nyquist, "z=1 @ nyquist", "b")])
    plt.show()