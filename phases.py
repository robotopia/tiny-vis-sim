import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ITRS, CartesianRepresentation
from astropy.constants import c
import matplotlib.pyplot as plt
import argparse

def main():

    parser = argparse.ArgumentParser(description="Mini simulation for SKA job interview. This code places calibrator source at 10h00m0s -26d0m0s.")
    parser.add_argument("--phase_centre", help="Change phase centre for fringe stopping (in format '00h00m00s 00d00m00s')")
    parser.add_argument("--stop_tracking_after_x_mins", type=float, help="Stop station beam tracking after the specified number of minutes")
    args = parser.parse_args()

    f = np.linspace(149.5, 168, 400) * u.MHz
    λ = c/f
    #t = Time("2025-03-21 12:50:00", scale='utc')
    t = Time(np.linspace(0, 139, 500)*60 + 1426596618, scale='utc', format='gps')

    #F, Tbase = np.meshgrid(f, (t - t[0]).to('min'))

    # SKA Low's approx location
    loc = EarthLocation.of_site('MWA')

    # Put one station at above location, and the other 100 m to the east
    station1 = loc

    baseline_length = 3.5* u.km
    baseline_bearing = 127 * u.deg # degrees east of north
    R_earth = 6371*u.km
    delta_lon = (baseline_length * np.sin(baseline_bearing) / (R_earth * np.cos(loc.lat))) * u.rad
    delta_lat = (baseline_length * np.cos(baseline_bearing) / R_earth) * u.rad

    station2 = EarthLocation(lat=loc.lat + delta_lat, lon=loc.lon + delta_lon, height=loc.height)

    baseline = u.Quantity(station2.to_geocentric()) - u.Quantity(station1.to_geocentric())

    # Choose a source in the sky
    src = SkyCoord('10h00m0s -26d0m0s', frame='icrs')
    src_ecef = src.transform_to(ITRS(obstime=t, location=station1)).cartesian.xyz # Already normalised

    alpha = 0.5
    S_160 = 1.1
    V = S_160*(f[np.newaxis,:]/(160*u.MHz))**alpha * np.exp(-2j*np.pi * np.dot(baseline, src_ecef)[:,np.newaxis] / λ[np.newaxis,:])

    # Try different things that can go wrong

    pointing_ecef = src_ecef
    station_beam_ecef = src_ecef
    measured_baseline = baseline
    separations = np.zeros(t.shape) * u.rad

    if args.phase_centre is not None: # Pointing the telescope (i.e. change phase centre for fringe-stopping and pointing for station beams) somewhere else

        pointing = SkyCoord(args.phase_centre, frame='icrs')
        pointing_ecef = pointing.transform_to(ITRS(obstime=t, location=station1)).cartesian.xyz # Already normalised

        station_beam = SkyCoord(args.phase_centre, frame='icrs')
        station_beam_ecef = station_beam.transform_to(ITRS(obstime=t, location=station1)).cartesian.xyz # Already normalised

        separations = np.arccos(np.einsum('ij,ij->j', station_beam_ecef, src_ecef))

    if args.stop_tracking_after_x_mins is not None:

        stop_tracking_time = args.stop_tracking_after_x_mins * u.min

        t_idx = np.where((t - t[0]).to('min') > stop_tracking_time)[0][0]
        pointing_at_stopping_time = station_beam_ecef[:,t_idx]
        print(f'{pointing_at_stopping_time.shape = }')
        print(f'{src_ecef[:,t_idx:].shape = }')
        print(f'{separations[t_idx:].shape = }')
        separations[t_idx:] = np.arccos(np.dot(pointing_at_stopping_time, src_ecef[:,t_idx:]))

    # Choose a gaussian with FWHM = 3 degrees and apply the primary beam
    FWHM = 3 * u.deg
    σ = FWHM / 2.355
    pb_corr = np.exp(-0.5*(separations[:,np.newaxis]/σ)**2)

    Vhat = pb_corr * V *  np.exp(2j*np.pi * np.dot(measured_baseline, pointing_ecef)[:,np.newaxis] / λ[np.newaxis,:])

    fig = plt.figure(figsize=(9,7))
    #plt.pcolormesh(f.to('MHz').value, (t - t[0]).to('min').value, np.angle(V).value)
    plt.pcolormesh(f.to('MHz').value, (t - t[0]).to('min').value, np.imag(Vhat))
    cbar = plt.colorbar()
    cbar.set_label("Imag. part of vis")
    plt.xlabel("Frequency (MHz)")
    plt.ylabel("Time (min)")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig("phases.png")

if __name__ == '__main__':
    main()
