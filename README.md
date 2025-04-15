# Toy visibility simulation

## Requirements

- `numpy`
- `astropy`
- `matplotlib`
- `argparse`

## Usage

```
usage: phases.py [-h] [--phase_centre PHASE_CENTRE] [--stop_tracking_after_x_mins STOP_TRACKING_AFTER_X_MINS]

Mini simulation for SKA job interview. This code places calibrator source at 10h00m0s -26d0m0s.

options:
  -h, --help            show this help message and exit
  --phase_centre PHASE_CENTRE
                        Change phase centre for fringe stopping (in format '00h00m00s 00d00m00s')
  --stop_tracking_after_x_mins STOP_TRACKING_AFTER_X_MINS
                        Stop station beam tracking after the specified number of minutes
```

### Example

```
python phases.py --phase_centre "10h05m0s -26d0m0s" --stop_tracking_after_x_mins 35
```
