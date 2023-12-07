# analyze.py

A simple Python script to analyze the performance of a Blended-Wing Body aircraft.

For use during the MIT 16.100 Final Project.

# Cross-Platform Compatibility

`analyze.py` works on both Windows and Linux now!

## Windows Users

Windows users should use `analyzewin.py`.

## Linux Users

Linux users should use `analyze.py`.

# Usage

```
python3 analyze.py [AVL-file] [cylinders-file] [conditions-file]
```

```
python3 analyze.py [AVL-file] [cylinders-file] [conditions-file] --thick-scale 1.0
```

## AVL File

A *metric* AVL file is required.

## Cylinders File

A *metric* cylinders file is required.

A sixth column must be included to specify the classification of each tank.
A `0` in this column indicates a passenger cabin, and a `1` indicates a fuel tank.
No other values are permitted.

## Conditions File

See `conditions.txt` for an example of this file format.

`conditions.txt` contains five values:
- Propulsion system type (`int`).  This is `1` for a turbine engine and `2` for a fuel cell.
- Cruise mach number (`float`).
- Local speed of sound (`float`, meters per second).
- Free-stream density (`float`, kilograms per cubic meter).
- Kinematic viscosity (`float`, square meters per second).

Note that all lines beginning with a hash (`#`) in a flight conditions file will be ignored. 

## Thickness Scale

This is an *optional* argument specifying the thickness scale of the wing.

If not specified, this defaults to `1.0`.

# Setup Instructions

## Windows

Before using `analyzewin.py`, copy `avl.exe` and `xfoil.exe` into the current directory.

## Linux

Before using `analyze.py`, copy `avl` and `xfoil` into the current directory.

## macOS

Unfortunately, I don't know how to use macOS well enough to provide instructions.

# Changelog

**v0.2.0**: Updated `analyze.py` to allow for calculation of objective and range sensitivities.
**v0.1.0**: Initial commit.

# Thanks

Thanks to MIT Professor Qiqi Wang for:
- `design.py`
- `graph.py`
- `vehicle.py`
- `*.dat` airfoils
- `cylinders.txt`

Thanks to MIT Professor Mark Drela for:
- AVL (not included)
- XFOIL (not included)

Thanks to Boeing and Bob Liebeck for the original `bwb.avl`.