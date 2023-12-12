# analyze.py

A simple Python script to analyze the performance of a Blended-Wing Body aircraft.

For use during the MIT 16.100 Final Project.

# Cross-Platform Compatibility

`analyze.py` works on Windows and Linux!

# Usage

```
python3 analyze.py [AVL-file] [cylinders-file] [conditions-file] [lift-drag-ratio]
```

```
python3 analyze.py [AVL-file] [cylinders-file] [conditions-file] [lift-drag-ratio] --thick-scale 1.0
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

## Lift-Drag Ratio

This is a *required* argument specifying the lift-drag ratio, as calculated by CFD.

## Thickness Scale

This is an *optional* argument specifying the thickness scale of the wing.

If not specified, this defaults to `1.0`.

# Changelog

**v0.4.0**: Removed aerodynamic analysis capability due to transonic behaviors.

**v0.3.0**: Updated `analyze.py` to allow for thickness scaling.  Added `analyzewin.py` for Windows users.

**v0.2.0**: Updated `analyze.py` to allow for calculation of objective and range sensitivities.

**v0.1.0**: Initial commit.

# Thanks

Thanks to MIT Professor Qiqi Wang for:
- `design.py`
- `graph.py`
- `vehicle.py`
- `*.dat` airfoils
- `cylinders.txt`

Thanks to Boeing and Bob Liebeck for the original `bwb.avl`.