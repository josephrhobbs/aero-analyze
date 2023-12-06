# analyze.py

A simple Python script to analyze the performance of a Blended-Wing Body aircraft.

For use during the MIT 16.100 Final Project.

# Usage

```
python3 analyze.py [AVL-file] [cylinders-file] [conditions-file]
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

# Changelog

**v0.1.0**: Initial commit.

# Thanks

Thanks to MIT Professor Qiqi Wang for:
- `design.py`
- `graph.py`
- `vehicle.py`
- `*.dat` airfoils

Thanks to MIT Professor Mark Drela for:
- AVL (not included)
- XFOIL (not included)