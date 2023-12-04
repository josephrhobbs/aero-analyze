import argparse
import numpy as np
from graph import AvlInput
from math import sqrt, log, floor
import os

###############
# AVL & XFOIL #
###############

AVL = "./avl"
XFOIL = "xfoil"

#############
# DUMP FILE #
#############

DUMP = "temp.txt"

#############
# CONSTANTS #
#############

TARGET_RANGE = 14_140_000 # m
PAYLOAD = 53_000 # kg
RHO_H2 = 71 # kg / m3
D0 = 3 # m
PI = 3.141592653 # ul
PHI = 0.05 # ul
GRAVITY = 9.81 # m / s2
ENERGY = 141.86E+6 # J / kg
PAX_AREA = 0.8 # m2
PAX_MASS = 176.67 # kg

C1 = 0.003 # 1/m
C2 = 12 # kg m2
C3 = 0.2 # kg N^-0.89
C4 = 0.015 # kg / N
C5 = 0.0006 # kg / W

C6 = C4 + 100 * C5 # kg / N

def section_spar(sec, xSpar):
    x = (xSpar - sec.xyzLE[0]) / sec.chord
    iLE = (sec.airfoil[1:,0] < sec.airfoil[:-1,0]).sum()
    zUpper = np.interp(x, sec.airfoil[iLE::-1,0], sec.airfoil[iLE::-1,1])
    zLower = np.interp(x, sec.airfoil[iLE:,0], sec.airfoil[iLE:,1])
    return zLower * sec.chord + sec.xyzLE[2], zUpper * sec.chord + sec.xyzLE[2]

def compute_spar(wing, cylinders):
    xLE = np.array([s.xyzLE[0] for s in wing.sections])
    yLE = np.array([s.xyzLE[1] for s in wing.sections])
    chord = np.array([s.chord for s in wing.sections])
    xSparRoot = xLE[0] + chord[0] * 0.4
    xSparTip = xLE[-1] + chord[-1] * 0.4
    xSpar = xSparRoot + (xSparTip - xSparRoot) * yLE / yLE.max()
    zLower, zUpper = np.transpose([section_spar(s, x) for s, x in zip(wing.sections, xSpar)])
    
    b = 2 * (yLE[-1] - yLE[0])
    t = zUpper[0] - zLower[0]
    L = pow(pow(yLE[-1] - yLE[0], 2) + pow(xLE[-1] - xLE[0], 2), 0.5)

    return b, L, t

def compute_planform(avl_fname):
    iSec = 0
    avlInput = AvlInput(args.avl_fname)
    wing = avlInput.surfaces[0]
    with open(avl_fname) as fin:
        sections = []
        for line in fin:
            if line.strip() == '!Sref     Cref     Bref':
                fin.readline()
            elif line.strip() == '!     Xle      Yle      Zle    Chord     Ainc':
                sec = wing.sections[iSec]
                nextsecY = wing.sections[iSec + 1].xyzLE[1] if iSec + 1 < len(wing.sections) else 0
                nextsecX = wing.sections[iSec + 1].xyzLE[0] if iSec + 1 < len(wing.sections) else 0
                s = (nextsecY - sec.xyzLE[1]) * (nextsecX + sec.xyzLE[0])
                sections.append(s)
                fin.readline()
                iSec += 1
    wing = avlInput.surfaces[0]
    return sum(sections), sections

def compute_reference_length(avl_fname):
    _S, sections = compute_planform(avl_fname)
    return sections[0]

def fuel_cylinder_masses(diameter, length):
    """
    Returns (empty tank mass, fuel mass) for a tank
        of the provided dimensions
    """
    eta = pow(diameter, 2) / sqrt(pow(diameter, 4) + pow(D0, 4))
    m_fuel = RHO_H2 * (1/4 * PI * pow(diameter, 2) * (length - diameter) + 1/6 * PI * pow(diameter, 3))
    m_empty = (1/eta - 1) * m_fuel
    return (m_empty, m_fuel)

def compute_masses(avl_fname, cylinders_fname, proptype):
    avlInput = AvlInput(args.avl_fname)

    cylinders = np.loadtxt(args.cylinders_fname, skiprows=1)
    wing = avlInput.surfaces[0]

    # Find planform area, spar length, span, and spar root thickness
    S, _sections = compute_planform(avl_fname)
    b, L, t = compute_spar(wing, cylinders)

    # Find total cabin area
    total_cabin_area = 0
    for row in cylinders:
        if row[5] == 0:
            total_cabin_area += row[2] * (row[4] - row[3])
    pax = total_cabin_area / PAX_AREA

    # Find position, type, and masses of each cylinder
    cyl = []
    for row in cylinders:
        yi = abs(row[0])
        cyltype = "fuel" if row[5] else "cabin"
        empty_mass = 0
        landing_mass = 0
        takeoff_mass = 0
        if cyltype == "cabin":
            empty_mass = 0
            landing_mass = takeoff_mass = row[2] * (row[4] - row[3]) / PAX_AREA * PAX_MASS
        else:
            (me, mf) = fuel_cylinder_masses(row[2], row[4] - row[3])
            empty_mass = me
            landing_mass = me + PHI * mf
            takeoff_mass = me + mf
        cyl.append((1 - 2*yi/b, empty_mass, landing_mass, takeoff_mass, cyltype))
    
    # COMPUTE STRUCTURAL MASS
    m_struct = C2 * S
    for k, me, ml, mto, cyltype in cyl:
        m_struct += C1 * pow(L, 2) / t * k * mto

    # COMPUTE PROPULSION SYSTEM MASS
    fuel_mto = 0
    for k, me, ml, mto, cyltype in cyl:
        fuel_mto += mto
    appx_takeoff_mass = 1.2 * (fuel_mto + PAYLOAD + m_struct)
    m_turbine = 2 * C3 * pow(appx_takeoff_mass * GRAVITY / 10, 0.89)
    m_fuelcell = appx_takeoff_mass * GRAVITY / 5 * C6

    # GET PROPULSION SYSTEM TYPE
    if proptype == 1:
        m_propulsion = m_turbine
        eta = 0.4
    else:
        m_propulsion = m_fuelcell
        eta = 0.45

    # COMPUTE TOTAL MASSES
    empty_mass = landing_mass = takeoff_mass = m_struct + PAYLOAD
    for k, me, ml, mto, cyltype in cyl:
        empty_mass += me
        landing_mass += ml
        takeoff_mass += mto

    return (empty_mass, landing_mass, takeoff_mass, pax, eta)

def compute_range(ml, mto, eta, cl, cd):
    return ENERGY * eta / GRAVITY * cl / cd * log(mto / ml)

def process_conditions(conditions_fname):
    """
    Read in following conditions:
        - Propulsion system type
        - Cruise Mach number
        - Local speed of sound
        - Free-stream density
        - Free-stream kinematic viscosity
    """
    with open(conditions_fname, "r") as f:
        data = [float(row) for row in f.read().split("\n") if len(row) > 0 and row[0] != "#"]
    assert len(data) == 5
    return data

def viscous_drag_array(re, aoa):
    print("VISCOUS DRAG MODEL")
    coefficients = []
    for i in range(15):
        airfoil = f"orig{i}foilmod.dat"
        X = f"""load {airfoil}

oper
iter 500
visc {re}
alfa 0
alfa {float(aoa/2)}
pacc
save.txt
dump.txt
alfa {float(aoa)}

quit
            """
        with open("x", "w") as f:
            f.write(X)
        os.system(f"{XFOIL} < x > {DUMP}")
        cdv = 0
        try:
            with open("save.txt", "r") as f:
                file = f.read()
                cdv = float(file.split("\n")[-2].strip().split()[2])
            os.remove("save.txt")
            os.remove("dump.txt")
            os.remove("x")
        except:
            print("XFOIL calculation diverged :(")
            with open("x", "r") as f:
                print("INPUT FILE")
                print(f.read())
                print()
            with open("save.txt", "r") as f:
                print("OUTPUT FILE")
                print(f.read())
                print()
            os.remove("save.txt")
            os.remove("dump.txt")
            os.remove("x")
            return [0] * 15
        coefficients.append(cdv)
    return coefficients

def compute_aoa_and_induced_drag(cl):
    print("INVISCID VORTEX LATTICE MODEL")
    A = f"""load bwb.avl

oper
A
C {cl}
X
W
forces.txt
W
test.txt

quit
    """
    with open("a", "w") as f:
        f.write(A)
    os.system(f"{AVL} < a > {DUMP}")
    try:
        with open("forces.txt", "r") as f:
            data = [[v.strip() for v in d.strip().split("=") if len(v.strip()) != 0] for d in f.read().split()]
        os.remove("forces.txt")
        os.remove("test.txt")
        os.remove("a")
    except:
        print("AVL calculation diverged :(")
        with open("a", "r") as f:
            print("INPUT FILE")
            print(f.read())
            print()
        with open("forces.txt", "r") as f:
            print("OUTPUT FILE")
            print(f.read())
            print()
        os.remove("forces.txt")
        os.remove("test.txt")
        os.remove("a")
        return 0, 0
    processed = []
    for row in data:
        if len(row) > 0:
            processed.append(row[0])
    aoa = float(processed[1 + processed.index("Alpha")])
    cdi = float(processed[1 + processed.index("CDind")])
    return aoa, cdi

def compute_viscous_drag(cdv_array, sections):
    avlInput = AvlInput(args.avl_fname)

    viscous_drag = 0

    for i in range(15):
        viscous_drag += cdv_array[i] * sections[i]

    return viscous_drag / sum(sections)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='substitute',
            description='substitute new planform into AVL file')
    parser.add_argument('avl_fname')
    parser.add_argument('cylinders_fname')
    parser.add_argument('conditions_fname')
    args = parser.parse_args()

    # READ IN FLIGHT CONDITIONS

    print(f"FLIGHT CONDITIONS")
    proptype, mach, ainf, rho, kinvisc = process_conditions(args.conditions_fname)
    uinf = mach * ainf
    l = compute_reference_length(args.avl_fname)
    re = uinf * l / kinvisc
    print(f"Uinf = {round(uinf, 4)} m/s")
    print(f"Re = {round(re / 1E+6, 4)}E+6")

    # COMPUTE FLIGHT TIME

    time = TARGET_RANGE / uinf
    print(f"Time = {round(time / 3_600, 4)} h")
    print()

    # COMPUTE MASSES AND PLANFORM FROM INPUT FILES

    print("AIRCRAFT MASS AND PASSENGER CAPACITY")
    me, ml, mto, pax, eta = compute_masses(args.avl_fname, args.cylinders_fname, proptype)
    S, sections = compute_planform(args.avl_fname)
    print(f"Empty mass = {round(me, 4)} kg")
    print(f"Landing mass = {round(ml, 4)} kg")
    print(f"Take-off mass = {round(mto, 4)} kg")
    print(f"Passengers = {floor(pax)} pax")
    print()

    # COMPUTE OBJECTIVE

    print("AIRCRAFT OBJECTIVE")
    obj = time * me
    print(f"OBJ = {round(obj / 1E+9, 4)}E+9 kg-s")
    print()

    # COMPUTE CL

    print("COEFFICIENT OF LIFT")
    q = 1/2 * rho * pow(uinf, 2)
    cl = mto * GRAVITY / q / S
    print(f"CL = {round(cl, 4)}")
    print()

    # COMPUTE AOA AND CDi FROM AVL

    aoa, cdi = compute_aoa_and_induced_drag(cl)
    print(f"AOA = {round(aoa, 4)}")
    print(f"CDi = {round(cdi, 4)}")
    print()

    # COMPUTE CDv FROM XFOIL

    cdv_array = viscous_drag_array(re, aoa)
    cdv = compute_viscous_drag(cdv_array, sections)
    print(f"CDv = {round(sum(cdv_array), 4)}")
    print()

    # COMPUTE CD FROM CDi AND CDv

    cd = cdv + cdi

    # COMPUTE RANGE

    print("AIRCRAFT RANGE")
    R = compute_range(ml, mto, eta, cl, cd)
    print(f"Range = {round(R / 1000, 4)} km")