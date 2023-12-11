import argparse
import numpy as np
from graph import AvlInput
from math import sqrt, log, floor, cos, acos
import os

###############
# AVL & XFOIL #
###############

PWD = os.getcwd()

AVL = rf"{PWD}/avl"
XFOIL = rf"{PWD}/xfoil"

#############
# CONSTANTS #
#############

TARGET_RANGE = 14_140_000 # m
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

# Wave drag model coefficients
W1 = -0.028143
W1 =  0.674524

def section_spar(sec, xSpar):
    x = (xSpar - sec.xyzLE[0]) / sec.chord
    iLE = (sec.airfoil[1:,0] < sec.airfoil[:-1,0]).sum()
    zUpper = np.interp(x, sec.airfoil[iLE::-1,0], sec.airfoil[iLE::-1,1])
    zLower = np.interp(x, sec.airfoil[iLE:,0], sec.airfoil[iLE:,1])
    return zLower * sec.chord + sec.xyzLE[2], zUpper * sec.chord + sec.xyzLE[2]

def compute_spar(avl_fname, cylinders_fname):
    avlInput = AvlInput(avl_fname)
    cylinders = np.loadtxt(cylinders_fname, skiprows=1)
    wing = avlInput.surfaces[0]

    xLE = np.array([s.xyzLE[0] for s in wing.sections])
    yLE = np.array([s.xyzLE[1] for s in wing.sections])
    chord = np.array([s.chord for s in wing.sections])
    xSparRoot = xLE[0] + chord[0] * 0.4
    xSparTip = xLE[-1] + chord[-1] * 0.4
    xSpar = xSparRoot + (xSparTip - xSparRoot) * yLE / yLE.max()
    zLower, zUpper = np.transpose([section_spar(s, x) for s, x in zip(wing.sections, xSpar)])
    
    b = 2 * (yLE[-1] - yLE[0])
    t = zUpper[0] - zLower[0]
    L = pow(pow(yLE[-1] - yLE[0], 2) + pow(xSparTip - xSparRoot, 2), 0.5)

    theta = acos(b/2 / L)

    return b, theta, t

def compute_planform(avl_fname):
    avlInput = AvlInput(avl_fname)
    wing = avlInput.surfaces[0]
    sections = []
    for i in range(len(wing.sections)):
         c = wing.sections[i].chord
         y0 = wing.sections[max(0, i-1)].xyzLE[1]
         y1 = wing.sections[min(len(wing.sections)-1, i+1)].xyzLE[1]
         dy = (y1 - y0) / 2
         s = 2 * c * dy
         sections.append(s)
    return sum(sections), sections

def compute_reference_dimensions(avl_fname):
    avlInput = AvlInput(avl_fname)
    wing = avlInput.surfaces[0]
    cref = avlInput.Cref
    bref = 2 * wing.sections[-1].xyzLE[1]
    return bref, cref

def get_cref(avl_fname):
    avlInput = AvlInput(avl_fname)
    return avlInput.Cref

def fuel_cylinder_masses(diameter, length):
    """
    Returns (empty tank mass, fuel mass) for a tank
        of the provided dimensions
    """
    eta = pow(diameter, 2) / sqrt(pow(diameter, 4) + pow(D0, 4))
    m_fuel = RHO_H2 * (1/4 * PI * pow(diameter, 2) * (length - diameter) + 1/6 * PI * pow(diameter, 3))
    m_empty = (1/eta - 1) * m_fuel
    return (m_empty, m_fuel)

def compute_masses(avl_fname, cylinders_fname, proptype, b, theta, t, verbose=False):
    avlInput = AvlInput(avl_fname)

    cylinders = np.loadtxt(cylinders_fname, skiprows=1)
    wing = avlInput.surfaces[0]

    # Find planform area and spar length
    S, _sections = compute_planform(avl_fname)
    L = b/2 / cos(theta)

    # Find total cabin area
    total_cabin_area = 0
    for row in cylinders:
        if row[5] == 0:
            total_cabin_area += row[2] * (row[4] - row[3])
    pax = total_cabin_area / PAX_AREA
    payload = total_cabin_area * 53_000 / 240

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
    appx_takeoff_mass = 1.2 * (fuel_mto + payload + m_struct)
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
    fuel_mass = 0
    tank_mass = 0
    empty_mass = landing_mass = takeoff_mass = m_struct + payload + m_propulsion
    for k, me, ml, mto, cyltype in cyl:
        empty_mass += me
        landing_mass += ml
        takeoff_mass += mto
        tank_mass += me
        fuel_mass += mto - me

    if verbose:
        print("MASS BREAKDOWN")
        print(f"Structural mass = {round(m_struct, 4)} kg")
        print(f"Payload mass = {round(payload, 4)} kg")
        print(f"Propulsion system mass = {round(m_propulsion, 4)} kg")
        print(f"Fuel mass = {round(fuel_mass, 4)} kg")
        print(f"Empty tank mass = {round(tank_mass, 4)} kg")
        print()

    return (empty_mass, landing_mass, takeoff_mass, pax, payload, eta)

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

def viscous_drag_array(re, aoa, mach, thick_scale, verbose=False):
    if verbose:
        print("VISCOUS DRAG MODEL")
    coefficients = []
    for i in range(15):
        airfoil = f"orig{i}foilmod.dat"
        hard = "hard" if i == 0 else ""
        X = f"""load {airfoil}

gdes
tfac {thick_scale} {thick_scale}
exec

oper
visc {re}
mach {mach}
alfa 0
alfa {float(aoa/2)}
alfa {float(aoa)}
{hard}

quit
            """
        with open("x", "w") as f:
            f.write(X)
        os.system(f"{XFOIL} < x > save.txt")
        cdv = 0
        try:
            with open("save.txt", "r") as f:
                file = f.read()
                cdv = float(file.split("CD = ")[-1].strip()[:6])
            os.remove("save.txt")
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
            os.remove("x")
            return [0] * 15
        coefficients.append(cdv)
    return coefficients

def compute_aoa_and_induced_drag(avl_fname, cl, verbose=False):
    if verbose:
        print("INVISCID VORTEX LATTICE MODEL")
    A = f"""load {avl_fname}

oper
A
C {cl}
X

quit
    """
    with open("a", "w") as f:
        f.write(A)
    os.system(f"{AVL} < a > save.txt")
    aoa = 0
    cdi = 0
    try:
        with open("save.txt", "r") as f:
            file = f.read()
            aoa = float(file.split("Alpha =  ")[-1].strip()[:8])
            cdi = float(file.split("CDind = ")[-1].strip()[:8])
        os.remove("save.txt")
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
        os.remove(DUMP)
        os.remove("a")
    return aoa, cdi

def compute_viscous_drag(cdv_array, sections):
    avlInput = AvlInput(args.avl_fname)

    viscous_drag = 0

    for i in range(15):
        viscous_drag += cdv_array[i] * sections[i]

    return viscous_drag / sum(sections)

def compute_wave_drag(mach, aoa):
    # Compute critical Mach number
    Mcr = W1 * aoa + W2

    # Compute coefficient of wave drag
    if mach >= Mcr:
        return 20 * pow(mach - Mcr, 4)
    else:
        return 0

def analyze(args, b, theta, t, cref, thick_scale, verbose=False):
    # READ IN FLIGHT CONDITIONS

    proptype, mach, ainf, rho, kinvisc = process_conditions(args.conditions_fname)
    uinf = mach * ainf
    bref = b
    Sref = bref * cref
    re = uinf * cref / kinvisc

    if verbose:
        print("FLIGHT CONDITIONS")
        print(f"bref = {bref} m")
        print(f"cref = {cref} m")
        print(f"Sref = {round(Sref, 4)} m2")
        print(f"Uinf = {round(uinf, 4)} m/s")
        print(f"Re = {round(re / 1E+6, 4)}E+6")

    # COMPUTE FLIGHT TIME

    time = TARGET_RANGE / uinf

    if verbose:
        print(f"Time = {round(time / 3_600, 4)} h")
        print()

    # COMPUTE MASSES AND PLANFORM FROM INPUT FILES

    me, ml, mto, pax, payload, eta = compute_masses(args.avl_fname, args.cylinders_fname, proptype, b, theta, t, verbose)
    S, sections = compute_planform(args.avl_fname)

    if verbose:
        print("AIRCRAFT MASS")
        print(f"Empty mass = {round(me, 4)} kg")
        print(f"Landing mass = {round(ml, 4)} kg")
        print(f"Take-off mass = {round(mto, 4)} kg")
        print()
        print("PASSENGER CAPACITY")
        print(f"Passengers = {floor(pax)} pax")
        print()

    # COMPUTE OBJECTIVE

    obj = time * me / payload

    if verbose:
        print("AIRCRAFT OBJECTIVE")
        print(f"OBJ A = {round(obj * payload / 1E+9, 4)}E+9 kg s")
        print(f"OBJ B = {round(obj / 3_600.0, 4)} h")
        print()

    # COMPUTE CL

    q = 1/2 * rho * pow(uinf, 2)
    cl = mto * GRAVITY / q / Sref

    if verbose:
        print("COEFFICIENT OF LIFT")
        print(f"CL = {round(cl, 4)}")
        print()

    # COMPUTE AOA AND CDi FROM AVL

    aoa, cdi = compute_aoa_and_induced_drag(args.avl_fname, cl, verbose)

    if verbose:
        print(f"AOA = {round(aoa, 4)}")
        print(f"CDi = {round(cdi, 4)}")
        print()

    # COMPUTE CDv FROM XFOIL

    cdv_array = viscous_drag_array(re, aoa, mach, thick_scale, verbose)
    cdv = compute_viscous_drag(cdv_array, sections)

    if verbose:
        print(f"CDv = {round(cdv, 4)}")
        print()

    # COMPUTE CDw FROM MODEL

    cdw = compute_wave_drag(mach, aoa)

    # COMPUTE CD FROM CDi AND CDv

    cd = cdv + cdi

    if verbose:
        print("COEFFICIENT OF DRAG")
        print(f"CD = {round(cd, 4)}")
        print()

    # COMPUTE RANGE

    R = compute_range(ml, mto, eta, cl, cd)
    
    if verbose:
        print("AIRCRAFT RANGE")
        print(f"L/D = {round(cl/cd, 4)}")
        print(f"Range = {round(R / 1000, 4)} km")

    return obj, R

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='substitute',
            description='substitute new planform into AVL file')
    parser.add_argument('avl_fname')
    parser.add_argument('cylinders_fname')
    parser.add_argument('conditions_fname')
    parser.add_argument('--thick-scale', type=float, default=1.0)
    args = parser.parse_args()

    b, theta, t = compute_spar(args.avl_fname, args.cylinders_fname)
    cref = get_cref(args.avl_fname)

    obj, R = analyze(args, b, theta, args.thick_scale*t, cref, args.thick_scale, verbose=True)