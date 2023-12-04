import os
import argparse
import numpy as np
import pylab
import dataclasses
import vehicle

def section_spar(sec, xSpar):
    x = (xSpar - sec.xyzLE[0]) / sec.chord
    iLE = (sec.airfoil[1:,0] < sec.airfoil[:-1,0]).sum()
    zUpper = np.interp(x, sec.airfoil[iLE::-1,0], sec.airfoil[iLE::-1,1])
    zLower = np.interp(x, sec.airfoil[iLE:,0], sec.airfoil[iLE:,1])
    return zLower * sec.chord + sec.xyzLE[2], zUpper * sec.chord + sec.xyzLE[2]

def naca4(number, n=100):
    """
    Returns 2*n+1 points in [0 1] for the given 4 digit NACA number string
    """
    assert len(number) == 4
    if number[2:] <= '06':
        number = number[:2] + '06'

    m = float(number[0])/100.0
    p = float(number[1])/10.0
    t = float(number[2])/10.0 + float(number[3])/100.0

    a0 = +0.2969
    a1 = -0.1260
    a2 = -0.3516
    a3 = +0.2843
    a4 = -0.1015 # For finite thick TE

    beta = np.linspace(0.0,np.pi,n+1)
    x = [(0.5*(1.0-np.cos(xx))) for xx in beta]  # Half cosine based spacing

    yt = [5*t*(a0*np.sqrt(xx)+a1*xx+a2*np.power(xx,2)+a3*np.power(xx,3)+a4*np.power(xx,4)) for xx in x]

    xc1 = [xx for xx in x if xx <= p]
    xc2 = [xx for xx in x if xx > p]

    if p == 0:
        xu = x
        yu = yt

        xl = x
        yl = [-xx for xx in yt]

        xc = xc1 + xc2
        zc = [0]*len(xc)
    else:
        yc1 = [m/np.power(p,2)*xx*(2*p-xx) for xx in xc1]
        yc2 = [m/np.power(1-p,2)*(1-2*p+xx)*(1-xx) for xx in xc2]
        zc = yc1 + yc2

        dyc1_dx = [m/np.power(p,2)*(2*p-2*xx) for xx in xc1]
        dyc2_dx = [m/np.power(1-p,2)*(2*p-2*xx) for xx in xc2]
        dyc_dx = dyc1_dx + dyc2_dx

        theta = [np.arctan(xx) for xx in dyc_dx]

        xu = [xx - yy * np.sin(zz) for xx,yy,zz in zip(x,yt,theta)]
        yu = [xx + yy * np.cos(zz) for xx,yy,zz in zip(zc,yt,theta)]

        xl = [xx + yy * np.sin(zz) for xx,yy,zz in zip(x,yt,theta)]
        yl = [xx - yy * np.cos(zz) for xx,yy,zz in zip(zc,yt,theta)]

    X = xu[::-1] + xl[1:]
    Z = yu[::-1] + yl[1:]

    return np.transpose([X,Z])

class AvlSection(vehicle.AirfoilSection):
    def __init__(self, reader, path):
        self.path = path
        section = reader.readFloats()
        assert len(section) == 5 or len(section) == 7
        x, y, z, c, a = section[:5]
        self.xyzLE = (x, y, z)
        self.chord = c
        self.aInc = np.deg2rad(a)
        if len(section) == 7:
            self.nSpan, self.sSpace = section[5:]
        self.controls = {}

        keywords = ['NACA', 'AIRFOIL', 'AFIL', 'AFILE', 'DESIGN', 'CONTROL', 'CLAF', 'CDCL', ]
        while reader.nextLine and reader.nextLine.split()[0] in keywords:
            method = getattr(self, f'read_{reader.readline().split()[0]}')
            method(reader)

        if not hasattr(self, 'airfoil'):
            self.airfoil = naca4('0012')

    def read_CONTROL(self, reader):
        name, control = reader.readStringFloats()
        self.controls[name] = control

    def read_CLAF(self, reader):
        reader.readFloat()

    def read_CDCL(self, reader):
        reader.readFloat()

    def read_AIRFOIL(self, reader):
        xy = []
        while (len(reader.nextLine.strip().split())) == 2:
            xy.append(reader.readFloats(2))
        self.airfoil = np.array(xy)

    def read_DESIGN(self, reader):
        reader.readline()

    def read_AFILE(self, reader):
        line = reader.readline()
        if len(line.split()) == 1:
            fname = line.strip()
        else:
            assert len(line.split()) == 3
            fname = line.split()[0].strip()
        if fname.startswith('"') and fname.endswith('"'):
            fname = fname[1:-1]
        fname = fname if os.path.exists(fname) else os.path.join(self.path, fname)
        self.airfoil = np.loadtxt(fname, skiprows=1)

    def read_AFIL(self, reader):
        self.read_AFILE(reader)

    def read_NACA(self, reader):
        line = reader.readline().strip()
        self.airfoil = naca4(line)

def AvlSurface(reader, path) -> tuple[vehicle.Surface]:
    def read_INDEX(s, reader, path):
        reader.readline()

    def read_CDCL(s, reader, path):
        reader.readline()

    def read_YDUPLICATE(s, reader, path):
        s['yDuplicate'] = reader.readFloat()

    def read_ANGLE(s, reader, path):
        s['angle'] = reader.readFloat()

    def read_AINC(s, reader, path):
        read_ANGLE(s, reader)

    def read_TRANSLATE(s, reader, path):
        s['translate'] = reader.readFloats(3)

    def read_SCALE(s, reader, path):
        s['scale'] = reader.readFloats(3)

    def read_NOWAKE(s, reader, path):
        pass

    def read_SECTION(s, reader, path):
        if 'sections' not in s:
            s['sections'] = []
        s['sections'].append(AvlSection(reader, path))

    def read_COMPONENT(s, reader, path):
        s['component'] = int(reader.readFloat())

    line = reader.readline()
    assert line == 'SURFACE', line
    name = reader.readline().replace(' ', '_')
    vortexSpacing = reader.readFloats()

    keywords = ['COMPONENT', 'INDEX', 'YDUPLICATE',
            'SCALE', 'TRANSLATE', 'ANGLE', 'AINC', 'NOWAKE', 'NOALBE',
            'NOLOAD', 'CDCL', 'SECTION']
    surf = {}
    while reader.nextLine and reader.nextLine.split()[0] in keywords:
        method = locals()[f'read_{reader.readline().split()[0]}']
        method(surf, reader, path)

    translate = np.array(surf.get('translate', (0, 0, 0)))
    scale = np.array(surf.get('scale', (1, 1, 1)))
    angle = surf.get('angle', 0)
    nSecs = len(surf['sections'])
    assert nSecs > 1, f"{name} has {len(surf['sections'])} sections"
    aRot = [surf['sections'][i].x_rot_rad(surf['sections'][i+1])
            for i in range(nSecs-1)]
    aRot = [vehicle.mean_angle(aRot[max(0,i-1):min(nSecs,i+1)]) for i in range(nSecs)]
    sections = [vehicle.AirfoilSection(
        xyzLE = tuple(np.array(s.xyzLE) * scale + translate),
        chord = s.chord * scale[0],
        aInc = s.aInc + surf.get('angle', 0), aRot = a,
        airfoil = s.airfoil) for s, a in zip(surf['sections'], aRot)]
    surface = vehicle.Surface(name, sections)
    if vehicle.is_collocate(surface.sections[0], surface.sections[-1]):
        surfaces = vehicle.splitPeriodicSurface(surface)
    else:
        surfaces = (surface,)

    if 'yDuplicate' in surf:
        surfaces = tuple(s for s in surfaces) \
                 + tuple(vehicle.mirrorSurface(s, surf['yDuplicate'])
                         for s in surfaces)
    return surfaces


def AvlBody(reader, path) -> tuple[vehicle.CenterBody]:
    def read_YDUPLICATE(b, reader, path):
        b['yDuplicate'] = reader.readFloat()

    def read_TRANSLATE(b, reader, path):
        b['translate'] = reader.readFloats(3)

    def read_SCALE(b, reader, path):
        b['scale'] = reader.readFloats(3)

    def read_ANGLE(b, reader, path):
        b['angle'] = reader.readFloat()

    def read_BFIL(b, reader, path):
        line = reader.readline()
        if len(line.split()) == 1:
            fname = line.strip()
        else:
            assert len(line.split()) == 3
            fname = line.split()[0].strip()
        if fname.startswith('"') and fname.endswith('"'):
            fname = fname[1:-1]
        fname = fname if os.path.exists(fname) else os.path.join(path, fname)
        b['profile'] = np.loadtxt(fname, skiprows=1)

    line = reader.readline()
    assert line == 'BODY', line
    name = reader.readline().replace(' ', '_')
    vortexSpacing = reader.readFloats()

    keywords = ['YDUPLICATE', 'SCALE', 'TRANSLATE', 'BFIL', 'ANGLE']
    body = {}
    while reader.nextLine in keywords:
        method = locals()[f'read_{reader.readline()}']
        method(body, reader, path)

    return vehicle.CenterBody(name = name,
            profile = body.get('profile', naca4('0012')),
            translate = body.get('translate', (0, 0, 0)),
            angle = body.get('angle', 0),
            scale = body.get('scale', 1))

class AvlInput:
    def __init__(self, fname):
        self.path = os.path.dirname(os.path.realpath(fname))
        self.favl = open(fname)
        self.nextLine = None
        self.name = self.readline().replace(' ', '_')
        self.Mach = self.readFloat()
        self.Isyms = self.readFloats(3)
        assert self.Isyms[0] == 0
        assert self.Isyms[1] == 0
        assert self.Isyms[2] == 0
        self.Sref, self.Cref, self.Bref = self.readFloats(3)
        self.xyzref = self.readFloats(3)

        keywords = ['SURFACE', 'BODY']
        if self.nextLine not in keywords:
            self.CDp = self.readFloat()
        else:
            self.CDp = 0.0

        self.surfaces = []
        self.bodies = []
        while self.nextLine in keywords:
            method = getattr(self, f'read_{self.nextLine}')
            method()

    def read_SURFACE(self):
        self.surfaces.extend(AvlSurface(self, self.path))

    def read_BODY(self):
        self.bodies.append(AvlBody(self, self.path))

    def readline(self):
        thisLine = self.nextLine
        while True:
            line = self.favl.readline()
            if line == '':
                self.nextLine = None
                break
            line = line.strip()
            if not (line.startswith('!') or line.startswith('#') or len(line) == 0):
                if '!' in line:
                    line = line.split('!')[0]
                self.nextLine = line
                break
        if thisLine is not None:
            return thisLine
        else:
            return self.readline()

    def readStringFloats(self, n=None):
        numbers = self.readline().split()
        if n is not None:
            assert len(numbers) >= n + 1
        else:
            n = len(numbers) - 1
        return numbers[0], [float(x) for x in numbers[1:n+1]]

    def readFloats(self, n=None):
        numbers = self.readline().split()
        if n is not None:
            assert len(numbers) >= n
        else:
            n = len(numbers)
        return [float(x) for x in numbers[:n]]

    def readFloat(self):
        return self.readFloats(1)[0]

    def to_egads(self, context, config, verbosity=0):
        return {b.name : b.to_egads(context, config, verbosity)
                for b in self.surfaces + self.bodies}

def plot_cylinder(y, z, d, x0, x1, dxPlot, dxPlot0, dzPlot, color):
    y0, y1 = y - d/2, y + d/2
    z0, z1 = z - d/2, z + d/2
    dCos = np.cos(np.linspace(0, np.pi, 100)) * d/2
    dSin = np.sin(np.linspace(0, np.pi, 100)) * d/2
    xLine = np.hstack([x0, x1 + dSin, x0 - dSin])
    yLine = np.hstack([y0, y - dCos, y + dCos])
    zLine = np.hstack([z0, z - dCos, z + dCos])
    pylab.plot(xLine, yLine, color, lw=2)
    pylab.plot(dxPlot + xLine, dzPlot + zLine, color, lw=2)

    dy = np.cos(np.linspace(0, 2*np.pi, 100)) * d/2
    dz = np.sin(np.linspace(0, 2*np.pi, 100)) * d/2
    pylab.plot(dxPlot0 + z + dz, y + dy, color, lw=2)

def interpolate_sections(sec0, sec1, a1):
    airfoil = (1 - a1) * sec0.airfoil + a1 * sec1.airfoil
    xyzLE = (1 - a1) * np.array(sec0.xyzLE) + a1 * np.array(sec1.xyzLE)
    chord = (1 - a1) * sec0.chord + a1 * sec1.chord
    return dataclasses.replace(sec0,
        airfoil = airfoil,
        xyzLE = tuple(xyzLE),
        chord = chord)

def plot_section(sec, dxPlot, dzPlot, color):
    xz = sec.airfoil[:,:2] * sec.chord + np.array([sec.xyzLE[0], sec.xyzLE[2]])
    pylab.plot(dxPlot + xz[:,0], dzPlot + xz[:,1], color)

def plot_all(wing, cylinders, root_spar_x_over_c, tip_spar_x_over_c):
    xLE = np.array([s.xyzLE[0] for s in wing.sections])
    yLE = np.array([s.xyzLE[1] for s in wing.sections])
    zLE = np.array([s.xyzLE[2] for s in wing.sections])
    chord = np.array([s.chord for s in wing.sections])
    zMax = zLE + np.array([s.airfoil[:,1].max() * s.chord for s in wing.sections])
    zMin = zLE + np.array([s.airfoil[:,1].min() * s.chord for s in wing.sections])

    pylab.plot(xLE, yLE, 'k')
    pylab.plot(xLE + chord, yLE, 'k')
    pylab.plot([xLE[-1], xLE[-1] + chord[-1]], [yLE[-1], yLE[-1]], 'k')
    pylab.plot(xLE, -yLE, 'k')
    pylab.plot(xLE + chord, -yLE, 'k')
    pylab.plot([xLE[-1], xLE[-1] + chord[-1]], [-yLE[-1], -yLE[-1]], 'k')

    pylab.plot([xLE, xLE + chord], [yLE, yLE], ':k')

    xSparRoot = xLE[0] + chord[0] * root_spar_x_over_c
    xSparTip = xLE[-1] + chord[-1] * tip_spar_x_over_c
    xSpar = xSparRoot + (xSparTip - xSparRoot) * yLE / yLE.max()
    pylab.plot(xSpar, yLE, '--k', lw=2)
    pylab.plot(xSpar, -yLE, '--k', lw=2)

    dxPlot0 = xLE.min() - chord.max() / 4 - zMax.max()
    pylab.plot(dxPlot0 + zMax, yLE, 'k')
    pylab.plot(dxPlot0 + zMin, yLE, 'k')
    pylab.plot(dxPlot0 + zMax, -yLE, 'k')
    pylab.plot(dxPlot0 + zMin, -yLE, 'k')

    zLower, zUpper = np.transpose([section_spar(s, x) for s, x in zip(wing.sections, xSpar)])
    pylab.plot(dxPlot0 + zLower, yLE, '--k')
    pylab.plot(dxPlot0 + zUpper, yLE, '--k')
    pylab.plot(dxPlot0 + zLower, -yLE, '--k')
    pylab.plot(dxPlot0 + zUpper, -yLE, '--k')

    pylab.plot([dxPlot0 + zLower, dxPlot0 + zUpper], [yLE, yLE], ':r')

    dxPlot = 1.2 * ((xLE + chord).max() - xLE.min())
    dzPlot = (zMax - zMin).max() * 1.2 * np.arange(len(wing.sections)) - yLE.max() * 0.9
    colors = ['r', 'g', 'b', 'c', 'm', 'y']

    print('Spar_X         Spar_Y        Spar_Height')
    for sec, xs, dzp, zl, zu in zip(wing.sections, xSpar, dzPlot, zLower, zUpper):
        plot_section(sec, dxPlot, dzp, 'k')
        pylab.plot([dxPlot + xs, dxPlot + xs], [zl + dzp, zu + dzp], '--k', lw=2)
        print(f'{xs:6.2f}        {sec.xyzLE[1]:6.2f}          {zu-zl:7.4f}')

    dxPlot += 1.2 * chord.max()
    dzPlot = (zMax - zMin).max() * 1.5 * np.arange(len(cylinders)) - yLE.max() * 0.9
    for i, (y, z, d, x0, x1, _) in enumerate(cylinders):
        color = colors[i % len(colors)]
        plot_cylinder(y, z, d, x0, x1, dxPlot, dxPlot0, dzPlot[i], color)
        iSec = np.interp(abs(y), yLE, np.arange(yLE.size))
        sec = wing.sections[int(iSec)] if i % 1 == 0 \
                else interpolate_sections(wing.sections[int(iSec)], wing.sections[int(iSec)+1], iSec%1)
        plot_section(sec, dxPlot, dzPlot[i], 'k')

    pylab.axis('scaled')
    pylab.grid()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='graph',
            description='Fit cylinders into wing')
    parser.add_argument('input_avl_fname')
    parser.add_argument('cylinders_fname')

    args = parser.parse_args()
    avlInput = AvlInput(args.input_avl_fname)
    cylinders = np.loadtxt(args.cylinders_fname, skiprows=1)

    plot_all(avlInput.surfaces[0], cylinders, 0.4, 0.4)
    pylab.show()
