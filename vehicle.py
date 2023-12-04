import random
import itertools
import collections
import dataclasses
import pylab
import numpy as np

Config = collections.namedtuple('Config', [
    'max_span_tethickness_ratio',
    'max_body_butt_ratio',
    'wing_tip'])

def thicken_te(xz, target_thickness):
    xz = np.array(xz).copy()
    thickness = xz[0,1] - xz[-1,1]
    increment = max(0, target_thickness - thickness)
    is_upper = np.hstack([xz[1:,0] < xz[:-1,0], [False]])
    is_lower = np.hstack([[False], xz[1:,0] > xz[:-1,0]])
    zeta = (xz[:,0] - xz[:,0].min()) / (xz[:,0].max() - xz[:,0].min())
    xz[is_upper,1] += increment * zeta[is_upper]**2 / 2
    xz[is_lower,1] -= increment * zeta[is_lower]**2 / 2
    return xz

def getLEidx(xyz):
    TE = 0.5*(xyz[0] + xyz[-1])
    d2 = ((xyz - TE)**2).sum(1)
    return d2.argmax()

def egads_line(context, node0, node1):
    nodes = [node0, node1]
    points = [np.array(n.getTopology()[3]) for n in nodes]
    direction = points[1] - points[0]
    data = list(points[0]) + list(direction)
    line = context.makeGeometry(egads.CURVE, egads.LINE, data)
    length = np.sqrt((direction**2).sum())
    return context.makeTopology(egads.EDGE, egads.TWONODE, geom=line,
                                children=nodes, reals=[0,length])

@dataclasses.dataclass
class AirfoilSection:
    '''
    An airfoil section
    '''
    xyzLE : tuple[float]
    chord : float
    aInc : float
    airfoil : np.ndarray
    aRot : float = 0.0

    def __repr__(self):
        return f'Section at {self.xyzLE}, chord {self.chord}, aInc {self.aInc}, aRot {self.aRot}, airfoil with {self.airfoil.shape[0]} points\n'

    def xyz(self, te_thickness):
        xz = thicken_te(self.airfoil, te_thickness / self.chord) * self.chord
        c, s = np.cos(self.aInc), np.sin(self.aInc)
        xz = np.transpose([xz[:,0] * c + xz[:,1] * s,
                           xz[:,1] * c - xz[:,0] * s])
        c, s = np.cos(self.aRot), np.sin(self.aRot)
        xyz = np.transpose([xz[:,0],
                           -xz[:,1] * s,
                           +xz[:,1] * c])
        return xyz + np.array(self.xyzLE)

    def x_rot_rad(self, neighbor):
        dy = neighbor.xyzLE[1] - self.xyzLE[1]
        dz = neighbor.xyzLE[2] - self.xyzLE[2]
        return np.arctan2(dz, dy)

    def to_egads(self, context, te_thickness, verbosity=0):
        xyz = self.xyz(te_thickness)
        bspline = context.approximate([len(xyz), 0], xyz)
        tLE = bspline.invEvaluate(self.xyzLE)[0]
        points = np.array([bspline.evaluate(t)[0] for t in [0, tLE, 1]])
        nodes = [context.makeTopology(egads.NODE, reals=list(p)) for p in points]
        upper = context.makeTopology(egads.EDGE, egads.TWONODE, geom=bspline,
                                     children=[nodes[0], nodes[1]], reals=[0,tLE])
        lower = context.makeTopology(egads.EDGE, egads.TWONODE, geom=bspline,
                                     children=[nodes[1], nodes[2]], reals=[tLE,1])
        te = egads_line(context, nodes[0], nodes[2])
        loop = context.makeTopology(egads.LOOP, egads.CLOSED, children=[upper, lower, te],
                                    senses=[egads.SFORWARD, egads.SFORWARD, egads.SREVERSE])
        c, s = np.cos(self.aRot), np.sin(self.aRot)
        plane_data = list(self.xyzLE) + [0., -s, c, 1., 0., 0.]
        plane = context.makeGeometry(egads.SURFACE, egads.PLANE, plane_data)
        airfoil = context.makeTopology(egads.FACE, egads.SFORWARD, geom=plane, children=[loop])
        return airfoil, nodes

def is_collocate(sec0, sec1):
    return sec0.xyzLE == sec1.xyzLE

def mean_angle(angles):
    return np.arctan2(np.sin(angles).mean(), np.cos(angles).mean())

def averageSections(s0, s1):
    aRot = mean_angle([s0.aRot, s1.aRot])
    thicken = 1 / np.cos(aRot - s0.aRot)
    airfoil = s0.airfoil * thicken
    return dataclasses.replace(s0,
        airfoil = airfoil,
        xyzLE = tuple((np.array(s0.xyzLE) + np.array(s1.xyzLE)) / 2),
        chord = (s0.chord + s1.chord) / 2,
        aInc = mean_angle([s0.aInc, s1.aInc]),
        aRot = aRot)

def node_pair_dist(nodes0, nodes1):
    def dist_ordered_sq(n0, n1):
        dist2 = [None,None]
        for i in (0,1):
            x0 = n0[i].getTopology()[3]
            x1 = n1[i].getTopology()[3]
            dist2[i] = ((np.array(x0) - np.array(x1))**2).sum()
        return max(dist2[0], dist2[1])
    dist2 = min(dist_ordered_sq(nodes0, nodes1),
                dist_ordered_sq(nodes0, [nodes1[1], nodes1[0]]))
    return np.sqrt(dist2)

def label_egads_edges(context, body, nodes, name):
    edges = body.getBodyTopos(egads.EDGE)
    for e in edges:
        ne = e.getTopology()[4]
        if len(ne) != 2:
            continue
        if node_pair_dist(ne, nodes) < 1E-7:
            e.attributeAdd('edgeName', name)
            return True
    return False

@dataclasses.dataclass
class Surface:
    '''
    An aerodynamic surface
    '''
    name : str
    sections : list[AirfoilSection]
    connected0 = None
    connected1 = None

    def estimate_span(self):
        span = None
        for sec in self.sections:
            if span is None:
                prev_sec = sec
                span = 0
            else:
                dxyz = np.array(sec.xyzLE) - np.array(prev_sec.xyzLE)
                span += np.sqrt((dxyz**2).sum())
        return span

    def to_egads(self, context, config, verbosity=0):
        return SurfaceChain([self], [False]).to_egads(context, config, verbosity)

def splitPeriodicSurface(s : Surface) -> tuple[Surface]:
    assert len(s.sections) >= 3
    iSplit = len(s.sections) // 2
    s0 = dataclasses.replace(s, sections=s.sections[:iSplit+1], name=s.name+'_L')
    s1 = dataclasses.replace(s, sections=s.sections[iSplit:], name=s.name+'_R')
    return (s0, s1)

def mirrorSurface(s : Surface, yDuplicate : float) -> Surface:
    sections = []
    for sec in s.sections:
        sections.append(dataclasses.replace(sec, aRot=-sec.aRot,
            xyzLE=(sec.xyzLE[0],
                   yDuplicate + (yDuplicate - sec.xyzLE[1]),
                   sec.xyzLE[2])))
    return dataclasses.replace(s,
            name = s.name + '_mirrored', sections = sections)

@dataclasses.dataclass
class SurfaceChain:
    '''
    A chain of connected aerodynamic surfaces
    '''
    surfaces : list[Surface]
    isReversed : list[bool]
    isPeriodic : bool = False

    @property
    def nSurfaces(self):
        return len(self.surfaces)

    @property
    def name(self):
        return '_'.join([s.name for s in self.surfaces])

    def first_section(self, iSurf=0):
        idx = -1 if self.isReversed[iSurf] else 0
        return self.surfaces[iSurf].sections[idx]

    def last_section(self, iSurf=-1):
        idx = 0 if self.isReversed[iSurf] else -1
        return self.surfaces[iSurf].sections[idx]

    def reversedCopy(self):
        idx = list(reversed(range(self.nSurfaces)))
        return dataclasses.replace(self,
                surfaces = [self.surfaces[i] for i in idx],
                isReversed=[not self.isReversed[i] for i in idx])

    def estimate_span(self):
        span = 0
        for s in self.surfaces:
            span += s.estimate_span()
        return span

    def iter_sections(self):
        for iSurf in range(self.nSurfaces):
            if iSurf == 0 and self.isPeriodic:
                prevSec = self.last_section()
            elif iSurf == 0:
                prevSec = self.first_section()
            else:
                prevSec = self.last_section(iSurf-1)

            if iSurf == self.nSurfaces - 1 and self.isPeriodic:
                nextSec = self.first_section()
            elif iSurf == self.nSurfaces - 1:
                nextSec = self.last_section()
            else:
                nextSec = self.first_section(iSurf+1)

            surfSecs = self.surfaces[iSurf].sections
            if self.isReversed[iSurf]:
                idx = reversed(range(len(surfSecs)))
                surfSecs = [surfSecs[i] for i in idx]
            yield [averageSections(prevSec, surfSecs[0]),
                   *surfSecs[1:-1],
                   averageSections(nextSec, surfSecs[-1])]

    def egads_sections(self, context, config):
        te_thickness = self.estimate_span() / config.max_span_tethickness_ratio
        for surface, sections in zip(self.surfaces, self.iter_sections()):
            face_list, nodes_list = [], []
            for sec in sections:
                face, nodes = sec.to_egads(context, te_thickness, verbosity=0)
                face_list.append(face)
                assert len(nodes) == 3
                nodes_list.append(nodes)
            yield surface.name, face_list, [[nodes_list[0][i], nodes_list[-1][i]] for i in range(3)]

    def to_egads(self, context, config, verbosity=0):
        names, bodies = [], []
        for name, faces, (te1_nodes, le_nodes, te2_nodes) in self.egads_sections(context, config):
            if self.isPeriodic or config.wing_tip == 0:
                rc0 = rc1 = None
            else:
                rc0 = [0, config.wing_tip] if name == self.surfaces[0].name else None
                rc1 = [0, config.wing_tip] if name == self.surfaces[-1].name else None
            body = egads.blend(faces, rc0, rc1)
            assert label_egads_edges(context, body, te1_nodes, 'trailingEdge')
            assert label_egads_edges(context, body, te2_nodes, 'trailingEdge')
            assert label_egads_edges(context, body, le_nodes, 'leadingEdge')
            for i in [0,1]:
                # assert label_egads_edges(context, body, [te1_nodes[i], te2_nodes[i]], 'surfaceJoint')
                assert label_egads_edges(context, body, [te1_nodes[i], le_nodes[i]], 'surfaceJoint')
                assert label_egads_edges(context, body, [te2_nodes[i], le_nodes[i]], 'surfaceJoint')
            for face in body.getBodyTopos(egads.FACE):
                face.attributeAdd('faceName', name)
            names.append(name)
            bodies.append(body.copyObject())
        if len(bodies) == 1:
            return bodies[0]
        bodyFaces = [b.getBodyTopos(egads.FACE) for b in bodies]
        for i0 in range(len(bodies)):
            i1 = (i0 + 1) % len(bodies)
            matches = bodies[i0].matchBodyFaces(bodies[i1])
            if verbosity:
                print(f'matching face {[m[0]-1 for m in matches]} from {names[i0]} and face {[m[1]-1 for m in matches]} from {names[i1]}')
            for j0, j1 in matches:
                bodyFaces[i0][j0-1] = None
                bodyFaces[i1][j1-1] = None
        all_faces = []
        for faces in bodyFaces:
            all_faces.extend([f for f in faces if f is not None])
        model = egads.sewFaces(all_faces)
        model.saveModel('tmp.egads', True)
        sewnBodies = model.getTopology()[4]
        assert len(sewnBodies) == 1, f'Sewing {self.name} got {len(sewnBodies)} bodies'
        assert sewnBodies[0].getTopology()[1] == egads.SOLIDBODY
        return sewnBodies[0].copyObject()

def tryConnect01(c0 : SurfaceChain, c1 : SurfaceChain):
    if is_collocate(c0.last_section(), c1.first_section()):
        assert not(c0.isPeriodic or c1.isPeriodic)
        return SurfaceChain(
                surfaces = c0.surfaces + c1.surfaces,
                isReversed = c0.isReversed + c1.isReversed,
                isPeriodic = is_collocate(c0.first_section(), c1.last_section()))
    return None

def tryConnect(c0, c1):
    if isinstance(c0, Surface):
        c0 = SurfaceChain([c0], [True])
    if isinstance(c1, Surface):
        c1 = SurfaceChain([c1], [True])

    conn = tryConnect01(c0, c1)
    if conn is not None: return conn
    conn = tryConnect01(c0, c1.reversedCopy())
    if conn is not None: return conn
    conn = tryConnect01(c0.reversedCopy(), c1)
    if conn is not None: return conn
    conn = tryConnect01(c0.reversedCopy(), c1.reversedCopy())
    if conn is not None: return conn
    return None

def tryConnectSurfaces(surfaces):
    while True:
        nSurf = len(surfaces)
        for i, j in itertools.combinations(range(nSurf), 2):
            s0, s1 = surfaces[i], surfaces[j]
            connected = tryConnect(s0, s1)
            if connected is not None:
                success = True
                print(f'{s0.name} and {s1.name} are connected')
                surfaces = [connected] + surfaces[:i] \
                        + surfaces[i+1:j] + surfaces[j+1:]
                break
        if len(surfaces) == 1 or len(surfaces) == nSurf:
            return surfaces

@dataclasses.dataclass
class FuselageSection:
    '''
    An fuselate section
    '''
    xyzCtr : tuple[float]
    width : float
    height : float
    shape : float

    def xyz(self):
        theta = np.linspace(0, 1, 4 * 12 + 1) * np.pi * 2
        y, z = np.sin(theta), -np.cos(theta)
        y = self.width/2 * np.power(np.abs(y), 2./self.shape) * np.sign(y)
        z = -self.height/2 * np.power(np.abs(z), 2./self.shape) * np.sign(z)
        return np.transpose([y * 0, y, z]) + np.array(self.xyzCtr)

    def to_egads(self, context, config):
        if self.width == 0 or self.height == 0:
            return context.makeTopology(egads.NODE, reals=list(self.xyzCtr))
        xyz = self.xyz()
        bspline = context.approximate([len(xyz), 0], xyz)
        tTop = bspline.invEvaluate(xyz[len(xyz)//2])[0]
        points = np.array([bspline.evaluate(t)[0] for t in [0, tTop]])
        nodes = [context.makeTopology(egads.NODE, reals=list(p)) for p in points]
        starboard = context.makeTopology(egads.EDGE, egads.TWONODE, geom=bspline,
                                         children=[nodes[0], nodes[1]], reals=[0,tTop])
        port = context.makeTopology(egads.EDGE, egads.TWONODE, geom=bspline,
                                    children=[nodes[1], nodes[0]], reals=[tTop,1])
        loop = context.makeTopology(egads.LOOP, egads.CLOSED, children=[starboard, port],
                                    senses=[egads.SFORWARD, egads.SFORWARD])
        plane_data = list(self.xyzCtr) + [0., 1., 0., 0., 0., 1.]
        plane = context.makeGeometry(egads.SURFACE, egads.PLANE, plane_data)
        xsec = context.makeTopology(egads.FACE, egads.SFORWARD, geom=plane, children=[loop])
        return xsec

@dataclasses.dataclass
class Fuselage:
    '''
    A fuselage with streamwise cross sections
    '''
    name : str
    sections : list[FuselageSection]

    def egads_sections(self, context, config):
        face_list = []
        for sec in self.sections:
            face_list.append(sec.to_egads(context, config))
        return face_list

    def to_egads(self, context, config):
        faces = self.egads_sections(context, config)
        body = egads.blend(faces)
        for face in body.getBodyTopos(egads.FACE):
            face.attributeAdd('faceName', self.name)
        return body

@dataclasses.dataclass
class CenterBody:
    '''
    Axisymmetric center body
    '''
    name : str
    profile : np.ndarray
    translate : tuple[float]
    angle : float
    scale : float

    def to_egads(self, context, config):
        c, s = np.cos(self.angle), np.sin(self.angle)
        mat = (( c*self.scale, 0         , +s*self.scale, self.translate[0]),
               ( 0           , self.scale,  0           , self.translate[1]),
               (-s*self.scale, 0         ,  c*self.scale, self.translate[2]))
        transform = context.makeTransform(mat)

        length = self.profile[:,0].max() - self.profile[:,0].min()
        te_thickness = length / config.max_body_butt_ratio
        xz = thicken_te(self.profile, te_thickness)
        xyz = np.transpose([xz[:,0], xz[:,0] * 0, xz[:,1]])
        bspline = context.approximate([len(xyz), 0], xyz)
        zCenter = (xz[0,1] + xz[-1,1]) / 2
        tLE = bspline.invEvaluate([xz[:,0].min() - length * 100, 0, zCenter])[0]
        points = np.array([[xz[:,0].max(), 0, zCenter],
                           bspline.evaluate(0)[0],
                           [xz[:,0].min(), 0, zCenter]])
        nodes = [context.makeTopology(egads.NODE, reals=list(p)) for p in points]
        upper = context.makeTopology(egads.EDGE, egads.TWONODE, geom=bspline,
                                     children=[nodes[1], nodes[2]], reals=[0,tLE])
        te = egads_line(context, nodes[0], nodes[1])

        loop = context.makeTopology(egads.LOOP, egads.OPEN,
                                    children=[upper, te],
                                    senses=[egads.SREVERSE]*2)
        rot180 = context.makeTransform(((1, 0, 0,0),
                                        (0,-1, 0,0),
                                        (0, 0,-1,2 * zCenter)))
        loops = [loop, loop.copyObject(rot180)]
        sheet_bodies = [l.rotate(-180, [(0,0,zCenter), (1,0,0)]) for l in loops]
        # return [b.copyObject(transform) for b in sheet_bodies]
        fuse = sheet_bodies[0].fuseSheets(sheet_bodies[1])
        body = context.makeTopology(egads.BODY, egads.SOLIDBODY,
                                    children=fuse.getBodyTopos(egads.SHELL))
        for face in body.getBodyTopos(egads.FACE):
            face.attributeAdd('faceName', self.name)
        return body.copyObject(transform)

def mirrorCenterBody(b : CenterBody, yDuplicate : float) -> CenterBody:
    translateY = yDuplicate + (yDuplicate - b.translate[1])
    return dataclasses.replace(b,
            name = b.name + '_mirrored',
            translate = (b.translate[0], translateY, b.translate[2]))

def unionBodies(bodies, verbosity):
    while len(bodies) > 1:
        n_bodies = len(bodies)
        names = list(bodies.keys())
        random.shuffle(names)
        if verbosity:
            print(f'{n_bodies} bodies remaining')
        for i in range(n_bodies):
            for j in range(i):
                if verbosity:
                    print(f'Trying {names[i]} and {names[j]}')
                model = bodies[names[j]].copyObject().generalBoolean(bodies[names[i]].copyObject(), egads.FUSION)
                new_bodies = model.getTopology()[4]
                if len(new_bodies) == 1:
                    new_name = names[i] + '_' + names[j]
                    print(f'Union {names[i]} and {names[j]} suceeds')
                    new_bodies = {new_name : new_bodies[0].copyObject()}
                    for k in range(n_bodies):
                        if k != i and k != j:
                            new_bodies[names[k]] = bodies[names[k]]
                    bodies = new_bodies
                    break
            else:
                continue
            break
        else:
            break
    return bodies

