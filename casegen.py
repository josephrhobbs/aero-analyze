# Generates Flow360 Configuration Cases

import argparse
from graph import AvlInput
from math import sqrt

GAMMA = 1.4
R = 287.052874

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

def compute_reference_dimensions(avl_fname):
    avlInput = AvlInput(avl_fname)
    wing = avlInput.surfaces[0]
    cref = avlInput.Cref
    bref = 2 * wing.sections[-1].xyzLE[1]
    return bref, cref

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='casegen',
            description='generate Flow360 configuration files')
    parser.add_argument('avl_fname')
    parser.add_argument('flow360_fname')
    args = parser.parse_args()

    bref, cref = compute_reference_dimensions(args.avl_fname)
    Sref = bref * cref

    mach, aoa, temp, kinvisc, rampsteps = process_conditions(args.flow360_fname)
    uinf = mach * sqrt(GAMMA * R * temp)
    re = uinf * cref / kinvisc

    DEFAULT = f'''
{{
	"geometry": {{
		"comments": {{
			"meshUnit": "m"
		}},
		"refArea": {Sref},
		"momentCenter": [
			0,
			0,
			0
		],
		"momentLength": [
			{bref},
			{cref},
			{bref}
		]
	}},
	"freestream": {{
		"Reynolds": {re},
		"Mach": {mach},
		"Temperature": {temp},
		"alphaAngle": {aoa},
		"betaAngle": 0
	}},
	"volumeZones": {{
		"stationaryBlock": {{
			"modelType": "FluidDynamics"
		}}
	}},
	"boundaries": {{
		"stationaryBlock/farfield": {{
			"type": "Freestream"
		}},
		"stationaryBlock/unspecified": {{
			"type": "NoSlipWall"
		}}
	}},
	"actuatorDisks": [],
	"BETDisks": [],
	"timeStepping": {{
		"maxPseudoSteps": 5000,
		"CFL": {{
			"type": "ramp",
			"initial": 1,
			"final": 100,
			"rampSteps": {rampsteps}
		}},
		"physicalSteps": 1,
		"timeStepSize": "inf",
		"comments": {{}}
	}},
	"navierStokesSolver": {{
		"absoluteTolerance": 1e-9,
		"linearIterations": 35,
		"kappaMUSCL": -1,
		"orderOfAccuracy": 2
	}},
	"turbulenceModelSolver": {{
		"modelType": "SpalartAllmaras",
		"rotationCorrection": false,
		"absoluteTolerance": 1e-8,
		"linearIterations": 25,
		"kappaMUSCL": -1,
		"orderOfAccuracy": 2
	}},
	"heatEquationSolver": {{}},
	"volumeOutput": {{
		"outputFormat": "tecplot",
		"primitiveVars": false,
		"vorticity": false,
		"Cp": true,
		"Mach": true,
		"qcriterion": true
	}},
	"surfaceOutput": {{
		"outputFormat": "tecplot",
		"Cp": true,
		"Cf": false,
		"CfVec": true
	}},
	"aeroacousticOutput": {{
		"observers": []
	}}
}}
'''

    print(DEFAULT)