from openmm.app import *
from openmm import *                    
from openmm.unit import *
from openmmplumed import PlumedForce
import time
import os


t0 = time.time()

rst7 = AmberInpcrdFile("sim1.rst7")
prmtop = AmberPrmtopFile("complex_solvated.prmtop")
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)

system.addForce(MonteCarloBarostat(1.01325*bar, 300*kelvin, 25))

script = "plumed.dat"
with open(script, "r") as plumed_file:
    plumed = plumed_file.read()
system.addForce(PlumedForce(plumed))

platform = Platform.getPlatformByName('CUDA')
properties = {'Precision': 'mixed'}

integrator = LangevinIntegrator(300*kelvin, 5/picosecond, 0.002*picoseconds)


simulation = Simulation(prmtop.topology, system, integrator, platform, properties)


simulation.reporters.append(DCDReporter('plum.dcd', stride_size))

simulation.reporters.append(StateDataReporter('plum.log', stride_size, step=True, potentialEnergy=True, temperature=True))

simulation.context.setPositions(rst7.positions)
# Simulation

simulation.step(step_size)  

simulation.saveState('plum.xml')


t1 = time.time()
print('Simulation complete in', t1 - t0, 'seconds')

