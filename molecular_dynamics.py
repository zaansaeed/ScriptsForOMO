# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 16:26:23 2021


@author: wanxuan
"""

"""Explicit Function"""

from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
# import mdtraj as md
import random
import time
from openmm.unit import nanoseconds, picoseconds, kelvin, picosecond, nanometer


def explicit_function(PDB_input_file, final_pdb_file, DCD_output_file,
                      equilibration_time=0.01 * nanoseconds,
                      simulation_time=10 * nanoseconds,
                      t_step=0.001 * picoseconds,
                      t_report=10 * picoseconds,
                      temp=300, pH=7.0):
    print("Uploading Model...")
    pdb = PDBFile(PDB_input_file)
    forcefield = ForceField('UNK_DD9ED3.xml', 'tip3p.xml')  # TIP3P model

    print("Building Model...")
    modeller = Modeller(pdb.topology, pdb.positions)

    print('Adding hydrogens...')
    modeller.addHydrogens(forcefield, pH=pH)  # Adding missing hydrogen ends

    print("Adding Water Molecules...")
    modeller.addSolvent(forcefield, model='tip3p',
                        padding=1.0 * nanometer)  # Adding water molecules surrounding the protein strain

    print("Creating System...")
    random.seed(round(time.time() * 1000))  # randomize starting point of random numbers
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0 * nanometer,
                                     constraints=HBonds)
    integrator = LangevinIntegrator(temp * kelvin, 1 / picosecond, 0.001 * picoseconds)  # step size: 0.001*picoseconds

    print("Creating Simulation Context...")
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    # Minimizing System
    print("Minimizing Energy...")
    simulation.minimizeEnergy(maxIterations=25)

    # Equilibrate the System
    print("Equilibrating...")
    # simulation.step(10000000) # equilibrate for 10 ns (use on the cluster)
    simulation.step(int(equilibration_time / t_step))  #

    # Adding Reporters
    print("Collecting Data...")
    Nsteps_report = int(t_report / t_step)
    simulation.reporters.append(DCDReporter(DCD_output_file, Nsteps_report))  # Every 100ps record one traj data
    simulation.reporters.append(
        StateDataReporter(stdout, Nsteps_report, step=True, potentialEnergy=True, temperature=True, separator=' '))
    simulation.reporters.append(PDBReporter(final_pdb_file, Nsteps_report))

    # Running simulation
    # simulation.step(100000000)
    simulation.step(int(simulation_time / t_step))


def gyration(dcd_file, pdb_file):
    import mdtraj as md
    print("Calculating gyration radius...")
    t = md.load_dcd(dcd_file, top=pdb_file)
    gy = md.compute_rg(t)
    print(gy)  # print gyration radius of protein strain
    print(t)


if __name__ == '__main__':  # explicit_function(PDB_input_file,final_pdb_file,DCD_output_file,temp)
    geo_list = "UNK_DD9ED3"
    start_pdb = "UNK_DD9ED3.pdb"
    final_pdb = "UNK_DD9ED3_final.pdb"
    traj_dcd = geo_list + "_traj.dcd"

    # amino_to_pdb(geo_list, start_pdb, box_scaling = 3)
    #explicit_function(start_pdb, final_pdb, traj_dcd, temp=300)
    print("simulation complete")

    # gyration('Explicit_final_25_100ns.dcd','Explicit_final_20_500.pdb') #gyration(DCD_input_file,PBD_input_file)

    import mdtraj as md
    import matplotlib.pyplot as plt
    import numpy as np
    import MDAnalysis as mda
    from MDAnalysis.coordinates.DCD import DCDWriter

    traj = md.load("trajectory_no_water.dcd", top="structure_no_water.pdb")

    atom_indices = [a.index for a in traj.topology.atoms if a.element.symbol != "H"]
    distances = np.empty((traj.n_frames, traj.n_frames))
    for i in range(traj.n_frames):
        distances[i] = md.rmsd(traj, traj, i, atom_indices=atom_indices)
    beta = 1
    index = np.exp(-beta * distances / distances.std()).sum(axis=1).argmax()
    centroid = traj[100]
    centroid.save_pdb("centroid.pdb")





