"""
Find how far are the waters

"""

# python ~/path/scripts/mda/find_3_waters.py -nwater 1 -top morph_solv.prmtop -traj prod.dcd &

from pathlib import Path

import MDAnalysis as mda

from MDAnalysis.analysis.distances import distance_array

import numpy as np

import time

import argparse


def find_waters(topo, traj, nwater, solvated_complex, prot_sel='not resname LIG and not resname WAT', lig_sel='resname LIG', wat_sel='resname WAT',

                output_file='closest_water_ids.dat', output_summary_file='cumulative_water_time.dat',

                traj_dcd_prefix='com_wat'):

    start_time = time.time()

    LIG = []

    with open(solvated_complex,'r') as f:
        for lines in f.readlines():
            if 'LIG' in lines:
                LIG.append(lines)
    peptide_ligand = (len(LIG) == 0)

    if peptide_ligand:
        OXT_lines = []
        with open(solvated_complex,'r') as f:
            lines = f.readlines()
            filtered_lines = []
            prev_line = None

            for idx,line in enumerate(lines):
                if 'TER' in line:
                    if line.strip() != prev_line.strip():
                        prev_line = line
                        filtered_lines.append(line)
                else:
                    prev_line = line
                    filtered_lines.append(line)

            for line_no,line in enumerate(filtered_lines):
                if 'OXT' in line and 'LIG' not in line:
                    OXT_lines.append(line_no)
        resid_start = filtered_lines[OXT_lines[-2]+2].split()[4]
        resid_end = filtered_lines[OXT_lines[-1]].split()[4]

        lig_sel = "resid {}:{}".format(resid_start,resid_end)
        prot_sel = "not {} and not resname WAT".format(lig_sel)

    print(lig_sel)
    print(prot_sel)


    u = mda.Universe(topo, traj)

    print('Loaded', u.trajectory)

      

    prot = u.select_atoms(prot_sel)

    print('Prot: ', len(prot))

    lig = u.select_atoms(lig_sel)

    print('Lig', len(lig))

    water_o = u.select_atoms(wat_sel)

    print('Wat', len(water_o))


    closest_wat_number = nwater

    # create the trajectory writers, which will

    # mda_writers = [mda.Writer(f'{traj_dcd_prefix}{watNum}.dcd', n_atoms=len(lig) + len(prot) + 3*watNum)

    #                for watNum in range(0, closest_wat_number + 1)]

    out_traj = mda.Writer(f'com_wat%d.dcd'% nwater, n_atoms=len(lig) + len(prot) + 3 * closest_wat_number)



    cum_sums = np.zeros(len(water_o.residues))



    every_nth_frame = 1

    number_of_frames = 3

    start_frame = 0 # 1 means the second frame

    water_ids_per_frame = []

    for ts in u.trajectory[start_frame::every_nth_frame]:

        print("Time", ts.time)

        # get distances from water to X

        wat_prot = distance_array(water_o.positions, prot.positions, box=ts.dimensions)

        wat_lig = distance_array(water_o.positions, lig.positions, box=ts.dimensions)


        # take the shortest distances

        wat_prot_shortest = np.min(wat_prot, axis=1)

        wat_lig_shortest = np.min(wat_lig, axis=1)


        # get the minimum out of the 2 hydrogens and the oxygen

        wat_prot_res_min = np.min(wat_prot_shortest.reshape([len(water_o.residues), 3]), axis=1)

        wat_lig_res_min = np.min(wat_lig_shortest.reshape([len(water_o.residues), 3]), axis=1)


        # add up the shortest distances

        wat_others_sum = wat_prot_res_min + wat_lig_res_min

        cum_sums += wat_others_sum



        # PER FRAME SAVE

        # prepare for sorting

        chained = [(wid, dst) for wid, dst in zip(water_o.residues.resids, wat_others_sum)]

        # order by the distance

        chained_ord = sorted(chained, key=lambda x: x[1])

        # extract only the water IDs

        water_ids_per_frame.append([ts.frame] + [wid for wid, dst in chained_ord])



        # select the best waters, get their IDs and convert them to indices

        com_wat_resids = [chained_ord[x][0] - 1 for x in range(closest_wat_number)]

        com_wat = u.residues[com_wat_resids]

        #print('Writing atoms to the traj:', len(prot + lig + com_wat.atoms))

        out_traj.write(prot + lig + com_wat.atoms)



    # save the water molecules

    #np.savetxt(output_file, water_ids_per_frame, fmt='%d', header='Frame closest_watID1 closest_watID2 ..')



    # bind the two fields together for sorting

    chained = [(wid, cum) for wid, cum in zip(water_o.residues.resids, cum_sums)]

    # order by the cumulative sum

    chained_ord = sorted(chained, key=lambda x:x[1])

    to_save = np.array(chained_ord).T

    print("Final shape", to_save.shape)

    # sort according to the cumulative time

    #np.savetxt(output_summary_file, to_save, fmt='%.1f')

    print(f'Altogether the analysis took {time.time() - start_time:.0f} seconds')





if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generating dcd with 3 closest water molecules')

    parser.add_argument('-traj', type=str, help='Trajectory file', dest='traj')

    parser.add_argument('-top', type=str, help='Topology file', dest='top')
    
    parser.add_argument('-nwater', type=int, help='number of water molecules', dest='nwater')

    parser.add_argument('-solvated_complex', type=str, help='complex solvated pdb file', dest='solvated_complex',default=None)

    args = parser.parse_args()



    find_waters(args.top, args.traj, args.nwater,args.solvated_complex)



    # example:

    # /home/dresio/ucl/validation/resp_validation/mcl1/l1_l8/complex/lambda_0.00/rep1

    # python ~/path/scripts/mda/find_3_waters.py -top morph_solv.prmtop -traj prod.dcd &

