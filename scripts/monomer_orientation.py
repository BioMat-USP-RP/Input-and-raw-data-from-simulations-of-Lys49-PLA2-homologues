import MDAnalysis as mda
import numpy as np
import argparse

# Script used for calculation of azimuthal and tilt angle, as described in the paper

def calculate_angle(vec1, vec2):

    dot_product = np.inner(vec1, vec2)

    # Calculate the magnitudes (norms) of the vectors
    norm_vec1 = np.linalg.norm(vec1)
    norm_vec2 = np.linalg.norm(vec2)

    # Calculate the cosine of the angle
    cos_angle = dot_product / (norm_vec1 * norm_vec2)

    # Calculate the angle in radians
    angle_radians = np.arccos(cos_angle)

    # Convert the angle to degrees
    # angle_degrees = np.degrees(angle_radians)

    return angle_radians


def get_vectors(protein, resolution='cg'):
    select_atomistic = {
    'A1': 'chainID A and resid 108 97 and name CA',
    'A2': 'chainID A and resid 108 43 and name CA',
    'B1': 'chainID B and resid 108 97 and name CA',
    'B2': 'chainID B and resid 108 43 and name CA',
    }
    select_cg = {
    'A1': 'chainID A and resid 98 87 and name BB',
    'A2': 'chainID A and resid 98 42 and name BB',
    'B1': 'chainID B and resid 98 87 and name BB',
    'B2': 'chainID B and resid 98 42 and name BB',
    }

    if resolution == 'atomistic': 
        selection = select_atomistic
    elif resolution == 'cg':
        selection = select_cg

    vectors = {}

    for k, v in selection.items():
        s = protein.select_atoms(v).positions
        s = s[1] - s[0]
        s = s / np.linalg.norm(s)
        vectors[k] = s 

    return vectors.values()


def print_progress(i, traj_size):
    progress_percent = 100 * i/traj_size
    print(f"Progress: {progress_percent:.1f}%", end='\r')


# Parse command line
parser = argparse.ArgumentParser()

parser.add_argument("-f",  help="Trajectory file", dest='trajectory', type=str)
parser.add_argument("-s", help="Topology file (Should be a .pdb or .tpr)", dest='topology', type=str)
parser.add_argument("-res", help="Resolution (Options: 'cg' or 'atomistics')", dest='resolution', type=str)
parser.add_argument("-o", help="Output file", dest='out_file', type=str)

args = parser.parse_args()


topology_file = args.topology  # f'{dir}conf.pdb'
trajectory_file = args.trajectory  # f'{dir}center.xtc'
u = mda.Universe(topology_file, trajectory_file)

protein = u.select_atoms('not resname DPPG DPPC W')


traj = u.trajectory[0:-1]
traj_size = len(traj)
time_series = np.zeros((2, traj_size))


for i, ts in enumerate(traj):
    A1, A2, B1, B2 = get_vectors(protein, args.resolution)

    angle_A = calculate_angle(B1, A1)
    angle_T = calculate_angle(np.cross(B2, B1), np.cross(A1, A2))

    time_series[0][i] = np.rad2deg(angle_A)
    time_series[1][i] = np.rad2deg(angle_T)

    print_progress(i, traj_size)


np.savetxt(
    args.out_file, 
    time_series.T, 
    header=f'Theta A\tTheta T', 
    delimiter=f'\t',
    fmt='%.5e'
)


