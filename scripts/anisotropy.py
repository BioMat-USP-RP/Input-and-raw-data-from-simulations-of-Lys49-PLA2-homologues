import MDAnalysis as mda
import numpy as np
import pandas as pd

# Script used for analyis of geometric parameters later used for calculation of FRET efficience
# Orientation factor
# Trp distance
# Angles theta t, d and a, as described in the paper

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
    #angle_degrees = np.degrees(angle_radians)

    return angle_radians

def transition_dipole(trp):
    #Defines the transitions dipoles as rotations of the long axis vector
    #about the normal vector of the indol plane with origin in the COG of the indol group

    indol = trp.select_atoms("name NE1 CD1 CG CD2 CE2 CE3 CZ3 CH2 CZ2")
    COG = indol.center_of_geometry()

    #The vector goes from cd1 to the COG
    #Translate to origin

    cd1 = trp.select_atoms("name CD1")
    long_axis = cd1.positions-COG

    #Define normal do the indol group
    #First find two vectors in the indol plane

    cz3 = trp.select_atoms("name CZ3")
    ch2 = trp.select_atoms("name CH2")
    plane_vector1 = cz3.positions - COG
    plane_vector2 = ch2.positions - COG

    #take cross product and normalize

    normal = np.cross(plane_vector1, plane_vector2)
    normal = normal/np.linalg.norm(normal)

    #Rotate the long axis vector to give the two trasition dipoles using the Rodrigues formula

    angle1 = np.deg2rad(-38)
    angle2 = np.deg2rad(54)

    La = long_axis*np.cos(angle1) \
      + np.cross(normal, long_axis)*np.sin(angle1) \
      + normal*(np.inner(normal, long_axis))*(1-np.cos(angle1))

    Lb = long_axis*np.cos(angle2) \
      + np.cross(normal, long_axis)*np.sin(angle2) \
      + normal*(np.inner(normal, long_axis))*(1-np.cos(angle2))

    return La, Lb

def separation_vector(trp1, trp2):
    #calculates de separation vector between the two indol groups

    indol1 = trp1.select_atoms("name NE1 CD1 CG CD2 CE2 CE3 CZ3 CH2 CZ2")
    indol2 = trp2.select_atoms("name NE1 CD1 CG CD2 CE2 CE3 CZ3 CH2 CZ2")

    sep_vector = indol2.center_of_geometry() - indol1.center_of_geometry()

    return sep_vector

#Load trajectory and topology

topology_file = '../input/conf.pdb'
trajectory_file = 'center.xtc'
u = mda.Universe(topology_file, trajectory_file)

TRP_donnor = u.select_atoms("chainID A and resname TRP")
TRP_acceptor = u.select_atoms("chainID B and resname TRP")

timeseries = []
# vectors_donnor = transition_dipole(TRP_donnor)

for ts in u.trajectory:

  vectors_donnor = transition_dipole(TRP_donnor)
  vectors_acceptor = transition_dipole(TRP_acceptor)

  # Angle between trasition moments
  angle_trans = calculate_angle(vectors_donnor[0], vectors_acceptor[0])

  # Angles between trasition moments and separation vector
  sep_vec = separation_vector(TRP_donnor, TRP_acceptor)
  angle_Rd = calculate_angle((vectors_donnor[0]), sep_vec)
  angle_Ra = calculate_angle((vectors_acceptor[0]), sep_vec)

  # Orientation factor
  k = np.cos(angle_trans) - (3 * np.cos(angle_Rd) * np.cos(angle_Ra))

  # Depolarization
  depol = 0.18*((1.5*np.power(np.cos(angle_trans),2)) - 0.5)
  # depol = 0.18 - depol

  timeseries.append([ts.frame, k[0][0]*k[0][0], np.linalg.norm(sep_vec), depol[0][0], angle_trans[0][0], angle_Rd[0], angle_Ra[0]])

  data = np.array(timeseries, dtype="object")

data_df = pd.DataFrame(data, columns=['Frame','Orientation factor', 'COG distance', 'Depolarization', 'angle_trans', 'angle_Rd', 'angle_Ra'])

t = pd.DataFrame(data_df, columns=['COG distance', 'Orientation factor', 'angle_trans', 'angle_Rd', 'angle_Ra'])
t.to_csv('aniso-La-La.dat', index=False, sep=' ')


