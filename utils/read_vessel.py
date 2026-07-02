"""
reads the vessel file (.part) from Erik
and turns it into something useful for zoidberg
(set of rzlines or a boundary object or something)
"""

import numpy as np
from zoidberg.rzline import RZline, line_from_points
import dill
from tqdm import tqdm

file = "G1600-12-89-simplified_vessel.part"

with open(file, "r") as f:
    lines =f.readlines()

label = lines[0]

num_phi, num_points, _, __, ___ = lines[1].split(" ")
num_phi = int(num_phi)
num_points = int(num_points)

idx = 1
phi_slices = []
phi_coords = []
for i in range(num_phi):
    idx += 1
    phi_line = lines[idx].split(" ")
    assert len(phi_line) == 1, "mismatch"
    phi_coords.append(float(phi_line[0]))

    points = []
    for j in range(num_points):
        idx += 1
        R_ij, Z_ij = lines[idx].replace("\n", "").replace(",", "").split(" ")
        points.append([float(R_ij), float(Z_ij)])
        
    phi_slices.append(points)

phi_coords = np.array(phi_coords)*np.pi/180 # to rad

rzlines = []
for phi_slice in tqdm(phi_slices, desc="Crafting RZlines..."):
    rs, zs = np.array(phi_slice).T
    rzlines.append(line_from_points(rs, zs))

output = phi_coords, rzlines

output_file = "vessel_rzlines.dill"
with open(output_file, "wb") as f:
    dill.dump(output, f)

print("Done.")