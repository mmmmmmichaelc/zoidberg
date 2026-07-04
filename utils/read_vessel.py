"""
reads the vessel file (.part) from Erik
and turns it into something useful for zoidberg
(set of rzlines or a boundary object or something)
"""

import numpy as np
from zoidberg.rzline import line_from_points
import dill
from tqdm import tqdm
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("--vessel", type=str, required=True,
    help=r"Vessel file (possibly with extension .part). Format:\n" \
         r"\tlabel\n" \
         r"\tn  m  symmetry  {R0  Z0}\n" \
         r"\tfor i = 0...n-1\n" \
             r"\t\tphi_i [deg]\n" \
             r"\t\tfor j = 0...m-1\n" \
                 r"\t\t\tR_ij   Z_ij")
parser.add_argument("--name", type=str, default="vessel.dill",
    help="Name for the output (.dill) file. Default: vessel.dill Note \nthat" \
         "this should not include the file path. It will be stored in the same \n"\
          "directory as the input vessel file.")
parser.add_argument("--nfp", type=int, default=1,
    help="Number of field periods in the stellarator. Will truncate the\n" \
         "vessel file output to range from phi=0 to phi=2pi/nfp. Default: 1")

file = parser.parse_args().vessel
assert os.path.exists(file), "Invalid file path specified for vessel."

with open(file, "r") as f:
    lines =f.readlines()

file_dir = "/".join(file.split("/")[:-1]) + "/"

assert os.path.exists(file_dir), "Issue interpreting file directory."

try:
    label = lines[0]

    num_phi, num_points, _, __, ___ = lines[1].replace(" \n", "").split(" ")
    num_phi = int(num_phi)
    num_points = int(num_points)

    idx = 1
    cutoff = 360/parser.parse_args().nfp
    phi_slices = []
    phi_coords = []
    for i in range(num_phi):
        idx += 1
        phi_line = lines[idx].replace(" \n", "").split(" ")
        assert len(phi_line) == 1, "mismatch"
        if float(phi_line[0]) < cutoff and not np.isclose(float(phi_line[0]), cutoff, rtol=1e-3):
            phi_coords.append(float(phi_line[0]))

            points = []
            for j in range(num_points):
                idx += 1
                R_ij, Z_ij = lines[idx].replace(" \n", "").replace(",", "").split(" ")
                points.append([float(R_ij), float(Z_ij)])
                
            phi_slices.append(points)
        else:
            break

    phi_coords = np.array(phi_coords)*np.pi/180 # to rad
except:
    print("Encountered an exception when parsing the file. Ensure \n" \
    "that your file given to --vessel has the following format:" \
    "\tlabel\n" \
    "\tn  m  symmetry  {R0  Z0}\n" \
    "\tfor i = 0...n-1\n" \
        "\t\tphi_i [deg]\n" \
        "\t\tfor j = 0...m-1\n" \
            "\t\t\tR_ij   Z_ij")
    raise ValueError("Ensure your file format matches the style printed above.\n" \
                    r"Watch out for invisible text such as '\n'.")

rzlines = []
for phi_slice in tqdm(phi_slices, desc="Crafting RZlines..."):
    rs, zs = np.array(phi_slice).T
    rzlines.append(line_from_points(rs, zs))

output = phi_coords, rzlines

# print(f"Vessel phi (y) range: {phi_coords[0]} to {phi_coords[-1]}")

with open(file_dir + parser.parse_args().name, "wb") as f:
    dill.dump(output, f)

print("Done writing " + file_dir + parser.parse_args().name + ".")