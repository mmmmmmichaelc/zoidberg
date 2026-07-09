""" 
Generates a flux-coordinate independent (FCI) grid file (.nc) to be used with an FCI-supporting
branch of Hermes-3. Note that this code has been tested on only two-field-period stellarators
and may therefore need adjustment in your case. Run this script with the --help tag to see
more information about the necessary inputs.
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import dill
import argparse
import os
from zoidberg.field import MGRID
from zoidberg.grid import Grid
from zoidberg.rzline import line_from_points
from zoidberg.fieldtracer import trace_poincare
from zoidberg.poloidal_grid import grid_elliptic
from zoidberg.zoidberg import make_maps, write_maps


parser = argparse.ArgumentParser()
parser.add_argument('--nx', type=int, default=64,
    help="Grid resolution in the x (radial) direction. Default: 64")
parser.add_argument('--ny', type=int, default=64,
    help="Grid resolution in the y (toroidal) direction. Default: 64")
parser.add_argument('--nz', type=int, default=64,
    help="Grid resolution in the z (poloidal) direction. Default: 64")
parser.add_argument('--vmgrid', type=str, required=True,
    help="Stellopt .mgrid file specifying the vacuum magnetic field")
parser.add_argument('--pmgrid', type=str, required=True,
    help="Stellopt .mgrid file specifying the finite beta magnetic field.\n" \
    "That is, the vmgrid and pmgrid fields should add to the actual field.")
parser.add_argument('--vessel', type=str, required=True,
    help="File specifying the outer boundary. Probably the output of \n" \
    "read_vessel.py.")
parser.add_argument('--start_r', type=float, required=True,
    help="Location in R (cylindrical coordinates) to begin field line tracing\n" \
    "to determine the geometry of the inner boundary. Tracing begins at y=0.")
parser.add_argument('--start_z', type=float, required=True,
    help="Location in Z (cylindrical coordinates) to begin field line tracing\n" \
    "to determine the geometry of the inner boundary. Tracing begins at y=0.")
parser.add_argument('--gridname', type=str, default="fci_grid.nc",
    help="A name for the output grid file. Default: fci_grid.nc. Note this is stored \n"\
         "in the same directory as the vessel file (do not specify a path).")
parser.add_argument('--plot', action=argparse.BooleanOptionalAction,
    help="Shows each of the grids using plt.show() as they're generated. Also saves\n" \
         "visualizations of the grid to the same folder that constains the vessel file.\n" \
         "Default: False")
parser.set_defaults(feature=False)

args = parser.parse_args()

plotting = args.plot

vacuum_mgrid = args.vmgrid
plasma_mgrid = args.pmgrid
vessel_file = args.vessel
file_dir = "/".join(vessel_file.split("/")[:-1]) + "/"

assert os.path.exists(file_dir), "Issue interpreting vessel file directory."
assert os.path.exists(vacuum_mgrid), "Invalid file path specified for input vmgrid"
assert os.path.exists(plasma_mgrid), "Invalid file path specified for input pmgrid"
assert os.path.exists(vessel_file), "Invalid file path specified for input vessel"

magnetic_field = MGRID(vacuum_mgrid, plasma_mgrid)

with open(vessel_file, "rb") as f:
    phi_coords, outer_lines = dill.load(f)

num_original_phi = len(phi_coords)

start_r = args.start_r
start_z = args.start_z
start_y = 0.0

nslices = args.ny
assert len(outer_lines) % nslices == 0, f"An ny which does not divide the toroidal resolution {len(outer_lines)}" \
                                        "of the vessel file is unsupported"

ycoords = np.array([phi_coord for i, phi_coord in enumerate(phi_coords) if ((i % (len(phi_coords) // nslices)) == 0)])

outer_lines = [outer_line for i, outer_line in enumerate(outer_lines) if ((i % (num_original_phi // nslices)) == 0)]

nx = args.nx
nz = args.nz

print("Tracing particles...")
coord, y_slices = trace_poincare(magnetic_field, start_r, start_z, y_slices=ycoords, yperiod=2*np.pi) # TODO hardcoded for Eos, Helios

inner_lines = []

for i in range(len(y_slices)):
    r = coord[:, i, 0, 0]
    z = coord[:, i, 0, 1]
    line = line_from_points(r, z, smooth=True)
    line = line.equallySpaced(n=nz)
    inner_lines.append(line)


inner_lines = inner_lines[:len(outer_lines)]

print("Generating grid...")

pol_slices = [
    grid_elliptic(inner_line, outer_line, nx, nz, show=plotting, tol=1e-5)
    for inner_line, outer_line in zip(inner_lines, outer_lines)
]

grid = Grid(pol_slices, ycoords, np.pi, yperiodic=True) # TODO hardcoded for Eos, Helios

if plotting:
    print("Plotting...")
    fig = plt.figure()
    ax = fig.add_subplot()
    y_index = 0
    poloidal_grid, y_val = grid.getPoloidalGrid(y_index)
    poloidal_grid.plot(axis=ax)
    fig.savefig(file_dir + f"y_{y_index}_" + args.gridname.replace(".nc", ".pdf"))

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    grid.plot3D(axis=ax)
    fig.savefig(file_dir + args.gridname.replace(".nc", ".pdf"))

print("Making maps...")
maps = make_maps(grid, magnetic_field)

print("Writing maps...")
write_maps(grid, magnetic_field, maps, file_dir + args.gridname, metric2d=False)

print(f"Done writing {file_dir + args.gridname}.")

print("Making a Jacobian free copy...")
src = nc.Dataset(f"{file_dir + args.gridname}", "r")
dst = nc.Dataset(f"{(file_dir + args.gridname).replace(".nc", "_noJ.nc")}", "w", format=src.file_format)

for name, dim in src.dimensions.items():
    dst.createDimension(name, (len(dim) if not dim.isunlimited() else None))

for name, var in src.variables.items():
    if name == "J":
        continue
    out_var = dst.createVariable(name, var.datatype, var.dimensions)
    out_var.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
    out_var[:] = var[:]

dst.setncatts({k: src.getncattr(k) for k in src.ncattrs()})

src.close()
dst.close()
print(f"Done writing {(file_dir + args.gridname).replace(".nc", "_noJ.nc")}")
print("The '_noJ' version is the one compatible with BSTING.")

