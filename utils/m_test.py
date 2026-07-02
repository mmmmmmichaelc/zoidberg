from netCDF4 import Dataset
from zoidberg.plot import plot_forward_map, plot_backward_map, plot_poincare
import dill
import numpy as np
from zoidberg.rzline import RZline
from zoidberg.field import VMEC, MGRID

def dedill_config(filename):
    # michael addition

    with open(filename, 'rb') as f:
        grid, map, magnetic_field = dill.load(f)

    return grid, map, magnetic_field

filename = "rotating_ellipse.dill.dill"

grid, maps, field = dedill_config(filename)

y_index = 0
plot_forward_map(grid, maps, yslice=y_index)
plot_backward_map(grid, maps, yslice=y_index)
# poloidal_grid, y_val = grid.getPoloidalGrid(y_index)
# poloidal_grid.plot()


# can also plot poloidal grid objects directly -- use getPoloidalGrid(y_index) method
# of Grid object and then plot() method of PoloidalGrid object

vacuum_mgrid = "equilibria/mgrid_G1600-12-89_QA2e-1_Bxdl25-128x128x128_vecpot.nc"
plasma_mgrid = "equilibria/bmw_G1600-12-89-QA2e-1_Bxdl25_128x128x128.nc"

mgrid = MGRID(vacuum_mgrid, plasma_mgrid)



# and there's field line tracing

xpos = 8.3
zpos = 0.0
yperiod = 2*np.pi
plot_poincare(
    mgrid,
    xpos,
    zpos,
    yperiod,
    nplot=1,
    y_slices=None,
    revs=40,
    nover=20,
    interactive=True,
)






