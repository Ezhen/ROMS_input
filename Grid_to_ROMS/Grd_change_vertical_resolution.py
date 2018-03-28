import pyroms
d = pyroms.grid.get_ROMS_grid('Parent')

theta_b = 3
theta_s = 7
Tcline=10
N = 15
vgrd = pyroms.vgrid.s_coordinate_4(d.vgrid.h, theta_b, theta_s, Tcline, N, hraw = d.vgrid.hraw)

grd = pyroms.grid.ROMS_Grid(d.name, d.hgrid, vgrd)

print 'write grid to netcdf file'
pyroms.grid.write_ROMS_grid(grd, filename='Parent15.nc')
