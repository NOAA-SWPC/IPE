import numpy as np
from netCDF4 import Dataset
from multiprocessing import Pool
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import glob
from os import path

class Grid:
  def __init__(self, path):
    # nlp, nmp, nfluxtube, nlon_geo, nlat_geo, nheights_geo
    # flux_tube_max(lp), facfac_interface, ii[1-4]_interface, dd_interface
    g = Dataset(path).groups['apex_grid']

    self.nlp          = len(g.dimensions['phony_dim_0'])
    self.nfluxtube    = len(g.dimensions['phony_dim_1'])
    self.nmp          = len(g.dimensions['phony_dim_2'])
    self.nlon_geo     = len(g.dimensions['phony_dim_3'])
    self.nlat_geo     = len(g.dimensions['phony_dim_4'])
    self.nheights_geo = len(g.dimensions['phony_dim_5'])

    self.flux_tube_max = g.variables['tube_max'][:]

    self.longitude_geo = np.zeros(self.nlon_geo)
    self.latitude_geo  = np.zeros(self.nlat_geo)
    self.altitude_geo  = np.zeros(self.nheights_geo)

    for i in range(self.nlon_geo):
      self.longitude_geo[i] = i * 360. / self.nlon_geo

    for i in range(self.nlat_geo):
      self.latitude_geo[i]  = -90. + i * 180. / (self.nlat_geo-1)

    for i in range(self.nheights_geo):
      self.altitude_geo[i]  = i*5. + 90.

    self.facfac_interface = np.zeros( (3, self.nlon_geo, self.nlat_geo, self.nheights_geo) )
    self.dd_interface     = np.zeros( (3, self.nlon_geo, self.nlat_geo, self.nheights_geo) )
    self.ii1_interface    = np.zeros( (3, self.nlon_geo, self.nlat_geo, self.nheights_geo) )
    self.ii3_interface    = np.zeros( (3, self.nlon_geo, self.nlat_geo, self.nheights_geo) )
    self.ii4_interface    = np.zeros( (3, self.nlon_geo, self.nlat_geo, self.nheights_geo) )

    self.facfac_interface[0] = g.variables['fac_1'][:]
    self.facfac_interface[1] = g.variables['fac_2'][:]
    self.facfac_interface[2] = g.variables['fac_3'][:]

    self.dd_interface[0] = 1. / g.variables['dd_1'][:]
    self.dd_interface[1] = 1. / g.variables['dd_2'][:]
    self.dd_interface[2] = 1. / g.variables['dd_3'][:]
    self.dtot_inv = np.sum(self.dd_interface,axis=0)

    self.ii1_interface[0] = g.variables['ii1_1'][:]
    self.ii1_interface[1] = g.variables['ii1_2'][:]
    self.ii1_interface[2] = g.variables['ii1_3'][:]

    self.ii3_interface[0] = g.variables['ii3_1'][:]
    self.ii3_interface[1] = g.variables['ii3_2'][:]
    self.ii3_interface[2] = g.variables['ii3_3'][:]

    self.ii4_interface[0] = g.variables['ii4_1'][:]
    self.ii4_interface[1] = g.variables['ii4_2'][:]
    self.ii4_interface[2] = g.variables['ii4_3'][:]

    self.npts = np.sum(self.flux_tube_max*self.nlp)

    npts = 0
    self.ilp_map = np.zeros( (2, self.npts) )
    for lp in range(self.nlp):
      for i in range(self.flux_tube_max[lp]):
        self.ilp_map[0,npts] = i
        self.ilp_map[1,npts] = lp
        npts += 1

  def interpolate_to_geogrid(self, apex_data):
    geo_data = np.zeros( (self.nlon_geo, self.nlat_geo, self.nheights_geo) )

    factor = self.facfac_interface


    for i in range(3):
      mp  = (self.ii1_interface[i]-1).astype(np.int32) # nlon, nlat, nheights
      in2 = (self.ii3_interface[i]-1).astype(np.int32) # nlon, nlat, nheights
      in1 = (self.ii4_interface[i]-1).astype(np.int32) # nlon, nlat, nheights

      iFlux2 = self.ilp_map[0,in2].astype(np.int32)    # nlon, nlat, nheights
      lp2    = self.ilp_map[1,in2].astype(np.int32)    # nlon, nlat, nheights

      iFlux1 = self.ilp_map[0,in1].astype(np.int32)    # nlon, nlat, nheights
      lp1    = self.ilp_map[1,in1].astype(np.int32)    # nlon, nlat, nheights

      geo_data += ((apex_data[mp,lp2,iFlux2] - apex_data[mp,lp1,iFlux1])*factor[i] + apex_data[mp,lp1,iFlux1]) * self.dd_interface[i]

    return np.swapaxes(geo_data / self.dtot_inv,0,2)

class Plasma:
  def __init__(self, filename, nlp, nmp, nfluxtube):
    self.ion_densities = np.zeros( (9, nmp, nlp, nfluxtube) )

    f = Dataset(filename).groups['apex']

    self.ion_densities[0] = f.variables['o_plus_density'][:]
    self.ion_densities[1] = f.variables['h_plus_density'][:]
    self.ion_densities[2] = f.variables['he_plus_density'][:]
    self.ion_densities[3] = f.variables['n_plus_density'][:]
    self.ion_densities[4] = f.variables['no_plus_density'][:]
    self.ion_densities[5] = f.variables['o2_plus_density'][:]
    self.ion_densities[6] = f.variables['n2_plus_density'][:]
    self.ion_densities[7] = f.variables['o_plus_2D_density'][:]
    self.ion_densities[8] = f.variables['o_plus_2P_density'][:]
#   self.tec              = np.sum(self.ion_densities,axis=0)
    self.electron_density = (np.sum(self.ion_densities,axis=0)) /1.e12

class IPE:
  def __init__(self, grid_filename):
    self.grid = Grid(grid_filename)

  def read_h5(self, h5_filename):
    self.plasma  = Plasma(h5_filename,  self.grid.nlp, self.grid.nmp, self.grid.nfluxtube)

  def write_netcdf(self, netcdf_filename):
    # Open
    o = Dataset(netcdf_filename, 'w', format='NETCDF4_CLASSIC')
    # Dimensions
    z_dim = o.createDimension('altitude',  self.grid.nheights_geo)

    y_dim = o.createDimension('latitude',  self.grid.nlat_geo)

    x_dim = o.createDimension('longitude', self.grid.nlon_geo)

    z_var = o.createVariable('altitude',  'f4', 'altitude',zlib=True,least_significant_digit=3)
    z_var.long_name = 'Altitude'
    z_var.units     = 'km'

    y_var = o.createVariable('latitude',  'f4', 'latitude',zlib=True,least_significant_digit=3)
    y_var.long_name = 'Latitude'
    y_var.units     = 'degrees_north'

    x_var = o.createVariable('longitude', 'f4', 'longitude',zlib=True,least_significant_digit=3)
    x_var.long_name = 'Longitude'
    x_var.units     = 'degrees_east'

    ne_var = o.createVariable('electron_density',       'f4', ('altitude','latitude','longitude',),zlib=True ,least_significant_digit=4) #,complevel=9)
    ne_var.long_name = "Electron Density"
    ne_var.units     = "m^{-3}"

    z_var[:] = self.grid.altitude_geo
    y_var[:] = self.grid.latitude_geo
    x_var[:] = self.grid.longitude_geo

    ne_var[:]  = self.grid.interpolate_to_geogrid(self.plasma.electron_density)

    o.close()

def load_and_write(i):
  print files[i]
  timestamp = files[i][-15:-3]
  print timestamp
  ipe.read_h5(files[i])
  ipe.write_netcdf(path.join(args.outdir,"IPE_Ne.geo.{}.nc4".format(timestamp)))

## input parsing options
parser = ArgumentParser(description='Diff two NetCDF files as defined in this script', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-g', '--gridfile', help='path to IPE_Grid.h5',      type=str, required=True)
parser.add_argument('-i', '--indir',    help='path to input directory',  type=str, default=".")
parser.add_argument('-o', '--outdir',   help='path to output directory', type=str, default="output")
args = parser.parse_args()

MAX_PROCS = 8

ipe = IPE(args.gridfile)
files = glob.glob(path.join(args.indir,"IPE_State.apex.*.h5"))
p = Pool(min([len(files),MAX_PROCS]))
p.map(load_and_write,range(len(files)))
#load_and_write(0)

