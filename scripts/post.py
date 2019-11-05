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
    self.electron_temperature = f.variables['electron_temperature'][:]
    self.ion_temperature      = f.variables['ion_temperature'][:]

class Neutral:
  def __init__(self, filename, nlp, nmp, nfluxtube):
    self.velocity_apex = np.zeros( (3, nmp, nlp, nfluxtube) )

    f = Dataset(filename).groups['apex']

    self.oxygen              = f.variables['o_density'][:]
    self.hydrogen            = f.variables['h_density'][:]
    self.helium              = f.variables['he_density'][:]
    self.nitrogen            = f.variables['n_density'][:]
    self.molecular_oxygen    = f.variables['o2_density'][:]
    self.molecular_nitrogen  = f.variables['n2_density'][:]
    self.neutral_temperature = f.variables['neutral_temperature'][:]
    self.velocity_apex[0]    = f.variables['neutral_apex1_velocity'][:]
    self.velocity_apex[1]    = f.variables['neutral_apex2_velocity'][:]
    self.velocity_apex[2]    = f.variables['neutral_apex3_velocity'][:]

class IPE:
  def __init__(self, grid_filename):
    self.grid = Grid(grid_filename)

  def read_h5(self, h5_filename):
    self.plasma  = Plasma(h5_filename,  self.grid.nlp, self.grid.nmp, self.grid.nfluxtube)
    self.neutral = Neutral(h5_filename, self.grid.nlp, self.grid.nmp, self.grid.nfluxtube)

  def write_netcdf(self, netcdf_filename):
    # Open
    o = Dataset(netcdf_filename, 'w', format='NETCDF4_CLASSIC')
    # Dimensions
    z_dim = o.createDimension('altitude',  self.grid.nheights_geo)

    y_dim = o.createDimension('latitude',  self.grid.nlat_geo)

    x_dim = o.createDimension('longitude', self.grid.nlon_geo)

    z_var = o.createVariable('altitude',  'f4', 'altitude')
    z_var.long_name = 'Altitude'
    z_var.units     = 'km'

    y_var = o.createVariable('latitude',  'f4', 'latitude')
    y_var.long_name = 'Latitude'
    y_var.units     = 'degrees_north'

    x_var = o.createVariable('longitude', 'f4', 'longitude')
    x_var.long_name = 'Longitude'
    x_var.units     = 'degrees_east'

    # Variables
    op_var   = o.createVariable('o_plus',               'f4', ('altitude','latitude','longitude',))
    op_var.long_name = "O+ number density"
    op_var.units     = "m^{-3}"

    hp_var   = o.createVariable('h_plus',               'f4', ('altitude','latitude','longitude',))
    hp_var.long_name = "H+ number density"
    hp_var.units     = "m^{-3}"

    hep_var  = o.createVariable('he_plus',              'f4', ('altitude','latitude','longitude',))
    hep_var.long_name = "He+ number density"
    hep_var.units     = "m^{-3}"

    np_var   = o.createVariable('n_plus',               'f4', ('altitude','latitude','longitude',))
    np_var.long_name = "N+ number density"
    np_var.units     = "m^{-3}"

    nop_var  = o.createVariable('no_plus',              'f4', ('altitude','latitude','longitude',))
    nop_var.long_name = "NO+ number density"
    nop_var.units     = "m^{-3}"

    o2p_var  = o.createVariable('o2_plus',              'f4', ('altitude','latitude','longitude',))
    o2p_var.long_name = "O2+ number density"
    o2p_var.units     = "m^{-3}"

    n2p_var  = o.createVariable('n2_plus',              'f4', ('altitude','latitude','longitude',))
    n2p_var.long_name = "N2+ number density"
    n2p_var.units     = "m^{-3}"

    op2d_var = o.createVariable('o_plus_2d',            'f4', ('altitude','latitude','longitude',))
    op2d_var.long_name = "O+(2D) number density"
    op2d_var.units     = "m^{-3}"

    op2p_var = o.createVariable('o_plus_2p',            'f4', ('altitude','latitude','longitude',))
    op2p_var.long_name = "O+(2P) number density"
    op2p_var.units     = "m^{-3}"

    he_var   = o.createVariable('helium',               'f4', ('altitude','latitude','longitude',))
    he_var.long_name = "Neutral He density"
    he_var.units     = "kg m^{-3}"

    o_var    = o.createVariable('oxygen',               'f4', ('altitude','latitude','longitude',))
    o_var.long_name = "Neutral O density"
    o_var.units     = "kg m^{-3}"

    o2_var   = o.createVariable('molecular_oxygen',     'f4', ('altitude','latitude','longitude',))
    o2_var.long_name = "Neutral O2 density"
    o2_var.units     = "kg m^{-3}"

    n2_var   = o.createVariable('molecular_nitrogen',   'f4', ('altitude','latitude','longitude',))
    n2_var.long_name = "Neutral N2 density"
    n2_var.units     = "kg m^{-3}"

    n_var    = o.createVariable('nitrogen',             'f4', ('altitude','latitude','longitude',))
    n_var.long_name = "Neutral N density"
    n_var.units     = "kg m^{-3}"

    h_var    = o.createVariable('hydrogen',             'f4', ('altitude','latitude','longitude',))
    h_var.long_name = "Neutral H density"
    h_var.units     = "kg m^{-3}"

    t_var    = o.createVariable('temperature',          'f4', ('altitude','latitude','longitude',))
    t_var.long_name = "Neutral temperature"
    t_var.units     = "K"

    u_var    = o.createVariable('u',                    'f4', ('altitude','latitude','longitude',))
    u_var.long_name = "Apex1 Velocity"
    u_var.units     = "m s^{-1}"

    v_var    = o.createVariable('v',                    'f4', ('altitude','latitude','longitude',))
    v_var.long_name = "Apex2 Velocity"
    v_var.units     = "m s^{-1}"

    w_var    = o.createVariable('w',                    'f4', ('altitude','latitude','longitude',))
    w_var.long_name = "Apex3 Velocity"
    w_var.units     = "m s^{-1}"

    it_var   = o.createVariable('ion_temperature',      'f4', ('altitude','latitude','longitude',))
    it_var.long_name = "Ion temperature"
    it_var.units     = "K"

    et_var   = o.createVariable('electron_temperature', 'f4', ('altitude','latitude','longitude',))
    et_var.long_name = "Electron temperature"
    et_var.units     = "K"

    z_var[:] = self.grid.altitude_geo
    y_var[:] = self.grid.latitude_geo
    x_var[:] = self.grid.longitude_geo

    op_var[:]   = self.grid.interpolate_to_geogrid(self.plasma.ion_densities[0])
    hp_var[:]   = self.grid.interpolate_to_geogrid(self.plasma.ion_densities[1])
    hep_var[:]  = self.grid.interpolate_to_geogrid(self.plasma.ion_densities[2])
    np_var[:]   = self.grid.interpolate_to_geogrid(self.plasma.ion_densities[3])
    nop_var[:]  = self.grid.interpolate_to_geogrid(self.plasma.ion_densities[4])
    o2p_var[:]  = self.grid.interpolate_to_geogrid(self.plasma.ion_densities[5])
    n2p_var[:]  = self.grid.interpolate_to_geogrid(self.plasma.ion_densities[6])
    op2d_var[:] = self.grid.interpolate_to_geogrid(self.plasma.ion_densities[7])
    op2p_var[:] = self.grid.interpolate_to_geogrid(self.plasma.ion_densities[8])
    he_var[:]   = self.grid.interpolate_to_geogrid(self.neutral.helium)
    o_var[:]    = self.grid.interpolate_to_geogrid(self.neutral.oxygen)
    o2_var[:]   = self.grid.interpolate_to_geogrid(self.neutral.molecular_oxygen)
    n2_var[:]   = self.grid.interpolate_to_geogrid(self.neutral.molecular_nitrogen)
    n_var[:]    = self.grid.interpolate_to_geogrid(self.neutral.nitrogen)
    h_var[:]    = self.grid.interpolate_to_geogrid(self.neutral.hydrogen)
    t_var[:]    = self.grid.interpolate_to_geogrid(self.neutral.neutral_temperature)
    u_var[:]    = self.grid.interpolate_to_geogrid(self.neutral.velocity_apex[0])
    v_var[:]    = self.grid.interpolate_to_geogrid(self.neutral.velocity_apex[1])
    w_var[:]    = self.grid.interpolate_to_geogrid(self.neutral.velocity_apex[2])
    it_var[:]   = self.grid.interpolate_to_geogrid(self.plasma.ion_temperature)
    et_var[:]   = self.grid.interpolate_to_geogrid(self.plasma.electron_temperature)

    o.close()

def load_and_write(i):
  print files[i]
  timestamp = files[i][-15:-3]
  print timestamp
  ipe.read_h5(files[i])
  ipe.write_netcdf(path.join(args.outdir,"IPE_Params.geo.{}.nc4".format(timestamp)))

## input parsing options
parser = ArgumentParser(description='Diff two NetCDF files as defined in this script', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-g', '--gridfile', help='path to IPE_Grid.h5',      type=str, required=True)
parser.add_argument('-i', '--indir',    help='path to input directory',  type=str, default=".")
parser.add_argument('-o', '--outdir',   help='path to output directory', type=str, default="output")
args = parser.parse_args()

MAX_PROCS = 4

ipe = IPE(args.gridfile)
files = glob.glob(path.join(args.indir,"IPE_State.apex.*.h5"))
p = Pool(min([len(files),MAX_PROCS]))
p.map(load_and_write,range(len(files)))
