import numpy as np
import h5py 

months          = '01','02','03','04','05','06','07','08','09','10','11','12'
sources        = 'SAVA','BORF','TEMF','DEFO','PEAT','AGRI'
#sources         = 'AGRI'
pollutant       = 'PM25'
outfilename     = 'AGRI_PM25_emissions.nc'

# Calculate air quality-relevant emissions from each of the 14 GFED basis
# regions for each year.

directory    = '..'


"""
Read in emission factors
"""
species = [] # names of the different gas and aerosol species
EFs     = np.zeros((41, 6)) # 41 species, 6 sources

k = 0
f = open(directory+'/GFED4_Emission_Factors.txt')
while 1:
    line = f.readline()
    if line == "":
        break
        
    if line[0] != '#':
        contents = line.split()
        species.append(contents[0])
        EFs[k,:] = contents[1:]
        k += 1
                
f.close()

# we are interested in CO for this example (4th row):
#EF_CO = EFs[3,]
# EF_NOx = EFs[7,:]
EF_PM25 = EFs[9,:]
# EF_SOx = EFs[14,:]
# EF_NH3 = EFs[33,:]
# (I also want VOCs)

start_year = 2015
end_year   = 2016



"""
make table with summed DM emissions for each region, year, and source
"""
CO_table = np.zeros((15, end_year - start_year + 1)) # region, year

for year in range(start_year, end_year+1):
    string = directory+'/GFED4.1s_'+str(year)+'.hdf5'
    f = h5py.File(string, 'r')
    
    
    if year == start_year: # these are time invariable    
        basis_regions = f['/ancill/basis_regions'][:]
        grid_area     = f['/ancill/grid_cell_area'][:]
    
    
    CO_emissions = np.zeros((720, 1440))
    for month in range(12):
        # read in DM emissions
        string = '/emissions/'+months[month]+'/DM'
        DM_emissions = f[string][:]
#        for source in range(len(sources)):
#            # read in the fractional contribution of each source
#            string = '/emissions/'+months[month]+'/partitioning/DM_'+sources[source]
#            contribution = f[string][:]
#            # calculate CO emissions as the product of DM emissions (kg DM per 
#            # m2 per month), the fraction the specific source contributes to 
#            # this (unitless), and the emission factor (g CO per kg DM burned)
#            CO_emissions += DM_emissions * contribution * EF_PM25[source]
        # if you only want to save out one source:
        # read in the fractional contribution of each source
        string = '/emissions/'+months[month]+'/partitioning/DM_'+'AGRI'
        contribution = f[string][:]
        # calculate CO emissions as the product of DM emissions (kg DM per 
        # m2 per month), the fraction the specific source contributes to 
        # this (unitless), and the emission factor (g CO per kg DM burned)
        CO_emissions += DM_emissions * contribution * EF_PM25[5]


print(CO_emissions)
CO_emissions_flip = np.flipud(CO_emissions)

import netCDF4
ny, nx = (720, 1440)
lon = np.linspace(-179.875,179.875,nx)
lat = np.linspace(89.875,-89.875,ny)

ncout = netCDF4.Dataset(outfilename,'w','NETCDF3') # using netCDF3 for output format 
ncout.createDimension('lon',nx)
ncout.createDimension('lat',ny)
lonvar = ncout.createVariable('lon','float32',('lon'))
lonvar[:] = lon
latvar = ncout.createVariable('lat','float32',('lat'))
latvar[:] = lat
myvar = ncout.createVariable('myvar','float32',('lat','lon'))
myvar.units = "g year-1"
myvar.long_name = "CO emissions from all sources"
myvar.setncattr('units','mm')
myvar[:] = CO_emissions_flip

ncout.close();

