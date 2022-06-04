import netCDF4 as nc 
import datetime
import numpy as np


this_year = 1980
file = 'C:/Users/zhiyuan/Code/nc/data/uwnd/uwnd.' + str(this_year) + '.nc'
dataset = nc.Dataset(file)
#dataset.variables.keys()

levels = dataset.variables['level']
lats = dataset.variables['lat']
lons = dataset.variables['lon']
timestamps = dataset.variables['time']
tmps = dataset.variables['air'][:].data # (time,level,lat,lon)
# b = datetime.datetime(1970, 1, 1, 0) - datetime.datetime(1800, 1, 1, 0) 
# dates = datetime.datetime.utcfromtimestamp((timestamps[0]-1490184)*3600)

steps = 10  # 时间步长
air_tendayly = []

for tendays in range(35):
    air_tendayly.append(sum(tmps[tendays*10:tendays*10+steps])/steps)

air_tendayly = np.array(air_tendayly)


steps = 365  # 时间步长
air_yearly = []
air_yearly.append(sum(tmps[0:365])/365)
air_yearly = np.array(air_yearly)

#旬
new_file = nc.Dataset('C:/Users/zhiyuan/Code/nc/output/tendayly/tmp/'+ str(this_year) +'_tmp_tendayly.nc', 'w', format='NETCDF4')
new_file.createDimension('level', len(levels))
new_file.createDimension('lat', len(lats))
new_file.createDimension('lon', len(lons))
new_file.createDimension('time', 35)
new_file.createDimension('tmp_tendayly', len(air_tendayly))

level = new_file.createVariable('level', np.float32, dimensions='level')
lat = new_file.createVariable('lat', np.float32, dimensions='lat')
lon = new_file.createVariable('lon', np.float32, dimensions='lon')
time = new_file.createVariable('time', np.float64, dimensions='time')
tmp_tendayly = new_file.createVariable('tmp_tendayly', np.float32, dimensions=('time', 'level', 'lat', 'lon'))

level.description = "高度层,hPa"
lat.description = "纬度，°"
lon.description = "经度，°"
time.long_name = "时间，旬的起始小时"
time.delta_t = "0000-00-10 00:00:00"
time.avg_period = "0001-00-00 00:00:00"
time.standard_name = "time"
time.axis = "T"
time.climatology = "climatology_bounds"
time.climo_period = "1980/01/01 - 1980/12/31"
time.units = 'hours since 1800-01-01 00:00:0.0'
time.start_time = "1800-01-01 00:00:00"
tmp_tendayly.description = "某旬某高度层某点平均气温，K"

new_file.variables['level'][:] = levels[:]
new_file.variables['lat'][:] = lats[:]
new_file.variables['lon'][:] = lons[:]
new_file.variables['time'][:] = timestamps[::10][0:35]
new_file.variables['tmp_tendayly'][:] = air_tendayly

new_file.close()




file = 'C:/Users/zhiyuan/Code/nc/output/tendayly/tmp/1980_tmp_tendayly.nc'
dataset = nc.Dataset(file)
dataset.variables['time']



import netCDF4 as nc 
file = 'C:/Users/zhiyuan/Code/nc/data/tmp/omega.sig995.mon.ltm.nc'
dataset = nc.Dataset(file)
dataset.variables.keys()
dataset.variables['time'].shape









#年
new_file = nc.Dataset('C:/Users/zhiyuan/Code/nc/output/yearly/tmp/'+ str(this_year) +'_tmp_tendayly.nc', 'w', format='NETCDF4')
new_file.createDimension('level', len(levels))
new_file.createDimension('lat', len(lats))
new_file.createDimension('lon', len(lons))
new_file.createDimension('tenday', 35)
new_file.createDimension('tmp_tendayly', len(air_tendayly))

level = new_file.createVariable('level', np.float32, dimensions='level')
lat = new_file.createVariable('lat', np.float32, dimensions='lat')
lon = new_file.createVariable('lon', np.float32, dimensions='lon')
tenday = new_file.createVariable('tenday', np.float32, dimensions='tenday')
tmp_tendayly = new_file.createVariable('tmp_tendayly', np.float32, dimensions=('tenday', 'level', 'lat', 'lon'))

level.description = "高度层,hPa"
lat.description = "纬度，°"
lon.description = "经度，°"
tenday.description = "旬数，当年第n旬"
tmp_tendayly.description = "某旬某高度层某点平均气温，K"

new_file.variables['level'][:] = levels[:]
new_file.variables['lat'][:] = lats[:]
new_file.variables['lon'][:] = lons[:]
new_file.variables['tenday'][:] = timestamps[::10][0:35]
new_file.variables['tmp_tendayly'][:] = air_tendayly

new_file.close()


