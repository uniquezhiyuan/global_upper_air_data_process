import netCDF4 as nc 
import datetime
import numpy as np


def data_export(file_name):
    file = '/root/nc/data/omega/' + file_name
    print(file + "  已读取。")
    dataset = nc.Dataset(file)
    # dataset.variables.keys()
    
    levels = dataset.variables['level']
    lats = dataset.variables['lat']
    lons = dataset.variables['lon']
    timestamps = dataset.variables['time']
    omegas = dataset.variables['omega'][:].data # (time,level,lat,lon)
    # b = datetime.datetime(1970, 1, 1, 0) - datetime.datetime(1800, 1, 1, 0) 
    # dates = datetime.datetime.utcfromtimestamp((timestamps[0]-1490184)*3600)
    
    #旬
    steps = 10  # 时间步长
    omega_tendayly = []
    
    for tendays in range(36):
        omega_tendayly.append(sum(omegas[tendays*10:tendays*10+steps])/steps)
    
    omega_tendayly = np.array(omega_tendayly)
    print(file + "  旬数据已生成。")
    
    #月
    steps = 30  # 时间步长
    omega_monthly = []
    
    for months in range(12):
        omega_monthly.append(sum(omegas[months*30:months*30+steps])/steps)
    
    omega_monthly = np.array(omega_monthly)
    print(file + "  月数据已生成。")
    
    #年
    steps = 365  # 时间步长
    omega_yearly = []
    omega_yearly.append(sum(omegas[0:365])/365)
    omega_yearly = np.array(omega_yearly)
    print(file + "  年数据已生成。")
    
    #旬
    tendayly_file = nc.Dataset('/root/nc/output/tendayly/omega/'+ file_name[6:10] +'_omega_tendayly.nc', 'w', format='NETCDF4')
    print(file_name + "  旬数据创建完成，准备写入。")
    tendayly_file.createDimension('level', len(levels))
    tendayly_file.createDimension('lat', len(lats))
    tendayly_file.createDimension('lon', len(lons))
    tendayly_file.createDimension('time', 36)
    tendayly_file.createDimension('omega_tendayly', len(omega_tendayly))
    
    level = tendayly_file.createVariable('level', np.float32, dimensions='level')
    lat = tendayly_file.createVariable('lat', np.float32, dimensions='lat')
    lon = tendayly_file.createVariable('lon', np.float32, dimensions='lon')
    time = tendayly_file.createVariable('time', np.float64, dimensions='time')
    omega_t = tendayly_file.createVariable('omega_tendayly', np.float32, dimensions=('time', 'level', 'lat', 'lon'))
    
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
    omega_t.description = "某旬某高度层某点位势高度，米"
    
    tendayly_file.variables['level'][:] = levels[:]
    tendayly_file.variables['lat'][:] = lats[:]
    tendayly_file.variables['lon'][:] = lons[:]
    tendayly_file.variables['time'][:] = timestamps[::10][0:36]
    tendayly_file.variables['omega_tendayly'][:] = omega_tendayly
    
    tendayly_file.close()
    print(file_name + " 旬数据已写入。")
    
    
    #月
    monthly_file = nc.Dataset('/root/nc/output/monthly/omega/'+ file_name[6:10] +'_omega_monthly.nc', 'w', format='NETCDF4')
    print(file_name + "  月数据创建完成，准备写入。")
    monthly_file.createDimension('level', len(levels))
    monthly_file.createDimension('lat', len(lats))
    monthly_file.createDimension('lon', len(lons))
    monthly_file.createDimension('time', 12)
    monthly_file.createDimension('omega_monthly', len(omega_monthly))
    
    level = monthly_file.createVariable('level', np.float32, dimensions='level')
    lat = monthly_file.createVariable('lat', np.float32, dimensions='lat')
    lon = monthly_file.createVariable('lon', np.float32, dimensions='lon')
    time = monthly_file.createVariable('time', np.float64, dimensions='time')
    omega_m = monthly_file.createVariable('omega_monthly', np.float32, dimensions=('time', 'level', 'lat', 'lon'))
    
    level.description = "高度层,hPa"
    lat.description = "纬度，°"
    lon.description = "经度，°"
    time.long_name = "时间，月的起始小时"
    time.delta_t = "0000-00-10 00:00:00"
    time.avg_period = "0001-00-00 00:00:00"
    time.standard_name = "time"
    time.axis = "T"
    time.climatology = "climatology_bounds"
    time.climo_period = "1980/01/01 - 1980/12/31"
    time.units = 'hours since 1800-01-01 00:00:0.0'
    time.start_time = "1800-01-01 00:00:00"
    omega_m.description = "某月某高度层某点位势高度，米"
    
    monthly_file.variables['level'][:] = levels[:]
    monthly_file.variables['lat'][:] = lats[:]
    monthly_file.variables['lon'][:] = lons[:]
    monthly_file.variables['time'][:] = timestamps[::30][0:12]
    monthly_file.variables['omega_monthly'][:] = omega_monthly
    
    monthly_file.close()
    print(file_name + " 月数据已写入。")
    
    #年
    yearly_file = nc.Dataset('/root/nc/output/yearly/omega/'+ file_name[6:10] +'_omega_yearly.nc', 'w', format='NETCDF4')
    yearly_file.createDimension('level', len(levels))
    yearly_file.createDimension('lat', len(lats))
    yearly_file.createDimension('lon', len(lons))
    yearly_file.createDimension('time', 1)
    yearly_file.createDimension('omega_yearly', len(omega_yearly))
    
    level = yearly_file.createVariable('level', np.float32, dimensions='level')
    lat = yearly_file.createVariable('lat', np.float32, dimensions='lat')
    lon = yearly_file.createVariable('lon', np.float32, dimensions='lon')
    time = yearly_file.createVariable('time', np.float64, dimensions='time')
    omega_y = yearly_file.createVariable('omega_yearly', np.float32, dimensions=('time', 'level', 'lat', 'lon'))
    
    level.description = "高度层,hPa"
    lat.description = "纬度，°"
    lon.description = "经度，°"
    time.long_name = "时间，年的起始小时"
    time.delta_t = "0000-00-10 00:00:00"
    time.avg_period = "0001-00-00 00:00:00"
    time.standard_name = "time"
    time.axis = "T"
    time.climatology = "climatology_bounds"
    time.climo_period = "1980/01/01 - 1980/12/31"
    time.units = 'hours since 1800-01-01 00:00:0.0'
    time.start_time = "1800-01-01 00:00:00"
    omega_y.description = "某年某高度层某点位势高度，米"
    
    yearly_file.variables['level'][:] = levels[:]
    yearly_file.variables['lat'][:] = lats[:]
    yearly_file.variables['lon'][:] = lons[:]
    yearly_file.variables['time'][:] = timestamps[::365][0]
    yearly_file.variables['omega_yearly'][:] = omega_yearly
    
    yearly_file.close()
    print(file_name + " 年数据已写入。")


file_name_list = ["omega." + str(1980+i) + ".nc" for i in range(42)]
for i in file_name_list:
    data_export(i)

