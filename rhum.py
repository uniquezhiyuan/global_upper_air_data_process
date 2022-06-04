import netCDF4 as nc 
import datetime
import numpy as np


def data_export(file_name):
    file = '/root/nc/data/rhum/' + file_name
    print(file + "  已读取。")
    dataset = nc.Dataset(file)
    # dataset.variables.keys()
    
    levels = dataset.variables['level']
    lats = dataset.variables['lat']
    lons = dataset.variables['lon']
    timestamps = dataset.variables['time']
    rhums = dataset.variables['rhum'][:].data # (time,level,lat,lon)
    # b = datetime.datetime(1970, 1, 1, 0) - datetime.datetime(1800, 1, 1, 0) 
    # dates = datetime.datetime.utcfromtimestamp((timestamps[0]-1490184)*3600)
    
    #旬
    steps = 10  # 时间步长
    rhum_tendayly = []
    
    for tendays in range(36):
        rhum_tendayly.append(sum(rhums[tendays*10:tendays*10+steps])/steps)
    
    rhum_tendayly = np.array(rhum_tendayly)
    print(file + "  旬数据已生成。")
    
    #月
    steps = 30  # 时间步长
    rhum_monthly = []
    
    for months in range(12):
        rhum_monthly.append(sum(rhums[months*30:months*30+steps])/steps)
    
    rhum_monthly = np.array(rhum_monthly)
    print(file + "  月数据已生成。")
    
    #年
    steps = 365  # 时间步长
    rhum_yearly = []
    rhum_yearly.append(sum(rhums[0:365])/365)
    rhum_yearly = np.array(rhum_yearly)
    print(file + "  年数据已生成。")
    
    #旬
    tendayly_file = nc.Dataset('/root/nc/output/tendayly/rhum/'+ file_name[5:9] +'_rhum_tendayly.nc', 'w', format='NETCDF4')
    print(file_name + "  旬数据创建完成，准备写入。")
    tendayly_file.createDimension('level', len(levels))
    tendayly_file.createDimension('lat', len(lats))
    tendayly_file.createDimension('lon', len(lons))
    tendayly_file.createDimension('time', 36)
    tendayly_file.createDimension('rhum_tendayly', len(rhum_tendayly))
    
    level = tendayly_file.createVariable('level', np.float32, dimensions='level')
    lat = tendayly_file.createVariable('lat', np.float32, dimensions='lat')
    lon = tendayly_file.createVariable('lon', np.float32, dimensions='lon')
    time = tendayly_file.createVariable('time', np.float64, dimensions='time')
    rhum_t = tendayly_file.createVariable('rhum_tendayly', np.float32, dimensions=('time', 'level', 'lat', 'lon'))
    
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
    rhum_t.description = "某旬某高度层某点位势高度，米"
    
    tendayly_file.variables['level'][:] = levels[:]
    tendayly_file.variables['lat'][:] = lats[:]
    tendayly_file.variables['lon'][:] = lons[:]
    tendayly_file.variables['time'][:] = timestamps[::10][0:36]
    tendayly_file.variables['rhum_tendayly'][:] = rhum_tendayly
    
    tendayly_file.close()
    print(file_name + " 旬数据已写入。")
    
    
    #月
    monthly_file = nc.Dataset('/root/nc/output/monthly/rhum/'+ file_name[5:9] +'_rhum_monthly.nc', 'w', format='NETCDF4')
    print(file_name + "  月数据创建完成，准备写入。")
    monthly_file.createDimension('level', len(levels))
    monthly_file.createDimension('lat', len(lats))
    monthly_file.createDimension('lon', len(lons))
    monthly_file.createDimension('time', 12)
    monthly_file.createDimension('rhum_monthly', len(rhum_monthly))
    
    level = monthly_file.createVariable('level', np.float32, dimensions='level')
    lat = monthly_file.createVariable('lat', np.float32, dimensions='lat')
    lon = monthly_file.createVariable('lon', np.float32, dimensions='lon')
    time = monthly_file.createVariable('time', np.float64, dimensions='time')
    rhum_m = monthly_file.createVariable('rhum_monthly', np.float32, dimensions=('time', 'level', 'lat', 'lon'))
    
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
    rhum_m.description = "某月某高度层某点位势高度，米"
    
    monthly_file.variables['level'][:] = levels[:]
    monthly_file.variables['lat'][:] = lats[:]
    monthly_file.variables['lon'][:] = lons[:]
    monthly_file.variables['time'][:] = timestamps[::30][0:12]
    monthly_file.variables['rhum_monthly'][:] = rhum_monthly
    
    monthly_file.close()
    print(file_name + " 月数据已写入。")
    
    #年
    yearly_file = nc.Dataset('/root/nc/output/yearly/rhum/'+ file_name[5:9] +'_rhum_yearly.nc', 'w', format='NETCDF4')
    yearly_file.createDimension('level', len(levels))
    yearly_file.createDimension('lat', len(lats))
    yearly_file.createDimension('lon', len(lons))
    yearly_file.createDimension('time', 1)
    yearly_file.createDimension('rhum_yearly', len(rhum_yearly))
    
    level = yearly_file.createVariable('level', np.float32, dimensions='level')
    lat = yearly_file.createVariable('lat', np.float32, dimensions='lat')
    lon = yearly_file.createVariable('lon', np.float32, dimensions='lon')
    time = yearly_file.createVariable('time', np.float64, dimensions='time')
    rhum_y = yearly_file.createVariable('rhum_yearly', np.float32, dimensions=('time', 'level', 'lat', 'lon'))
    
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
    rhum_y.description = "某年某高度层某点位势高度，米"
    
    yearly_file.variables['level'][:] = levels[:]
    yearly_file.variables['lat'][:] = lats[:]
    yearly_file.variables['lon'][:] = lons[:]
    yearly_file.variables['time'][:] = timestamps[::365][0]
    yearly_file.variables['rhum_yearly'][:] = rhum_yearly
    
    yearly_file.close()
    print(file_name + " 年数据已写入。")


file_name_list = ["rhum." + str(1980+i) + ".nc" for i in range(42)]
for i in file_name_list:
    data_export(i)

