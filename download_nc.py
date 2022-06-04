
#https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/air.1980.nc
#https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/hgt.1980.nc
#https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/uwnd.1980.nc
#https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/vwnd.1980.nc
#https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/omega.1980.nc
#https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/air.1980.nc
#https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/rhum.1980.nc

#hgt
import wget

file_name_list = ["hgt." + str(1980+i) + ".nc" for i in range(43)]
hgt_file_path = "/root/nc/data/hgt/"

for i in file_name_list:
    wget.download("https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/" + i, hgt_file_path + i) 



#uwnd
import wget

file_name_list = ["uwnd." + str(1980+i) + ".nc" for i in range(43)]
uwnd_file_path = "/root/nc/data/uwind/"

for i in file_name_list:
    wget.download("https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/" + i, uwnd_file_path + i) 



#vwnd
import wget

file_name_list = ["vwnd." + str(1980+i) + ".nc" for i in range(42)]
vwnd_file_path = "/root/nc/data/vwnd/"

for i in file_name_list:
    wget.download("https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/" + i, vwnd_file_path + i) 


#omega
import wget

file_name_list = ["omega." + str(1980+i) + ".nc" for i in range(42)]
omega_file_path = "/root/nc/data/omega/"

for i in file_name_list:
    wget.download("https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/" + i, omega_file_path + i) 



#air
import wget

file_name_list = ["air." + str(1980+i) + ".nc" for i in range(42)]
air_file_path = "/root/nc/data/air/"

for i in file_name_list:
    wget.download("https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/" + i, air_file_path + i) 



#rhum
import wget

file_name_list = ["rhum." + str(1980+i) + ".nc" for i in range(42)]
rhum_file_path = "/root/nc/data/rhum/"

for i in file_name_list:
    wget.download("https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/" + i, rhum_file_path + i) 


