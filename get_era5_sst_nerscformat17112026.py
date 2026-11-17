#
# This script gets ERA5 SST boundary condition data
#
#
# conversion to grib based on Ed's file: sst_conv.py

# importing libraries
import os
import cdsapi
import subprocess
import xarray as xr
import numpy as np
import iris
import glob
import iris.coord_categorisation as icc
import os.path
import itertools
import pandas as pd
import cftime
import subprocess
from iris.util import unify_time_units # used to concatenate
from iris.util import equalise_attributes # used to concatenate
from scipy.interpolate import griddata # used for array interpolation

# setting directories
outdir = '/home/users/vikkitho/20cr/ERA5'
tmpdir = '/home/users/vikkitho/20cr/ERA5/tmp'

# currently calculates for a single year
YYYY = '2025'
print('YEAR = ',YYYY)


# Download ERA5 data
print('beginning download')
# Produce list of months 
def generate_yrmonths():
	years = range(int(YYYY),int(YYYY)+1) # single year
	months = ['01','02','03','04','05','06','07','08','09','10','11','12']
	yrmonths = [ int("%s%s" % (year,month)) for year in years for month in months]
	return yrmonths

yrmonths = generate_yrmonths()
#yrmonths = ['202203'] # use this line to select specific YYYYMM (can be list)
c = cdsapi.Client()


for yrmonth in yrmonths:
    year = str(yrmonth)[:4]
    mon  = str(yrmonth)[4:6]
    print('Processing:',year,mon)
    tmpfile = os.path.join(tmpdir,'ERA_sst_3hrly_'+year+mon+'.nc')
    request = {'product_type': 'ensemble_members',
               'format': 'netcdf',
               #"download_format": "unarchived",
               'variable': 'sea_surface_temperature',
               'year': year,
               'month': mon,
               'day': ['01', '02', '03',
                       '04', '05', '06',
                       '07', '08', '09',
                       '10', '11', '12',
                       '13', '14', '15',
                       '16', '17', '18',
                       '19', '20', '21',
                       '22', '23', '24',
                       '25', '26', '27',
                       '28', '29', '30',
                       '31',],
               'time': ["00:00", "03:00", "06:00",
                        "09:00", "12:00", "15:00",
                        "18:00", "21:00"],}
    print(request)
    c.retrieve('reanalysis-era5-single-levels',request,tmpfile)
    print('Downloaded',tmpfile)
    # Daymean - cannot use cdo as multiple ensembles
    #ftas = os.path.join(outdir,'sst','ERA5_sst_day_'+year+mon+'.nc')
    #cdo_cmd = ['cdo','-O','-b','F32','daymean',tmpfile,ftas]
    #print(' '.join(cdo_cmd))
    #ret = subprocess.call(cdo_cmd)
    #if not ret==0:
    #    raise Exception('Error with cdo command')
    #os.remove(tmpfile)
    # instead using xarray
    x = xr.open_dataset(tmpfile, engine="netcdf4")
    z = x['sst']
    z.attrs['standard_name'] = 'sea_surface_temperature'
    sst = z.to_iris()
    print('Loaded as xarray')
    for ens in np.arange(10):
        ens0 = sst[ens]
        print('Extracted single ensemble member', ens)
        # save single ensemble as netcdf
        iris.save(ens0, tmpdir+"/ERA_sst_3hrly_"+year+mon+"_ens"+str(ens)+".nc")
        # convert to daily mean with cdo
        ftas = os.path.join(outdir,'sst','ERA5_sst_day_'+year+mon+'_ens'+str(ens)+'.nc')
        fin = os.path.join(tmpdir,'ERA_sst_3hrly_'+year+mon+'_ens'+str(ens)+'.nc')
        cdo_cmd = ['cdo','-O','-b','F32','daymean',fin,ftas]
        print(' '.join(cdo_cmd))
        ret = subprocess.call(cdo_cmd)
        if not ret==0:
            raise Exception('Error with cdo command')
        os.remove(tmpdir+"/ERA_sst_3hrly_"+year+mon+"_ens"+str(ens)+".nc")    
    os.remove(tmpfile)
    print('Finished data download')



# Create single SST file per year per ensemble member (10 files)
print('Starting NetCDF & grib creation')
for ens in np.arange(10):
    print('ENS: ',ens)
    EN = "%02d" % (ens + 1,)
    
    month_list = iris.cube.CubeList([])
    for mon in ['01','02','03','04','05','06','07','08','09','10','11','12']:
        print('MONTH: ',mon)
        era=iris.load('/home/users/vikkitho/20cr/ERA5/sst/ERA5_sst_day_'+YYYY+mon+'_ens'+str(ens)+'.nc')[0] # outputs of get_era5.py
        day_list = iris.cube.CubeList([])
        for day in np.arange(np.shape(era)[0]):
            # Changing from mask to NaN over land/ice
            day1 = era[day,...]
            day1.data = day1.data.data
            day1.data
            for i in np.arange(361): # lat
                for j in np.arange(720): # lon
                    if day1.data[i,j]>400:
                        day1.data[i,j] = np.nan
            # Interpolation over land to infill
            a = day1.data
            x, y = np.indices(a.shape)
            interp = np.array(a)
            interp[np.isnan(interp)] = griddata(
                (x[~np.isnan(a)], y[~np.isnan(a)]), # points we know
                a[~np.isnan(a)],                    # values we know
                (x[np.isnan(a)], y[np.isnan(a)]))   # points to interpolate
            # Data back into cube format
            day1.data = interp
            # interpolate over Antarctica - taking most southerly value and repeating for all nan
            for i in np.arange(720): #loop over longitude
                ind = np.argwhere(np.isnan(day1.data[:,i]))[0]-1
                for each in np.argwhere(np.isnan(day1.data[:,i])):
                    day1.data[each,i] = day1.data[ind,i]
            # regrid era5 to nersc grid
            nersc_grid = iris.load_cube('/home/users/vikkitho/20cr/NERSC2009/sst_01.nc')[0,...] # generated using >cdo -f nc copy sst_01.grib sst_01.nc
            day1_grid = day1.regrid(nersc_grid, iris.analysis.Linear())
            day_list.append(day1_grid)
        month_list.append(day_list.merge_cube())
    removed_attributes = equalise_attributes(month_list)
    year_cube = month_list.concatenate_cube()
    # swap ERA5 data into a copied NERSC cube (to match up the variables)
    ner = iris.load_cube('/home/users/vikkitho/20cr/NERSC2009/sst_01.nc')
    e2 = ner.copy() # cube with NERSC variables/parameters
    e2.data = year_cube.data # move ERA5 data into NERSC formatting
    # save as netcdf
    #iris.save(e2, '/gws/nopw/j04/uk20cr/boundary_conditions_ERA/ERA'+str(YYYY)+'/sst_'+EN+'.nc')
    iris.save(e2, '/home/users/vikkitho/20cr/ERA5/ERA'+str(YYYY)+'/sst_'+str(EN)+'.nc')
    print(f"  Finished ensemble member {EN}, saved NETCDF: {'ERA'+str(YYYY)+'/sst_'+EN+'.nc'}")
    
    # save as grib (my method)
    #cdo_cmd = ['cdo','-b','16','-f','grb','setmissval,0','ERA5/ERA'+year+'/sst_'+str(EN)+'.nc','ERA5/ERA'+year+'/sst_'+str(EN)+'.grib']
    #ret = subprocess.call(cdo_cmd)
    #if not ret==0:
    #    raise Exception('Error with cdo command')
    
    # save as grib (Ed's method)
    # --- 0. Configuration ---
    #nc_file_path = '/gws/nopw/j04/uk20cr/boundary_conditions_ERA/ERA'+str(YYYY)+'/sst_'+EN+'.nc'
    nc_file_path = '/home/users/vikkitho/20cr/ERA5/ERA'+str(YYYY)+'/sst_'+str(EN)+'.nc'
    output_dir = 'grib_ensemble_outputs'
    os.makedirs(output_dir, exist_ok=True)
    # --- 1. Load NetCDF ---
    ds_final = xr.open_dataset(nc_file_path)
    print(f"\nProcessing to grib")
    ds_final["var11"].attrs.update({
        "shortName": "sst",
        "paramId": 34,
        "units": "K",
        "long_name": "Sea Surface Temperature",
        "standard_name": "sea_surface_temperature"
    })
    ds_final.lat.attrs['standard_name'] = 'latitude'
    ds_final.lon.attrs['standard_name'] = 'longitude'
    # Ensure latitude descending for GRIB1
    ds_final = ds_final.sortby('lat', ascending=False)
    output_file = '/gws/nopw/j04/uk20cr/boundary_conditions_ERA/ERA'+str(YYYY)+'/sst_' + EN + '.grib'
    output_file = '/home/users/vikkitho/20cr/ERA5/ERA'+str(YYYY)+'/sst_'+str(EN)+'.grib'
    # --- Prepare GRIB output ---
    grib_file_path = os.path.join(
        output_file
    )
    if os.path.exists(grib_file_path):
        os.remove(grib_file_path)
    tmp_files = []
    for t_idx, time_val in enumerate(ds_final['time'].values):
        print(f"    Writing timestep {t_idx+1}/{len(ds_final['time'])}")
        nc_slice = os.path.join(output_dir, f"tmp_{t_idx:03d}.nc")
        tmp_grb = os.path.join(output_dir, f"tmp_{t_idx:03d}.grb")
        # Save single timestep
        ds_final.isel(time=t_idx).to_netcdf(nc_slice)
        # Extract date for settaxis 
        dt = np.datetime64(time_val).astype("M8[s]").astype("O")
        dt2013 = dt.replace(year=int(YYYY)) # 2013 used as originally it was for that year, irrelevant now.
        date_str = dt2013.strftime("%Y-%m-%d")
        # Convert to GRIB1 with 12-bit packing and correct date
        subprocess.run(
            [
                "cdo", "-b", "12", "-f", "grb1",
                f"settaxis,{date_str},00:00:00,1day",
                "-copy", nc_slice, tmp_grb
            ],
            check=True
        )
        os.remove(nc_slice)
        # Force TMP metadata and analysis step
        tmp_fixed = tmp_grb + ".fixed"
        subprocess.run(
            [
                "grib_set",
                "-s", "edition=1,centre=7,subCentre=0,table2Version=2,"
                      "indicatorOfParameter=11,levelType=1,step=0,changeDecimalPrecision=2",
                tmp_grb, tmp_fixed
            ],
            check=True
        )
        os.remove(tmp_grb)
        tmp_files.append(tmp_fixed)
    # Concatenate all timesteps into final GRIB1 file
    with open(grib_file_path, "wb") as outfile:
        for f in tmp_files:
            with open(f, "rb") as infile:
                outfile.write(infile.read())
            os.remove(f)
    print(f"  Finished ensemble member {EN}, saved GRIB1: {grib_file_path}")
