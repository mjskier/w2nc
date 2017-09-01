# w2nc
Convert radar data from W format to NetCDF, using this [Fortran code as reference](https://github.com/mmbell/HRD_Doppler_synthesis/blob/master/src/io/w2grads.f)


## Compile

* ```NETCDF_LIB_PATH=-L..... NETCDF_INCLUDE_PATH=-I..... cmake .```
* make

## Files

* w2nc.cpp

   Take a W file, and generate a NetCDF file

* rnc.cpp

   Program to test that the 3D array in w2nc.cpp works correctly



