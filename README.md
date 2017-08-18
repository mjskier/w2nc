# w2nc
Convert radar data from W format to NetCDF

## TODO

* Convert to cmake

## Compile

* Fix location of NetCDF in the Makefile
* make

## Files

* w2nc.cpp

   Take a W file, and generate a NetCDF file

* w2grads.f
   Original Fortran program. Used as reference, and to compare the reading format.

* rnc.cpp
   Program to test that the 3D array in w2nc.cpp works correctly



