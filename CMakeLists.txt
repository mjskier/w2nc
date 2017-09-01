cmake_minimum_required (VERSION 2.8.9)
project (w2nc)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set (NETCDF_CXX 1)

# The NETCDF variables set by the NetCDF package don't work in my environment.
# (MacOS, netcdf installed via homebrew)
# The code uses the new C++ interface to netcdf4, so it needs libnetcdf-cxx4

# I set these on the command line when invoking "cmake ."
# NETCDF_LIB_PATH=.... NETCDF_INCLUDE_PATH=... cmake .

# You could hard code them here instead...
# set (NETCDF_LIB_PATH "-L/usr/local/Cellar/netcdf/4.4.1.1_5/lib")
# set (NETCDF_INCLUDE_PATH "-I/usr/local/Cellar/netcdf/4.4.1.1_5/include")

find_package (NetCDF REQUIRED)

include_directories(${NETCDF_INCLUDE_PATH})

add_executable (w2nc w2nc.cpp)
target_link_libraries (w2nc ${NETCDF_LIB_PATH} -lnetcdf-cxx4 -lnetcdf )

add_executable (rnc rnc.cpp)
target_link_libraries (rnc ${NETCDF_LIB_PATH} -lnetcdf-cxx4 -lnetcdf )

# To see what variables were set by the various modules, uncomment the following
# get_cmake_property(_variableNames VARIABLES)
# foreach (_variableName ${_variableNames})
#     message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()

# To print environment variables use cmake's command mode:
# execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "environment")
