NETCDF_INCLUDE_DIR = /home/bpmelli/.linuxbrew/Cellar/netcdf/4.4.1.1_4/include/
NETCDF_LIB_DIR     = /home/bpmelli/.linuxbrew/Cellar/netcdf/4.4.1.1_4/lib/ \
                     -Wl,-rpath,/home/bpmelli/.linuxbrew/Cellar/netcdf/4.4.1.1_4/lib/
NETCDF_LIBS        = -lnetcdf-cxx4 -lnetcdf
CXX_FLAGS = -std=c++11 -gdwarf-4
F77 = f77
CXX = g++

all: w2nc rnc

fw: w2grads.f
	$(F77) -g -o fw w2grads.f

w2nc: w2nc.o
	$(CXX) $(CXX_FLAGS) -o w2nc w2nc.o -L $(NETCDF_LIB_DIR) $(NETCDF_LIBS)

w2nc.o: w2nc.cpp
	$(CXX) $(CXX_FLAGS) -I$(NETCDF_INCLUDE_DIR) -c w2nc.cpp

rnc: rnc.o
	$(CXX) $(CXX_FLAGS) -o rnc rnc.o -L $(NETCDF_LIB_DIR) $(NETCDF_LIBS)

rnc.o: rnc.cpp
	$(CXX) $(CXX_FLAGS) -I$(NETCDF_INCLUDE_DIR) -c rnc.cpp

clean:
	rm -f fw w2nc rnc a.out *.o *.ctl *.dat

