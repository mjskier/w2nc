// Test the 3D array mechanism in w2nc.cpp
// Read and dump the dbz table.
//
// If the commented code at the end of Data::write method in w2nc.cpp is
//   uncommented, there will be a known pattern in the dbz table.
//
// Bruno Melli 8/15/17

#include <unistd.h>
#include <stdlib.h>

#include <ncFile.h>
#include <ncDim.h>
#include <ncVar.h>

typedef unsigned long u_long;

void usage(const char *s) {
  std::cout << "Usage: " << s << "-i <input W file>"
	    << std::endl;
}

int main(int argc, char *argv[])
{
  int opt;
  char *in = NULL, *out = NULL;
  
  while( (opt = getopt(argc, argv, "i:o:h")) != -1)
    switch(opt){
    case 'i':
      in = strdup(optarg);
      break;
    case 'h':
    case '?':
      usage(argv[0]);
      exit(0);
    }

  // Make sure we have the right options

  if ( in == NULL ) {
    usage(argv[0]);
    exit(1);
  }

  netCDF::NcFile file(in, netCDF::NcFile::read);
  // How can you tell if the file didn't exist???
  
  netCDF::NcDim xDim = file.getDim("xDim");
  netCDF::NcDim yDim = file.getDim("yDim");
  netCDF::NcDim zDim = file.getDim("zDim");
  
  int iDim = xDim.getSize();
  int jDim = yDim.getSize();
  int kDim = zDim.getSize();

  std::cout << "iDim: " << iDim << ", jDim: " << jDim << ", kDim: " << kDim << std::endl;

  netCDF::NcVar dbz = file.getVar("dbz");
  float val;
  for(u_long i = 0; i < iDim; i++)
    for(u_long j = 0; j < jDim; j++)
      for(u_long k = 0; k < kDim; k++) {
	dbz.getVar({k, j, i}, &val);
	std::cout << i << ", " << j << ", " << k << ": " << val << std::endl;
      }
}

