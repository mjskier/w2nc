// Convert a W format file to a netcdf
// Bruno Melli 8/15/17

#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

#include <ncFile.h>
#include <ncDim.h>
#include <ncVar.h>

#include <iostream>
#include <fstream>
#include <ios>
#include <vector>

class Array3D {
  
  
public:
  
  Array3D(int xDim, int yDim, int zDim);
  
  float& operator() (int i, int j, int k)             { return *(table + i + dimX * (j + dimY * k)); }
  const float& operator() (int i, int j, int k) const { return *(table + i + dimX * (j + dimY * k)); }

  float *data() { return table; }
  
private:

  int dimX;
  int dimY;
  int dimZ;

  float *table;
};

Array3D::Array3D(int x, int y, int z)
{
  dimX = x;
  dimY = y;
  dimZ = z;
  
  table = new float[dimX * dimY * dimZ];
}

class Data {

public:
  
  static constexpr int NoData = 32767;
  static constexpr float DivScaleFactor = 100000.0;
  static constexpr float fillValue = -9999.0;
  static constexpr float Pi = 3.14159;
  
  Data(char *fname, bool debug = false);
  
  bool valid() { return valid_flag; }
  
  bool ldlnl();
  bool write(const char *fname);
  
private:

  bool fill_char(char *buf , size_t size);
  bool fill_int(int &val);
  bool fill_short(short &val);
  bool fill_float(float &val);
  bool allocate_tables();

  char *input_path;
  std::ifstream file_handle;

  bool valid_flag;
  bool fDebug;

  bool iuv; // true if keyword below is RUV or ruv;
  
  char keyword[5];
  char fltname[9];
  char stmname[13];
  char radar[5];
  char experiment[33];
  char creattime[33];
  char extra1[29];

  int imax, jmax, kmax, kount, nmosm, iunfld, iatten, flag, extra2,
    extra3, lfn;

  float stime, etime, olat, olon, sx, sy, sz, xz, yz, zz, rot, ra, co1,
    co2, azmcor, elcor, thresh, powert, biel, azbiel, etime1, stime2,
    extra6, extra7;
    
  Array3D *u_table, *v_table, *w_table, *dbz_table, *div_table;

  float *rwd_table, *rws_table, *rww_table, *rdv_table;
  short *wd_table, *ws_table, *db_table, *ww_table, *dv_table;
};

bool Data::fill_int(int &val)
{
  int ival;
  file_handle.read(reinterpret_cast<char *>(&ival), sizeof(int));
  val = ival;
  return true;
}

bool Data::fill_short(short &val)
{
  short ival;
  file_handle.read(reinterpret_cast<char *>(&ival), sizeof(short));
  val = ival;
  return true;
}

bool Data::fill_float(float &val)
{
  float fval;
  file_handle.read(reinterpret_cast<char *>(&fval), sizeof(float));
  val = fval;
  return true;
}

bool Data::fill_char(char *buf, size_t size)
{
  file_handle.read(buf, size);
  if (file_handle.gcount() != size) {
    std::cerr << "Byte reqested: " << size << ", read: " << file_handle.gcount() << std::endl;
    return false;
  }
  // trim trailing spaces
  for(char *s = buf + size; s != buf; s--)
    if(*s == ' ')
      *s = '\0';
  buf[size] = '\0';
  return true;
}

bool Data::allocate_tables()
{
  int size = imax * jmax *kmax;
  
  u_table   = new Array3D(imax, jmax, kmax);
  v_table   = new Array3D(imax, jmax, kmax);
  w_table   = new Array3D(imax, jmax, kmax);
  dbz_table = new Array3D(imax, jmax, kmax);
  div_table = new Array3D(imax, jmax, kmax);

  rwd_table = new float[imax];
  rws_table = new float[imax];
  rww_table = new float[imax];
  rdv_table = new float[imax];

  wd_table = new short[imax];
  ws_table = new short[imax];
  db_table = new short[imax];
  ww_table = new short[imax];
  dv_table = new short[imax];

  return true;
}

bool Data::ldlnl() 
{
  allocate_tables();

  short dummy;

  for(int k = 0; k < kmax; k++) {
    for(int j = 0; j < jmax; j++) {
      
      // Get past the Fortran padding
      for(int d = 0; d < 4; d++) 
	fill_short(dummy);

      for(int i = 0; i < imax; i++) {
	if (iuv) {
	  fill_float(rwd_table[i]);
	  fill_float(rws_table[i]);
	  fill_float(rww_table[i]);
	  fill_short(db_table[i]);
	  fill_float(rdv_table[i]);
	  if (fDebug) {
	    std::cout << i + 1 << " " << j + 1 << " " << k + 1 << " " 
		      << rwd_table[i] << " " << rws_table[i] << " " << rww_table[i] << " " 
		      << db_table[i] << " " << rdv_table[i] << std::endl;
	  }
	} else {
	  fill_short(wd_table[i]);
	  fill_short(ws_table[i]);
	  fill_short(ww_table[i]);
	  fill_short(db_table[i]);
	  fill_short(dv_table[i]);
	  if (fDebug) {
	    std::cout << i + 1 << " " << j + 1 << " " << k + 1 << " " 
		      << wd_table[i] << " " << ws_table[i] << " " << ww_table[i] << " " 
		      << db_table[i] << " " << dv_table[i] << std::endl;
	  }
	}
      }
      for(int i = 0; i < imax; i++) {
	if ( ! iuv) {
	  float wdr = wd_table[i] * 0.1;
	  float wsp = ws_table[i] * 0.1;
	  if ( wdr < 0) {
	    (*u_table)(i, j, k) = fillValue;
	    (*v_table)(i, j, k) = fillValue;
	  } else {
	    float xu = - sin(wdr * Pi / 180.0) * wsp;
	    float xv = - cos(wdr * Pi / 180.0) * wsp;
	    (*u_table)(i, j, k) = xu;
	    (*v_table)(i, j, k) = xv;
	  }
	  if (ww_table[i] > -9000) {
	    (*w_table)(i, j, k) = ww_table[i] * 0.01;
	  } else {
	    (*w_table)(i, j, k) = fillValue;
	  }
	  if (dv_table[i] < NoData) {
	    (*div_table)(i, j, k) = dv_table[i] / DivScaleFactor;
	  } else {
	    (*div_table)(i, j, k) = fillValue;
	  }
	} else { // iuv
	  if (rww_table[i] <= - DivScaleFactor) rwd_table[i] = fillValue;
	  if (rws_table[i] <= - DivScaleFactor) rws_table[i] = fillValue;
	  if (rww_table[i] <= - DivScaleFactor) rww_table[i] = fillValue;
	  if (rdv_table[i] <= - DivScaleFactor) rdv_table[i] = fillValue;
	  (*u_table)(i, j, k) = rwd_table[i];
	  (*v_table)(i, j, k) = rws_table[i];
	  (*w_table)(i, j, k) = rww_table[i];
	  (*div_table)(i, j, k) = rdv_table[i];
	}
	if (db_table[i] > -9000) {
	  (*dbz_table)(i, j, k) = db_table[i] * 0.1;
	} else {
	  (*dbz_table)(i, j, k) = fillValue;
	}
	// TODO debug. Fill dbz with known pattern
	// (*dbz_table)(i, j, k) = i * 100 + j * 10 + k;
      }
    }
  }
}

Data::Data(char *fname, bool debug) : file_handle(fname, std::ios::in | std::ios::binary), fDebug(debug)
{
  input_path = fname;
  valid_flag = true;

  // std::ifstream file(fname, std::ios::in | std::ios::binary);
  if (! file_handle.good()) {
    std::cerr << "Unable to read '" << fname << "'" << std::endl;
    valid_flag = false;
    return;
  }

  iuv = false;

  // read past first 4 bytes (Fortran padding)
  fill_char(keyword, 4);

  // read header fields
  
  fill_char(keyword, 4);
  if( (strcmp(keyword, "RUV") == 0) || (strcmp(keyword, "ruv") == 0) )
    iuv = true;

  fill_char(fltname, 8);
  fill_char(stmname, 12);
  fill_char(radar, 4);
  fill_char(experiment, 32);
  fill_char(creattime, 32);
  fill_char(extra1, 28);
  fill_int(imax);
  fill_int(jmax);
  fill_int(kmax);
  fill_int(kount);
  fill_int(nmosm);
  fill_int(iunfld);
  fill_int(iatten);
  fill_int(flag);
  fill_int(extra2);
  fill_int(extra3);

  if(fDebug) {
    std::cout << "imax: " << imax << std::endl;
    std::cout << "jmax: " << jmax << std::endl;
    std::cout << "kmax: " << kmax << std::endl;
    std::cout << "nmosm: " << nmosm << std::endl;
    std::cout << "extra3: " << extra3 << std::endl;
  }
  
  fill_float(stime);
  fill_float(etime);
  fill_float(olat);
  fill_float(olon);
  fill_float(sx);
  fill_float(sy);
  fill_float(sz);
  fill_float(xz);
  fill_float(yz);
  fill_float(zz);
  fill_float(rot);
  fill_float(ra);
  fill_float(co1);
  fill_float(co2);
  fill_float(azmcor);
  fill_float(elcor);
  fill_float(thresh);
  fill_float(powert);
  fill_float(biel);
  fill_float(azbiel);
  fill_float(etime1);
  fill_float(stime2);
  fill_float(extra6);
  fill_float(extra7);

  if (fDebug) {
    std::cout << "olat: " << olat << std::endl;
    std::cout << "olon: " << olon << std::endl;
    std::cout << "extra7: " << extra7 << std::endl;
  }
}

// Write data to a netcdf file

bool Data::write(const char *fname)
{
  netCDF::NcFile out(fname, netCDF::NcFile::replace);

  // Save header values as attribute
  out.putAtt("keyword",    netCDF::ncChar, strlen(keyword), keyword);
  out.putAtt("fltname",    netCDF::ncChar, strlen(fltname), fltname);
  out.putAtt("radar",      netCDF::ncChar, strlen(radar), radar);
  out.putAtt("stmname",    netCDF::ncChar, strlen(stmname), stmname);
  out.putAtt("experiment", netCDF::ncChar, strlen(experiment),experiment);
  out.putAtt("creattime",  netCDF::ncChar, strlen(creattime), creattime);

  out.putAtt("azmcor", netCDF::ncFloat, azmcor);
  out.putAtt("elcor",  netCDF::ncFloat, elcor);
  out.putAtt("thresh", netCDF::ncFloat, thresh);
  out.putAtt("rot",    netCDF::ncFloat, rot);
  out.putAtt("olat",    netCDF::ncFloat, olat);
  out.putAtt("olon",    netCDF::ncFloat, olon);
  out.putAtt("sx",    netCDF::ncFloat, sx);
  out.putAtt("sy",    netCDF::ncFloat, sy);
  out.putAtt("sz",    netCDF::ncFloat, sz);
  out.putAtt("stime2",    netCDF::ncFloat, stime2);
  out.putAtt("etime1",    netCDF::ncFloat, etime1);
  
  // define the dims
  netCDF::NcDim xDim = out.addDim("xDim", imax);
  netCDF::NcDim yDim = out.addDim("yDim", jmax);
  netCDF::NcDim zDim = out.addDim("zDim", kmax);

  std::vector<netCDF::NcDim> dimArray = { zDim, yDim, xDim };
  
  const std::string fillValName = "_FillValue";
  const std::string missValName = "missing_value";
  
  netCDF::NcVar uVar = out.addVar("u", netCDF::ncFloat, dimArray);
  uVar.putAtt(fillValName, netCDF::ncFloat, (float) fillValue);
  uVar.putAtt(missValName, netCDF::ncFloat, (float) fillValue);
  uVar.putVar(u_table->data());

  netCDF::NcVar vVar = out.addVar("v", netCDF::ncFloat, dimArray);
  vVar.putAtt(fillValName, netCDF::ncFloat, (float) fillValue);
  vVar.putAtt(missValName, netCDF::ncFloat, (float) fillValue);
  vVar.putVar(v_table->data());

  netCDF::NcVar wVar = out.addVar("w", netCDF::ncFloat, dimArray);
  wVar.putAtt(fillValName, netCDF::ncFloat, (float) fillValue);
  wVar.putAtt(missValName, netCDF::ncFloat, (float) fillValue);
  wVar.putVar(w_table->data());

  netCDF::NcVar dbzVar = out.addVar("dbz", netCDF::ncFloat, dimArray);
  dbzVar.putAtt(fillValName, netCDF::ncFloat, (float) fillValue);
  dbzVar.putAtt(missValName, netCDF::ncFloat, (float) fillValue);
  dbzVar.putVar(dbz_table->data());

  netCDF::NcVar divVar = out.addVar("div", netCDF::ncFloat, dimArray);
  divVar.putAtt(fillValName, netCDF::ncFloat, (float) fillValue);
  divVar.putAtt(missValName, netCDF::ncFloat, (float) fillValue);
  divVar.putVar(div_table->data());
  
  out.close();
}

void usage(const char *s) {
  std::cout << "Usage: " << s << "[-d] -i <input W file> -o <output NetCDF file>"
	    << std::endl
	    << "\t-d\t\tTurn on debug flag (dump the data that was read)"
    	    << std::endl
    	    << "\t-i <file>\tInput file in the W format"
	    << std::endl
	    << "\t-o <file>\tOutput file in the NetCDF format format"
    	    << std::endl;
}

int main(int argc, char *argv[])
{
  int opt;
  char *in = NULL, *out = NULL;
  bool debug = false;
  
  while( (opt = getopt(argc, argv, "i:o:hd")) != -1)
    switch(opt){
    case 'd':
      debug = true;
      break;
    case 'i':
      in = strdup(optarg);
      break;
    case 'o':
      out = strdup(optarg);
      break;
    case 'h':
    case '?':
      usage(argv[0]);
      exit(0);
    }

  // Make sure we have the right options

  if ( (in == NULL) || (out == NULL) ) {
    usage(argv[0]);
    exit(1);
  }

  Data data(in, debug);
  if ( ! data.valid() )
    exit(1);
  data.ldlnl();
  data.write(out);
  return 0;
}
