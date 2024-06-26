#include "crpropa/GridTools.h"
#include "crpropa/magneticField/MagneticField.h"

#include <fstream>
#include <sstream>

namespace crpropa {

void scaleGrid(ref_ptr<Grid1f> grid, double a) {
	for (int ix = 0; ix < grid->getNx(); ix++)
		for (int iy = 0; iy < grid->getNy(); iy++)
			for (int iz = 0; iz < grid->getNz(); iz++)
				grid->get(ix, iy, iz) *= a;
}

void scaleGrid(ref_ptr<Grid3f> grid, double a) {
	for (int ix = 0; ix < grid->getNx(); ix++)
		for (int iy = 0; iy < grid->getNy(); iy++)
			for (int iz = 0; iz < grid->getNz(); iz++)
				grid->get(ix, iy, iz) *= a;
}

Vector3f meanFieldVector(ref_ptr<Grid3f> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	Vector3f mean(0.);
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				mean += grid->get(ix, iy, iz);
	return mean / Nx / Ny / Nz;
}

double meanFieldStrength(ref_ptr<Grid3f> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	double mean = 0;
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				mean += grid->get(ix, iy, iz).getR();
	return mean / Nx / Ny / Nz;
}

double meanFieldStrength(ref_ptr<Grid1f> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	double mean = 0;
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				mean += grid->get(ix, iy, iz);
	return mean / Nx / Ny / Nz;
}

double rmsFieldStrength(ref_ptr<Grid3f> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	double sumV2 = 0;
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				sumV2 += grid->get(ix, iy, iz).getR2();
	return std::sqrt(sumV2 / Nx / Ny / Nz);
}

double rmsFieldStrength(ref_ptr<Grid1f> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	double sumV2 = 0;
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				sumV2 += pow(grid->get(ix, iy, iz), 2);
	return std::sqrt(sumV2 / Nx / Ny / Nz);
}

std::array<float, 3> rmsFieldStrengthPerAxis(ref_ptr<Grid3f> grid) {
    size_t Nx = grid->getNx();
    size_t Ny = grid->getNy();
    size_t Nz = grid->getNz();
    float sumV2_x = 0;
    float sumV2_y = 0;
    float sumV2_z = 0;
    for (int ix = 0; ix < Nx; ix++)
        for (int iy = 0; iy < Ny; iy++)
            for (int iz = 0; iz < Nz; iz++) {
                sumV2_x += pow(grid->get(ix, iy, iz).x, 2);
                sumV2_y += pow(grid->get(ix, iy, iz).y, 2);
                sumV2_z += pow(grid->get(ix, iy, iz).z, 2);
            }
    return {
        std::sqrt(sumV2_x / Nx / Ny / Nz),
        std::sqrt(sumV2_y / Nx / Ny / Nz),
        std::sqrt(sumV2_z / Nx / Ny / Nz)
    };
}

void fromMagneticField(ref_ptr<Grid3f> grid, ref_ptr<MagneticField> field) {
	Vector3d origin = grid->getOrigin();
	Vector3d spacing = grid->getSpacing();
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	for (size_t ix = 0; ix < Nx; ix++)
		for (size_t iy = 0; iy < Ny; iy++)
			for (size_t iz = 0; iz < Nz; iz++) {
				Vector3d pos = Vector3d(double(ix) + 0.5, double(iy) + 0.5, double(iz) + 0.5) * spacing + origin;
				Vector3d B = field->getField(pos);
				grid->get(ix, iy, iz) = B;
	}
}

void fromMagneticFieldStrength(ref_ptr<Grid1f> grid, ref_ptr<MagneticField> field) {
	Vector3d origin = grid->getOrigin();
	Vector3d spacing = grid->getSpacing();
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	for (size_t ix = 0; ix < Nx; ix++)
		for (size_t iy = 0; iy < Ny; iy++)
			for (size_t iz = 0; iz < Nz; iz++) {
				Vector3d pos = Vector3d(double(ix) + 0.5, double(iy) + 0.5, double(iz) + 0.5) * spacing + origin;
				double s = field->getField(pos).getR();
				grid->get(ix, iy, iz) = s;
	}
}

void loadGrid(ref_ptr<Grid3f> grid, std::string filename, double c) {
	std::ifstream fin(filename.c_str(), std::ios::binary);
	if (!fin) {
		std::stringstream ss;
		ss << "load Grid3f: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}

	// get length of file and compare to size of grid
	fin.seekg(0, fin.end);
	size_t length = fin.tellg() / sizeof(float);
	fin.seekg (0, fin.beg);

	size_t nx = grid->getNx();
	size_t ny = grid->getNy();
	size_t nz = grid->getNz();

	if (length != (3 * nx * ny * nz))
		throw std::runtime_error("loadGrid: file and grid size do not match");

	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				Vector3f &b = grid->get(ix, iy, iz);
				fin.read((char*) &(b.x), sizeof(float));
				fin.read((char*) &(b.y), sizeof(float));
				fin.read((char*) &(b.z), sizeof(float));
				b *= c;
			}
		}
	}
	fin.close();
}

void loadGrid(ref_ptr<Grid1f> grid, std::string filename, double c) {
	std::ifstream fin(filename.c_str(), std::ios::binary);
	if (!fin) {
		std::stringstream ss;
		ss << "load Grid1f: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}

	// get length of file and compare to size of grid
	fin.seekg(0, fin.end);
	size_t length = fin.tellg() / sizeof(float);
	fin.seekg (0, fin.beg);

	size_t nx = grid->getNx();
	size_t ny = grid->getNy();
	size_t nz = grid->getNz();

	if (length != (nx * ny * nz))
		throw std::runtime_error("loadGrid: file and grid size do not match");

	for (int ix = 0; ix < nx; ix++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int iz = 0; iz < nz; iz++) {
				float &b = grid->get(ix, iy, iz);
				fin.read((char*) &b, sizeof(float));
				b *= c;
			}
		}
	}
	fin.close();
}

void dumpGrid(ref_ptr<Grid3f> grid, std::string filename, double c) {
	std::ofstream fout(filename.c_str(), std::ios::binary);
	if (!fout) {
		std::stringstream ss;
		ss << "dump Grid3f: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				Vector3f b = grid->get(ix, iy, iz) * c;
				fout.write((char*) &(b.x), sizeof(float));
				fout.write((char*) &(b.y), sizeof(float));
				fout.write((char*) &(b.z), sizeof(float));
			}
		}
	}
	fout.close();
}

void dumpGrid(ref_ptr<Grid1f> grid, std::string filename, double c) {
	std::ofstream fout(filename.c_str(), std::ios::binary);
	if (!fout) {
		std::stringstream ss;
		ss << "dump Grid1f: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				float b = grid->get(ix, iy, iz) * c;
				fout.write((char*) &b, sizeof(float));
			}
		}
	}
	fout.close();
}

void loadGridFromTxt(ref_ptr<Grid3f> grid, std::string filename, double c) {
	std::ifstream fin(filename.c_str());
	if (!fin) {
		std::stringstream ss;
		ss << "load Grid3f: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	// skip header lines
	while (fin.peek() == '#')
		fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				Vector3f &b = grid->get(ix, iy, iz);
				fin >> b.x >> b.y >> b.z;
				b *= c;
				if (fin.eof())
					throw std::runtime_error("load Grid3f: file too short");
			}
		}
	}
	fin.close();
}

ref_ptr<Grid3f> loadGrid3fFromTxt(std::string filename, double c) {
	std::ifstream fin(filename.c_str());
	if (!fin) {
		std::stringstream ss;
		ss << "load Grid3f: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}

	// search in header lines for GridProperties
	while (fin.peek() == '#') {
		std::string line;
		std::getline(fin, line);

		// find gridproperties in the header 
		if (line.rfind("GridProperties:") == 2) {	
			GridProperties gp(Vector3d(0.), 1, 1, 1, 1.); // simple grid properties for default
			std::stringstream ss(line); 

			// skip first names and check type 
			std::string name, type;
			ss >> name >> name >> name >> type;
			if (type != "Grid3f") 
				throw std::runtime_error("Tried to load Grid3f, but Gridproperties assume grid type " + type);

			// grid origin
			double x, y, z;
			ss >> name >> x >> y >> z ; 
			gp.origin = Vector3d(x, y, z);

			// grid size
			ss >> name >> gp.Nx >> gp.Ny >> gp.Nz;

			// spacing
			double dX, dY, dZ;
			ss >> name >> dX >> dY >> dZ;
			gp.spacing = Vector3d(dX, dY, dZ);

			// reflective
			ss >> name >> gp.reflective;
			
			// clip volume
			bool clip; 
			ss >> name >> clip;
			gp.setClipVolume(clip);

			// interpolation type 
			ss >> name >> type;
			if (type == "TRICUBIC")
				gp.setInterpolationType(TRICUBIC);
			else if (type == "NEAREST_NEIGHBOUR")
				gp.setInterpolationType(NEAREST_NEIGHBOUR);
			else 
				gp.setInterpolationType(TRILINEAR);

			// create new grid
			ref_ptr<Grid3f> grid = new Grid3f(gp);
			fin.close();

			// load data for grid
			loadGridFromTxt(grid, filename, c); 

			return grid;
		}
	}
	throw std::runtime_error("could not find GridProperties in file " + filename);
}


void loadGridFromTxt(ref_ptr<Grid1f> grid, std::string filename, double c) {
	std::ifstream fin(filename.c_str());
	if (!fin) {
		std::stringstream ss;
		ss << "load Grid1f: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	
	// skip header lines
	while (fin.peek() == '#') 
		fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				float &b = grid->get(ix, iy, iz);
				fin >> b;
				b *= c;
				if (fin.eof())
					throw std::runtime_error("load Grid1f: file too short");
			}
		}
	}
	fin.close();
}

ref_ptr<Grid1f> loadGrid1fFromTxt(std::string filename, double c) {
	std::ifstream fin(filename.c_str());
	if (!fin) {
		std::stringstream ss;
		ss << "load Grid1f: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}

	// search in header lines for GridProperties
	while (fin.peek() == '#') {
		std::string line;
		std::getline(fin, line);

		// find gridproperties in the header 
		if (line.rfind("GridProperties:") == 2) {	
			GridProperties gp(Vector3d(0.), 1, 1, 1, 1.); // simple grid properties for default
			std::stringstream ss(line); 

			// skip first names and check type 
			std::string name, type;
			ss >> name >> name >> name >> type;
			if (type != "Grid1f") 
				throw std::runtime_error("Tried to load Grid1f, but Gridproperties assume grid type " + type);

			// grid origin
			double x, y, z;
			ss >> name >> x >> y >> z ; 
			gp.origin = Vector3d(x, y, z);

			// grid size
			ss >> name >> gp.Nx >> gp.Ny >> gp.Nz;

			// spacing
			double dX, dY, dZ;
			ss >> name >> dX >> dY >> dZ;
			gp.spacing = Vector3d(dX, dY, dZ);

			// reflective
			ss >> name >> gp.reflective;

			// clip volume
			bool clip; 
			ss >> name >> clip;
			gp.setClipVolume(clip);

			// interpolation type 
			ss >> name >> type;
			if (type == "TRICUBIC")
				gp.setInterpolationType(TRICUBIC);
			else if (type == "NEAREST_NEIGHBOUR")
				gp.setInterpolationType(NEAREST_NEIGHBOUR);
			else 
				gp.setInterpolationType(TRILINEAR);

			// create new grid
			ref_ptr<Grid1f> grid = new Grid1f(gp);
			fin.close();

			// load data for grid
			loadGridFromTxt(grid, filename, c); 

			return grid;
		}
	}
	throw std::runtime_error("could not find GridProperties in file " + filename);
}

void dumpGridToTxt(ref_ptr<Grid3f> grid, std::string filename, double c, bool saveProp) {
	std::ofstream fout(filename.c_str());
	if (!fout) {
		std::stringstream ss;
		ss << "dump Grid3f: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}

	// store the properties in the file as header information
	if (saveProp) {
		fout << "# GridProperties: Type Grid3f" 
			<< "\t" << "origin: " << grid -> getOrigin()
			<< "\t" << "gridsize: " << grid -> getNx() << " " << grid -> getNy() << " " << grid -> getNz()
			<< "\t" << "spacing: " << grid -> getSpacing ()
			<< "\t" << "reflective: " << grid -> isReflective()
			<< "\t" << "clipVolume: " << grid -> getClipVolume()
			<< "\t" << "interpolation: " << grid -> getInterpolationTypeName() << "\n";
	}

	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				Vector3f b = grid->get(ix, iy, iz) * c;
				fout << b << "\n";
			}
		}
	}
	fout.close();
}

void dumpGridToTxt(ref_ptr<Grid1f> grid, std::string filename, double c, bool saveProp) {
	std::ofstream fout(filename.c_str());
	if (!fout) {
		std::stringstream ss;
		ss << "dump Grid1f: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}

	// save properties as header information 
	if (saveProp) {
		fout << "# GridProperties: Type Grid1f" 
			<< "\t" << "origin: " << grid -> getOrigin()
			<< "\t" << "gridsize: " << grid -> getNx() << " " << grid -> getNy() << " " << grid -> getNz()
			<< "\t" << "spacing: " << grid -> getSpacing ()
			<< "\t" << "reflective: " << grid -> isReflective()
			<< "\t" << "clipVolume: " << grid -> getClipVolume()
			<< "\t" << "interpolation: " << grid -> getInterpolationTypeName() << "\n";
	}

	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				float b = grid->get(ix, iy, iz) * c;
				fout << b << "\n";
			}
		}
	}
	fout.close();
}

#ifdef CRPROPA_HAVE_FFTW3F

std::vector<std::pair<int, float>> gridPowerSpectrum(ref_ptr<Grid3f> grid) {

  double rms = rmsFieldStrength(grid);
  size_t n = grid->getNx(); // size of array

  // arrays to hold the complex vector components of the B(k)-field
  fftwf_complex *Bkx, *Bky, *Bkz;
  Bkx = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * n * n * n);
  Bky = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * n * n * n);
  Bkz = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * n * n * n);

  fftwf_complex *Bx = (fftwf_complex *)Bkx;
  fftwf_complex *By = (fftwf_complex *)Bky;
  fftwf_complex *Bz = (fftwf_complex *)Bkz;

  // save to temp
  int i;
  for (size_t ix = 0; ix < n; ix++) {
    for (size_t iy = 0; iy < n; iy++) {
      for (size_t iz = 0; iz < n; iz++) {
        i = ix * n * n + iy * n + iz;
        Vector3<float> &b = grid->get(ix, iy, iz);
        Bx[i][0] = b.x / rms;
        By[i][0] = b.y / rms;
        Bz[i][0] = b.z / rms;
      }
    }
  }

  // in-place, real to complex, inverse Fourier transformation on each component
  // note that the last elements of B(x) are unused now
  fftwf_plan plan_x =
      fftwf_plan_dft_3d(n, n, n, Bx, Bkx, FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_execute(plan_x);
  fftwf_destroy_plan(plan_x);

  fftwf_plan plan_y =
      fftwf_plan_dft_3d(n, n, n, By, Bky, FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_execute(plan_y);
  fftwf_destroy_plan(plan_y);

  fftwf_plan plan_z =
      fftwf_plan_dft_3d(n, n, n, Bz, Bkz, FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_execute(plan_z);
  fftwf_destroy_plan(plan_z);

  float power;
  std::map<size_t, std::pair<float, int>> spectrum;
  int k;

  for (size_t ix = 0; ix < n; ix++) {
    for (size_t iy = 0; iy < n; iy++) {
      for (size_t iz = 0; iz < n; iz++) {
        i = ix * n * n + iy * n + iz;
        k = static_cast<int>(
            std::floor(std::sqrt(ix * ix + iy * iy + iz * iz)));
        if (k > n / 2. || k == 0)
          continue;
        power = ((Bkx[i][0] * Bkx[i][0] + Bkx[i][1] * Bkx[i][1]) +
                 (Bky[i][0] * Bky[i][0] + Bky[i][1] * Bky[i][1]) +
                 (Bkz[i][0] * Bkz[i][0] + Bkz[i][1] * Bkz[i][1]));
        if (spectrum.find(k) == spectrum.end()) {
          spectrum[k].first = power;
          spectrum[k].second = 1;
        } else {
          spectrum[k].first += power;
          spectrum[k].second += 1;
        }
      }
    }
  }

  fftwf_free(Bkx);
  fftwf_free(Bky);
  fftwf_free(Bkz);

  std::vector<std::pair<int, float>> points;
  for (std::map<size_t, std::pair<float, int>>::iterator it = spectrum.begin();
       it != spectrum.end(); ++it) {
    points.push_back(
        std::make_pair(it->first, (it->second).first / (it->second).second));
  }

  return points;
}

#endif // CRPROPA_HAVE_FFTW3F


} // namespace crpropa
