#include "crpropa/magneticField/JF12Field_improved.h"
#include "crpropa/Units.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"

#include <iostream>

namespace crpropa {

double JF12Field_improved::logisticFunction(double x, double x0, double w) const{
	return 1. / (1. + exp(-2. * (fabs(x) - x0) / w));
}

JF12Field_improved::JF12Field_improved(bool usedisk, bool usetoroidalhalo, bool usex, double delta, double hx) {
	
	useRegular = true;
	useStriated = false;
	useTurbulent = false;
	
	useDisk = usedisk;
	useToroidalHalo = usetoroidalhalo;
	useX = usex;
	
	// set widths and heigths of field and transition zones
	Delta = delta;
	hX = hx;
	r1 = 5 * kpc;
	r2 = 20 * kpc;
	r1s = r1 + Delta;
	r2s = r2 - Delta;
	
	// set spiral arm parameters
	pitch = 11.5 * M_PI / 180;
	sinPitch = sin(pitch);
	cosPitch = cos(pitch);
	tanPitch = tan(pitch);
	cotPitch =  1. / tanPitch;
	tan90MinusPitch = tan(M_PI / 2 - pitch);

	rArms[0] = 5.1 * kpc;
	rArms[1] = 6.3 * kpc;
	rArms[2] = 7.1 * kpc;
	rArms[3] = 8.3 * kpc;
	rArms[4] = 9.8 * kpc;
	rArms[5] = 11.4 * kpc;
	rArms[6] = 12.7 * kpc;
	rArms[7] = 15.5 * kpc;
	
	
	// set angles of seperating spiral field lines at r1, see paper
	phi0 = 0;
	
	for (int i = 1;i < 9; i++){
		phi0Arms[i] = M_PI - cotPitch * log(rArms[i-1] / r1);
	}
	
	// cyclic closure
	phi0Arms[0] = phi0Arms[8] + 2 * M_PI;
	phi0Arms[9] = phi0Arms[1] - 2 * M_PI;
	phi0Arms[10] = phi0Arms[2] - 2 *M_PI;
	
	// determine index of phi0
	idx0 = 0;
	while (phi0Arms[idx0] > phi0){
		idx0 += 1;
	}
	
	// set regular field parameters
	bRing = 0.1 * muG;
	hDisk = 0.40 * kpc;
	wDisk = 0.27 * kpc;

	bDisk[1] = 0.1 * muG; //called b_1
	bDisk[2] = 3.0 * muG; //b_2
	bDisk[3] = -0.9 * muG;//etc.
	bDisk[4] = -0.8 * muG;
	bDisk[5] = -2.0 * muG;
	bDisk[6] = -4.2 * muG;
	bDisk[7] = 0.0 * muG;
	bDisk[8] = 2.7 * muG;
	
	// re-compute b_8 for flux correction
	double flux1to7 = 0.;
	for (int i = 1; i < 8; i++){
		flux1to7 += (phi0Arms[i-1] - phi0Arms[i]) * bDisk[i];
	}
	bDisk[8] = -flux1to7 / (phi0Arms[7] - phi0Arms[8]);
	bDisk[0] = bDisk[8];
	bDisk[9] = bDisk[1];
	bDisk[10] = bDisk[2];
	
	// set coefficients for phi integration
	
	phiCoeff[0] = 0;
	for (int i = 1; i < 10; i++){
		phiCoeff[i] = phiCoeff[i-1] + (bDisk[i-1] - bDisk[i]) * phi0Arms[i-1];
	}
	
	//correct for H(phi0) = 0
	corr = phiCoeff[idx0] + bDisk[idx0] * phi0;
	for (int i = 1; i < 10; i++){
		phiCoeff[i] = phiCoeff[i] - corr;
	} 
	
	// azimuthal halo parameters
	bNorth = 1.4 * muG;
	bSouth = -1.1 * muG;
	rNorth = 9.22 * kpc;
	rSouth = 17 * kpc;
	wHalo = 0.20 * kpc;
	z0 = 5.3 * kpc;
	
	// X-field parameters
	bX = 4.6 * muG;
	thetaX0 = 49.0 * M_PI / 180;
	sinThetaX0 = sin(thetaX0);
	cosThetaX0 = cos(thetaX0);
	tanThetaX0 = tan(thetaX0);
	cotThetaX0 = 1. / tan(thetaX0);
	rXc = 4.8 * kpc;
	rX = 2.9 * kpc;

	// set striated field parameter
	sqrtbeta = sqrt(1.36);

	// set turbulent field parameters
	bDiskTurb[0] = 10.81 * muG;
	bDiskTurb[1] = 6.96 * muG;
	bDiskTurb[2] = 9.59 * muG;
	bDiskTurb[3] = 6.96 * muG;
	bDiskTurb[4] = 1.96 * muG;
	bDiskTurb[5] = 16.34 * muG;
	bDiskTurb[6] = 37.29 * muG;
	bDiskTurb[7] = 10.35 * muG;

	bDiskTurb5 = 7.63 * muG;
	zDiskTurb = 0.61 * kpc;

	bHaloTurb = 4.68 * muG; 
	rHaloTurb = 10.97 * kpc;
	zHaloTurb = 2.84 * kpc;
	
	
	//debugging output
	//~ std::cout << "phi0Arms:\n";
	//~ for (int i = 0; i < 11; i++){
		//~ std::cout << "Phi[" << i << "] = " << phi0Arms[i] << " rad \n";
	//~ }
	//~ 
	//~ std::cout << "Index of phi0 is: " << idx0 << "\n";
	//~ 
	//~ std::cout << "bDisk:\n";
	//~ for (int i = 0; i < 11; i++){
		//~ std::cout << "bDisk[" << i << "] = " << bDisk[i]/muG << " muG \n";
	//~ }
	//~ 
	//~ std::cout << "corr = " << corr/muG << " muG\n";
	//~ 
	//~ std::cout << "phiCoeff:\n";
	//~ for (int i = 0; i < 10; i++){
		//~ std::cout << "phiCoeff[" << i << "] = " << phiCoeff[i]/muG << " muG \n";
	//~ }
	
	
}

void JF12Field_improved::randomStriated(int seed) {
	useStriated = true;
	int N = 100;
	striatedGrid = new ScalarGrid(Vector3d(0.), N, 0.1 * kpc);

	Random random;
	if (seed != 0)
		random.seed(seed);

	for (int ix = 0; ix < N; ix++)
		for (int iy = 0; iy < N; iy++)
			for (int iz = 0; iz < N; iz++) {
				float &f = striatedGrid->get(ix, iy, iz);
				f = round(random.rand()) * 2 - 1;
			}
}

#ifdef CRPROPA_HAVE_FFTW3F
void JF12Field_improved::randomTurbulent(int seed) {
	useTurbulent = true;
	// turbulent field with Kolmogorov spectrum, B_rms = 1 and Lc = 60 parsec
	turbulentGrid = new VectorGrid(Vector3d(0.), 256, 4 * parsec);
	initTurbulence(turbulentGrid, 1, 8 * parsec, 272 * parsec, -11./3., seed);
}
#endif

void JF12Field_improved::setStriatedGrid(ref_ptr<ScalarGrid> grid) {
	useStriated = true;
	striatedGrid = grid;
}

void JF12Field_improved::setTurbulentGrid(ref_ptr<VectorGrid> grid) {
	useTurbulent = true;
	turbulentGrid = grid;
}

ref_ptr<ScalarGrid> JF12Field_improved::getStriatedGrid() {
	return striatedGrid;
}

ref_ptr<VectorGrid> JF12Field_improved::getTurbulentGrid() {
	return turbulentGrid;
}

void JF12Field_improved::setUseRegular(bool use) {
	useRegular = use;
}

void JF12Field_improved::setUseStriated(bool use) {
	if ((use) and (striatedGrid)) {
		std::cout << "JF12Field_improved: No striated field set: ignored" << std::endl;
		return;
	}
	useStriated = use;
}

void JF12Field_improved::setUseTurbulent(bool use) {
	if ((use) and (turbulentGrid)) {
		std::cout << "JF12Field_improved: No turbulent field set: ignored" << std::endl;
		return;
	}
	useTurbulent = use;
}

bool JF12Field_improved::isUsingRegular() {
	return useRegular;
}

bool JF12Field_improved::isUsingStriated() {
	return useStriated;
}

bool JF12Field_improved::isUsingTurbulent() {
	return useTurbulent;
}

Vector3d JF12Field_improved::getRegularField(const Vector3d& pos) const {
	Vector3d b(0.);

	double r = sqrt(pos.x * pos.x + pos.y * pos.y); // in-plane radius
	double d = pos.getR(); // distance to galactic center

	double phi = pos.getPhi(); // azimuth
	double sinPhi = sin(phi);
	double cosPhi = cos(phi);
	
	////////////////////////////////////////////////////////////////
	
	// disk field
	b += getDiskField(r, pos.z, phi, sinPhi, cosPhi);
	
	
	////////////////////////////////////////////////////////////////
	
	// toroidal halo field
	b += getToroidalHaloField(r, pos.z, sinPhi, cosPhi);
	
	////////////////////////////////////////////////////////////////
	
	// poloidal halo field
	
	b += getXField(r, pos.z, sinPhi, cosPhi);
	
	////////////////////////////////////////////////////////////////

	return b;
}

Vector3d JF12Field_improved::getStriatedField(const Vector3d& pos) const {
	return (getRegularField(pos)
			* (1. + sqrtbeta * striatedGrid->closestValue(pos)));
}

double JF12Field_improved::getTurbulentStrength(const Vector3d& pos) const {
	if (pos.getR() > 20 * kpc)
		return 0;

	double r = sqrt(pos.x * pos.x + pos.y * pos.y); // in-plane radius
	double phi = pos.getPhi(); // azimuth

	// disk
	double bDisk = 0;
	if (r < 5 * kpc) {
		bDisk = bDiskTurb5;
	} else {
		// spiral region
		double r_negx = r * exp(-(phi - M_PI) / tan90MinusPitch);
		if (r_negx > rArms[7])
			r_negx = r * exp(-(phi + M_PI) / tan90MinusPitch);
		if (r_negx > rArms[7])
			r_negx = r * exp(-(phi + 3 * M_PI) / tan90MinusPitch);

		for (int i = 7; i >= 0; i--)
			if (r_negx < rArms[i])
				bDisk = bDiskTurb[i];

		bDisk *= (5 * kpc) / r;
	}
	bDisk *= exp(-0.5 * (pos.z / zDiskTurb) * (pos.z / zDiskTurb));

	// halo
	double bHalo = bHaloTurb * exp(-r / rHaloTurb)
			* exp(-0.5 * (pos.z / zHaloTurb) * (pos.z / zHaloTurb));

	// modulate turbulent field
	return sqrt(bDisk * bDisk + bHalo * bHalo);
}

Vector3d JF12Field_improved::getTurbulentField(const Vector3d& pos) const {
	return (turbulentGrid->interpolate(pos) * getTurbulentStrength(pos));
}

Vector3d JF12Field_improved::getField(const Vector3d& pos) const {
	Vector3d b(0.);
	if (useTurbulent)
		b += getTurbulentField(pos);
	if (useStriated)
		b += getStriatedField(pos);
	else if (useRegular)
		b += getRegularField(pos);
	return b;
}

//improved disk field with transition to ring region, delta1 and delta2 are the transition widths
Vector3d JF12Field_improved::getDiskField(const double r, const double z, const double phi, const double sinPhi, const double cosPhi) const {
	
	Vector3d b(0.);

	if (!useDisk){
		return b;
	}
	
	else {
		
		double lfDisk = logisticFunction(z, hDisk, wDisk);
		
		double hint = PhiIntegralH(r, phi);
		double mag1 = getSpiralStrength(r, phi);
		
		if ((r1 < r) && (r < r2)) {
			
			double pdelta = p(r);
			double qdelta = q(r);
			double br = pdelta * mag1 * sinPitch;
			double bphi = pdelta * mag1 * cosPitch - qdelta * hint * sinPitch;
		
			b.x += br * cosPhi - bphi * sinPhi;
			b.y += br * sinPhi + bphi * cosPhi;
			
			b *= (1 - lfDisk);
		}
		
		return b;
	}		
}

//same as initial JF12, moved in a seperate function for clearness
Vector3d JF12Field_improved::getToroidalHaloField(const double r, const double z, const double sinPhi, const double cosPhi) const {
	
	Vector3d b(0.);
	
	if (!useToroidalHalo){
		return b;
	}
	
	if ((r*r + z * z <=1 * kpc) || (r*r + z * z >= 400 * kpc * kpc)){
		return b;
	}
	
	else{
		
		double lfDisk = logisticFunction(z, hDisk, wDisk);
		double bMagH = exp(-fabs(z) / z0) * lfDisk;
		
		if (z >= 0)
			bMagH *= bNorth * (1 - logisticFunction(r, rNorth, wHalo));
		else
			bMagH *= bSouth * (1 - logisticFunction(r, rSouth, wHalo));
		b.x += -bMagH * sinPhi;
		b.y += bMagH * cosPhi;
		return b;
	}
}

//improved X-field with parabolic field lines at abs(z) < hX
Vector3d JF12Field_improved::getXField(const double r, const double z, const double sinPhi, const double cosPhi) const {
	
	
	Vector3d b(0.);
	
	if (!useX || (r * r + z * z >= 400 * kpc * kpc)){
		return b;
	}
	
	if (r == 0){
		b.z = bX / ((1 + abs(z) * cotThetaX0 / rXc) * (1 + abs(z) * cotThetaX0 / rXc));
	}
	
	else{
		
		double bMagX;
		double sinThetaX, cosThetaX;
		double rp; //radius where current intial field line passes z = 0
		double rc = rXc + fabs(z) / tanThetaX0; 
		double r0c = rXc + hX / tanThetaX0; //radius where field line through rXc passes z = hX
		double f, r0, br0, bz0;
		bool inner = true;//distinguish between inner and outer region
		
		// return intial field if z>=hX
		if (fabs(z) >= hX){
			
			if (r < rc) {	
			// inner varying elevation region
			
				rp = r * rXc / rc;
				bMagX = bX * exp(-1 * rp / rX) * (rXc / rc) * (rXc / rc);
				
				double thetaX = atan(fabs(z) / (r - rp));
				
				if (z == 0)
					thetaX = M_PI / 2.;
				
				sinThetaX = sin(thetaX);
				cosThetaX = cos(thetaX);
				
				double zsign = z < 0 ? -1 : 1;
				
				b.x += zsign * bMagX * cosThetaX * cosPhi;
				b.y += zsign * bMagX * cosThetaX * sinPhi;
				b.z += bMagX * sinThetaX;

				return b;
			}
			
			else {
			//outer constant elevation region
				rp = r - fabs(z) / tanThetaX0;
				bMagX = bX * exp(-rp / rX) * (rp / r);
				
				sinThetaX = sinThetaX0;
				cosThetaX = cosThetaX0;
				
				double zsign = z < 0 ? -1 : 1;
				b.x += zsign * bMagX * cosThetaX * cosPhi;
				b.y += zsign * bMagX * cosThetaX * sinPhi;
				b.z += bMagX * sinThetaX;
				
				return b;
			}
		}
		//parabolic field lines for z<hX
		else {
				//determine r at which parabolic field line through (r,z) passes z = hX
				r0 = r * 1. / (1.- 1./ (2. * (hX + rXc * tanThetaX0)) * (hX - z * z / hX));
				
				//determine correct region (inner/outer)
				if (r0 >= r0c){
					r0 = r + 1. / (2. * tanThetaX0) * (hX - z * z / hX);
					inner = false;
				}
				
				//field strength at that position
				 if (r0 < r0c){
					 
					 rp = r0 * rXc / r0c;
					 double thetaX = atan(hX / (r0 - rp));
					 
					 //field strength at (r0,hX) for inner region
					 br0 = bX * exp(- rp / rX) * (rXc/ r0c) * (rXc/ r0c) * cos(thetaX);
					 bz0 = bX * exp(- rp / rX) * (rXc/ r0c) * (rXc/ r0c) * sin(thetaX);
				 }
				 else {
					 
					 //field strength at (r0,hX) for outer region
					 rp = r0 - hX / tanThetaX0;
					 br0 =  bX * exp(- rp / rX) * (rp/r0) * cosThetaX0;
					 bz0 =  bX * exp(- rp / rX) * (rp/r0) * sinThetaX0;
					 
				 }
				 
				 //compute factor F for solenoidality
				 if (inner){
					 f = 1. / ((1. - 1./( 2. + 2. * (rXc * tanThetaX0/ hX)) * (1. - (z / hX) * (z / hX))) * (1. - 1./( 2. + 2. * (rXc * tanThetaX0/ hX)) * (1. - (z / hX) * (z / hX))));
				 }
				 else{
					 f = 1. + 1/ (2 * r * tanThetaX0/ hX) * (1. - (z / hX) * (z / hX)); 
				 }
				 
				 
				 double br = z / hX * f * br0;
				 double bz = bz0 * f;
				 
				 b.x += br * cosPhi;
				 b.y += br * sinPhi;
				 b.z += bz;
				 
				 return b;
		}
	}
}

//transition polynomial p_delta(r)
double JF12Field_improved::p(const double r) const {
	
	//0 disk field outside
	if ((r < r1) || (r > r2)) {
		return 0.;
	}	
	//unchanged field
	if ((r > r1s) && (r < r2s)) {
		return r1/r;
	}
	
	//transitions region parameters
	double r_a = r1;
	double r_s = r1s;
	
	if (r >= r2s) {
		r_a = r2;
		r_s = r2s;
	}
	
	double fakt = (r_a / r_s - 2.) / ((r_a - r_s) *  (r_a - r_s));
	return (r1/r_s) * (2. - r / r_s + fakt * (r-r_s) * (r-r_s));
}	

//transition polynomial derivative
double JF12Field_improved::q(const double r) const {
	
	//0 disk field outside
	if ((r < r1) || (r > r2)) {
		return 0.;
	}	
	//unchanged field
	if ((r > r1s) && (r < r2s)) {
		return 0.;
	}
	
	//transitions region parameters
	double r_a = r1;
	double r_s = r1s;
	
	if (r >= r2s) {
		r_a = r2;
		r_s = r2s;
	}
	
	double fakt = (r_a / r_s - 2.) / ((r_a - r_s) * (r_a - r_s));
	return (r1/r_s) * (2. - 2. * r/r_s + fakt * (3. * r * r - 4. * r * r_s + r_s * r_s));
}

//evaluate BphiIntegral 
double JF12Field_improved::PhiIntegralH(const double r, const double phi) const {
	
	double H_ret = 0.;
	int idx = 1;
	
	if ((r1 < r) && (r < r2)){
		double phi1 = phi - log(r/r1) * cotPitch;
		phi1 = atan2(sin(phi1) , cos(phi1)); // map to [-pi,+pi]
		while (phi1 < phi0Arms[idx]){
			idx += 1;
		}
		H_ret = phi1 * bDisk[idx] + phiCoeff[idx];
	}
	
	return H_ret;
	
}

//return correct field strength b_j in spiral arm
double JF12Field_improved::getSpiralStrength(const double r, const double phi) const {
	
	double b_ret = 0.;
	int idx = 1;
	
	if ((r1 < r) && (r < r2)){
		double phi1 = phi - log(r/r1) * cotPitch;
		phi1 = atan2(sin(phi1), cos(phi1)); // map to [-pi,+pi]
		while (phi1 < phi0Arms[idx]){
			idx += 1;
		}
		b_ret = bDisk[idx];
	}
	
	return b_ret;
}

} // namespace crpropa
