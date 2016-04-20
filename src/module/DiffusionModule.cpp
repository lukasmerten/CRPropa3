#include "crpropa/module/DiffusionModule.h"

#include <crpropa/Random.h>

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <stdexcept>


using namespace crpropa;

// Defining Cash-Karp coefficients
const double a[] = { 0., 0., 0., 0., 0., 0., 1. / 5., 0., 0., 0., 0.,
		0., 3. / 40., 9. / 40., 0., 0., 0., 0., 3. / 10., -9. / 10., 6. / 5.,
		0., 0., 0., -11. / 54., 5. / 2., -70. / 27., 35. / 27., 0., 0., 1631.
				/ 55296., 175. / 512., 575. / 13824., 44275. / 110592., 253.
				/ 4096., 0. };

const double b[] = { 37. / 378., 0, 250. / 621., 125. / 594., 0., 512.
		/ 1771. };

const double bs[] = { 2825. / 27648., 0., 18575. / 48384., 13525.
		/ 55296., 277. / 14336., 1. / 4. };



DiffusionModule::DiffusionModule(ref_ptr<MagneticField> field, double tolerance, 
				 double minStep, double maxStep, double kappa) :
  field(field)
{
  setMaximumStep(maxStep);
  setMinimumStep(minStep);
  setTolerance(tolerance);
  setKappaN(kappa);
  setKappaB(kappa);
  setScale(1.);
  setAlpha(1./3.);

}


void DiffusionModule::process(Candidate *candidate) const {
	// save the new previous particle state
	ParticleState &current = candidate->current;
	candidate->previous = current;
	
	double step = clip(candidate->getNextStep(), minStep, maxStep);

	// rectilinear propagation for neutral particles
	if (current.getCharge() == 0) {
		Vector3d dir = current.getDirection();
		current.setPosition(current.getPosition() + dir * step);
		candidate->setCurrentStep(step);
		candidate->setNextStep(maxStep);
		return;
	}

	Vector3d PosIn = current.getPosition();
	double z = candidate->getRedshift();
	double rig = current.getEnergy() / current.getCharge();

	double DifCoeff = scale * 6.1e24 * pow((std::abs(rig) / 4.0e9), alpha);
	double BTensor[] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	BTensor[0] = pow( 2  * DifCoeff, 0.5);
	BTensor[4] = pow(2 * kappaN * DifCoeff, 0.5);
	BTensor[8] = pow(2 * kappaB * DifCoeff, 0.5);
	
	double eta[] = {0., 0., 0.};
	for(size_t i=0; i < 3; i++) {
	  eta[i] =  Random::instance().randNorm();
	}

	double TStep = BTensor[0] * eta[0];
	double NStep = BTensor[4] * eta[1];
	double BStep = BTensor[8] * eta[2];

	double h = step / c_light;
	double hTry, r;

	Vector3d TVec(0.);
	Vector3d NVec(0.);
	Vector3d BVec(0.);
	Vector3d PosOut = Vector3d(0.);
	Vector3d DirOut = Vector3d(0.);
	Vector3d PosErr = Vector3d(0.);
	 

	do {
	  hTry = h;
	  double propStep =  TStep * pow(hTry, 0.5);
	  
	  tryStep(PosIn, PosOut, PosErr, TVec, NVec, BVec, z, propStep);
	  
	  // calculate the relative position error r and the next time step h
	  r = PosErr.getR() / std::abs(propStep) / tolerance;
	  h *= 0.95 * pow(r, -0.2);
	  // prevent h from too strong variations
	  h = clip(h, 0.1 * hTry, 5 * hTry);
	} while (r > 1 && hTry > minStep / c_light && TVec.getR()==TVec.getR());

	if (TVec.getR() != TVec.getR()) {
	  Vector3d dir = current.getDirection();
	  current.setPosition(current.getPosition() + dir * step);
	  candidate->setCurrentStep(step);
	  candidate->setNextStep(step);
	  return;
	}

	Vector3d PO = PosIn + (TVec * std::abs(TStep) + NVec * NStep + BVec * BStep) * pow(hTry, 0.5);
	

	DirOut = (PO -PosIn).getUnitVector();
	current.setPosition(PO);
	current.setDirection(DirOut);
	candidate->setCurrentStep(hTry * c_light);
	candidate->setNextStep(h * c_light);
}


void DiffusionModule::tryStep(const Vector3d &PosIn, Vector3d &POut, Vector3d &PosErr, Vector3d &TVec,Vector3d &NVec,Vector3d &BVec,double z, double propStep) const {

	Vector3d k[] = {Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.)};
	Vector3d Pos[] = {Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.)};
	std::cout << "propStep " << propStep << "\n";
	for (size_t i = 0; i < 6; i++) {
		Vector3d y_n = PosIn;
		for (size_t j = 0; j < i; j++)
		  y_n += k[j] * a[i * 6 + j] * propStep;
		Pos[i] = y_n;
		
		// update k_i
		Vector3d BField(0.);
		try {
		  BField = field->getField(y_n, z);
		} catch (std::exception &e) {
		  std::cerr << "PropagationCK: Exception in getField." << std::endl;
		  std::cerr << e.what() << std::endl;
		}
		
		k[i] = BField.getUnitVector();
		POut += k[i] * b[i] * propStep;
		PosErr +=  (k[i] * (b[i] - bs[i])) * propStep;
		
		// Calculate the Normal-vector as the derivative of the Tangent-vector
		if (i > 0){
		  NVec += ( k[i] - k[i-1] ) / ( Pos[i-1] - Pos[i] ).getR();
		}

	}
	TVec = POut.getUnitVector();
	
	// Choose a random perpendicular vector as the Normal-vector in case of non-curved field line
	if (NVec.getR() == 0.){
	  NVec = TVec.cross( Random::instance().randVector() );
	    }
	
	NVec = NVec.getUnitVector();
	BVec = (TVec.cross(NVec)).getUnitVector();
	std::cout << "Tripod " << TVec <<", " << NVec <<", " << BVec <<"\n";
}


void DiffusionModule::setMinimumStep(double min) {
	if (min < 0)
		throw std::runtime_error("DiffusionModule: minStep < 0 ");
	if (min > maxStep)
		throw std::runtime_error("DiffusionModule: minStep > maxStep");
	minStep = min;
}

void DiffusionModule::setMaximumStep(double max) {
	if (max < minStep)
		throw std::runtime_error("DiffusionModule: maxStep < minStep");
	maxStep = max;
}

void DiffusionModule::setTolerance(double tol) {
	if ((tol > 1) or (tol < 0))
		throw std::runtime_error(
				"DiffusionModule: tolerance error not in range 0-1");
	tolerance = tol;
}

void DiffusionModule::setKappaN(double kap) {
	if ((kap > 1) or (kap < 0))
		throw std::runtime_error(
				"DiffusionModule: kappaN not in range 0-1");
	kappaN = kap;
}

void DiffusionModule::setKappaB(double kap) {
	if ((kap > 1) or (kap < 0))
		throw std::runtime_error(
				"DiffusionModule: kappaB not in range 0-1");
	kappaB = kap;
}

void DiffusionModule::setAlpha(double a) {
	if ((a > 2.) or (a < 0))
		throw std::runtime_error(
				"DiffusionModule: alpha not in range 0-2");
	alpha = a;
}

void DiffusionModule::setScale(double s) {
	if (s < 0)
		throw std::runtime_error(
				"DiffusionModule: Scale error: Scale < 0");
	scale = s;
}


double DiffusionModule::getMinimumStep() const {
	return minStep;
}

double DiffusionModule::getMaximumStep() const {
	return maxStep;
}

double DiffusionModule::getTolerance() const {
	return tolerance;
}

double DiffusionModule::getKappaN() const {
	return kappaN;
}

double DiffusionModule::getKappaB() const {
	return kappaB;
}

double DiffusionModule::getAlpha() const {
	return alpha;
}

double DiffusionModule::getScale() const {
	return scale;
}

std::string DiffusionModule::getDescription() const {
	std::stringstream s;
	s << "minStep: " << minStep / kpc  << " kpc, ";
	s << "maxStep: " << maxStep / kpc  << " kpc, ";
	s << "tolerance: " << tolerance << "\n";
	
	if (kappaN != 0. or kappaB != 0.) {
	  s << "kappaN: " << kappaN << ", ";
	  s << "kappaB: " << kappaB << "\n";
	  }
	
	if (alpha != 1./3.) {
	  s << "alpha: " << alpha << "\n";
	  }

	if (scale != 1.) {
	  s << "scale: " << scale << "\n";
	  }

	return s.str();
}
