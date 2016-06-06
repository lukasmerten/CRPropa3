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
				 double minStep, double maxStep, double epsilon) :
  minStep(0)
{
  setTurbulent(false);
  setField(field);
  setMaximumStep(maxStep);
  setMinimumStep(minStep);
  setTolerance(tolerance);
  setEpsilon(epsilon);
  setScale(1.);
  setAlpha(1./3.);

  }

DiffusionModule::DiffusionModule(ref_ptr<MagneticField> field, ref_ptr<MagneticField> turbField,  double tolerance, double minStep, double maxStep, double epsilon) :
  minStep(0)
{
  setField(field);
  setTurbulentField(turbField);
  setMaximumStep(maxStep);
  setMinimumStep(minStep);
  setTolerance(tolerance);
  setEpsilon(epsilon);
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
	Vector3d DirIn = current.getDirection();
	double z = candidate->getRedshift();
	double rig = current.getEnergy() / current.getCharge();

	double BTensor[] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

	calculateBTensor(rig, BTensor, PosIn, DirIn, z);
	
	//std::cout << "BTensor diagonal = " << BTensor[0] <<"\t" <<BTensor[4] <<"\t" <<BTensor[8] <<"\n";
	
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
	
	// Exception: Rectilinear propagation in case of veanishing magnetic field.
	if (TVec.getR() != TVec.getR()) {
	  Vector3d dir = current.getDirection();
	  current.setPosition(current.getPosition() + dir * step);
	  candidate->setCurrentStep(step);
	  candidate->setNextStep(step);
	  return;
	}
	// Integration of the SDE with a Mayorama-Euler-method
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
	for (size_t i = 0; i < 6; i++) {
		Vector3d y_n = PosIn;
		for (size_t j = 0; j < i; j++)
		  y_n += k[j] * a[i * 6 + j] * propStep;
		Pos[i] = y_n;
		
		// update k_i = direction of the regular magnetic mean field
		Vector3d BField(0.);
		try {
		  BField = field->getField(y_n, z);
		} 
		catch (std::exception &e) {
		  std::cerr << "PropagationCK: Exception in getField." << std::endl;
		  std::cerr << e.what() << std::endl;
		}
		
		k[i] = BField.getUnitVector();
		POut += k[i] * b[i] * propStep;
		PosErr +=  (k[i] * (b[i] - bs[i])) * propStep;
			
	}
	TVec = POut.getUnitVector();
	
	// Choose a random perpendicular vector as the Normal-vector

	NVec = TVec.cross( Random::instance().randVector() );
	NVec = NVec.getUnitVector();

	// Calculate the Binormal-vector
	BVec = (TVec.cross(NVec)).getUnitVector();
}


void DiffusionModule::calculateBTensor(double r, double BTen[], Vector3d pos, Vector3d dir, double z) const {
  if (isTurbulent == false) {
    double DifCoeff = scale * 6.1e24 * pow((std::abs(r) / 4.0e9), alpha);
    BTen[0] = pow( 2  * DifCoeff, 0.5);
    BTen[4] = pow(2 * epsilon * DifCoeff, 0.5);
    BTen[8] = pow(2 * epsilon * DifCoeff, 0.5);
    return;
  }
  // Diffusion model from Snodin et al. "Global diffusion of cosmic rays in random magnetic fields" (2016)
  else {
    double turbStrength = (turbField->getField(pos, z)).getR();
    Vector3d regField = field->getField(pos, z);
    double regStrength = regField.getR();
    double sinTheta = std::sin(dir.getAngleTo(regField));
    //std::cout << "b0 = " << turbStrength << "\t B0 = " << regStrength << "\n";
    double mu = turbStrength*turbStrength / ( turbStrength*turbStrength +  regStrength*regStrength );
    //std::cout << "mu = " << mu <<"\n";
    double RL = std::abs(r) / regStrength / c_light * sinTheta; //Larmor-radius
    //std::cout << "Larmorradius = " << RL / parsec <<"\n";
    double L = 272*parsec; //from lmax of the JF12 turbulent field component
    double kappa_0 = (0.0017 + 0.75 * (RL/L) ) * c_light*L;
    //std::cout << "kappa_0 = " << kappa_0 <<"\n";
    double kappa_parallel = kappa_0 + 1./3. * pow((RL/L), 1./3.) * (1-mu)/mu * c_light*L;
    double chi = 2.35; // obtained from FLRW simulations. Analytically this is chi = 4*D_0/l_c .
    double kappa_perp = ( kappa_0 * mu + 0.19 * (1-mu) * pow((RL / L), 0.61) * c_light*L ) / ( 1 + chi * (1-mu)/ mu);
    BTen[0] = pow( 2  * kappa_parallel, 0.5);
    BTen[4] = pow( 2 * kappa_perp, 0.5);
    BTen[8] = pow( 2 * kappa_perp, 0.5);

  }
  
  
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

void DiffusionModule::setEpsilon(double e) {
	if ((e > 1) or (e < 0))
		throw std::runtime_error(
				"DiffusionModule: epsilon not in range 0-1");
	epsilon = e;
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

void DiffusionModule::setTurbulent(bool b) {
  isTurbulent = b;
}

void DiffusionModule::setField(ref_ptr<MagneticField> f) {
	field = f;
}

void DiffusionModule::setTurbulentField(ref_ptr<crpropa::MagneticField> field){
  try {
    turbField = field;
    setTurbulent(true);
  }
  catch (std::exception &e) {
    std::cerr << "DiffusionModule: Exception in setTurbulentField" << std::endl;
    std::cerr << e.what() << std::endl;
  }
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

double DiffusionModule::getEpsilon() const {
	return epsilon;
}


double DiffusionModule::getAlpha() const {
	return alpha;
}

double DiffusionModule::getScale() const {
	return scale;
}

bool DiffusionModule::getTurbulent() const {
  return isTurbulent;
}



std::string DiffusionModule::getDescription() const {
	std::stringstream s;
	s << "minStep: " << minStep / kpc  << " kpc, ";
	s << "maxStep: " << maxStep / kpc  << " kpc, ";
	s << "tolerance: " << tolerance << "\n";
	
	if (epsilon != 0.1) {
	  s << "epsilon: " << epsilon << ", ";
	  }
	
	if (alpha != 1./3.) {
	  s << "alpha: " << alpha << "\n";
	  }

	if (scale != 1.) {
	  s << "scale: " << scale << "\n";
	  }

	return s.str();
}
