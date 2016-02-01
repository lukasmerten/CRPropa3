#include "crpropa/module/DiffusionModule.h"

#include <crpropa/Random.h>

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <stdexcept>


using namespace crpropa;

const double cash_karp_a[] = { 0., 0., 0., 0., 0., 0., 1. / 5., 0., 0., 0., 0.,
		0., 3. / 40., 9. / 40., 0., 0., 0., 0., 3. / 10., -9. / 10., 6. / 5.,
		0., 0., 0., -11. / 54., 5. / 2., -70. / 27., 35. / 27., 0., 0., 1631.
				/ 55296., 175. / 512., 575. / 13824., 44275. / 110592., 253.
				/ 4096., 0. };

const double cash_karp_b[] = { 37. / 378., 0, 250. / 621., 125. / 594., 0., 512.
		/ 1771. };

const double cash_karp_bs[] = { 2825. / 27648., 0., 18575. / 48384., 13525.
		/ 55296., 277. / 14336., 1. / 4. };


void DiffusionModule::tryStep(const Y &y, Y &out, Y &error, double h,
		ParticleState &particle, double z) const {
	std::vector<Y> k;
	k.reserve(6);

	out = y;
	error = Y(0);

	// calculate the sum of b_i * k_i
	for (size_t i = 0; i < 6; i++) {

		Y y_n = y;
		for (size_t j = 0; j < i; j++)
			y_n += k[j] * a[i * 6 + j] * h;

		// update k_i
		k[i] = dYdt(y_n, particle, z);

		out.x += k[i].x * b[i] * h;
		out.u = (out.x - y.x) / h / c_light;
		std::cout <<"out = " << out.u <<"\n";
		error.x += k[i].x * (b[i] - bs[i]) * h;
		error.u = error.x / h / c_light;
		//error += k[i] * (b[i] - b[i]) * h;
	}
}
	/* // Runge Kutta 4th-order
const double cash_karp_a[] = {0., 0., 0., 0., 0.,
			      1./2., 1./2., 0., 0., 0.,
			      1./2., 0., 1./2., 0., 0.,
			      1., 0., 0., 1., 0.};

const double cash_karp_b[] = {1./6., 1./3., 1./3., 1./6. };

const double cash_karp_bs[] = {1./6., 1./3., 1./3., 1./6. };


void DiffusionModule::tryStep(const Y &y, Y &out, Y &error, double h,
		ParticleState &particle, double z) const {
	std::vector<Y> k;
	k.reserve(4);

	out = y;
	error = Y(0);

	// calculate the sum of b_i * k_i
	for (size_t i = 0; i < 4; i++) {

		Y y_n = y;
		for (size_t j = 0; j < i; j++)
			y_n += k[j] * a[i * 5 + j] * h;

		// update k_i
		k[i] = dYdt(y_n, particle, z);

		out += k[i] * b[i] * h;
		//std::cout <<"out = " << out.x / kpc <<"\n";
		error += k[i] * (b[i] - bs[i]) * h;
	}
	}*/

DiffusionModule::Y DiffusionModule::dYdt(const Y &y, ParticleState &p, double z) const {
	// normalize direction vector to prevent numerical losses
	Vector3d B(0, 0, 0);
	try {
		B = field->getField(y.x, z);
	} catch (std::exception &e) {
		std::cerr << "PropagationCK: Exception in getField." << std::endl;
		std::cerr << e.what() << std::endl;
	}
	Vector3d velocity = B.getUnitVector() * c_light;
	// change direction according to magnetic field line
	Vector3d dudt = B.getUnitVector();
	return Y(velocity, dudt);
}


DiffusionModule::DiffusionModule(ref_ptr<MagneticField> field, double tolerance, 
				 double minStep, double maxStep) :
  field(field)
{
  setMaximumStep(maxStep);
  setMinimumStep(minStep);
  setTolerance(tolerance);
  
// load CK-coefficients
  a.assign(cash_karp_a, cash_karp_a + 36);
  b.assign(cash_karp_b, cash_karp_b + 6);
  bs.assign(cash_karp_bs, cash_karp_bs + 6);
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
		candidate->setCurrentStep(fabs(step));
		candidate->setNextStep(maxStep);
		return;
	}
	
	Y yIn(current.getPosition(), current.getDirection());
	Y yOut, yErr;
	double h = step / c_light;
	//std::cout <<"h = " << h / 31557600.<<"\n";
	double hTry, r;
	double z = candidate->getRedshift();
	
	double rig = current.getEnergy() / current.getCharge();
	double DifCoeff = 6.1e24 * pow((rig / 4.0e9), 1./3.);
	double B_xx = pow(2 * DifCoeff, 0.5);
	double propStep = B_xx * Random::instance().randNorm();
	do {
		hTry = h;
		std::cout <<"hTry = " << hTry / 31557600. <<"\n";
		std::cout <<"propStep = " << propStep*pow(hTry, 0.5) / kpc <<"\n";
		tryStep(yIn, yOut, yErr, (propStep*pow(hTry, 0.5))/c_light, current, z);
		// determine absolute direction error relative to tolerance
		//r = yErr.u.getR() / tolerance;
		r = yErr.u.getR() / tolerance;
		std::cout <<"yErr.x " << yErr.x.getR() / kpc <<"\n";
		std::cout <<"r = " << r <<"\n";
		// new step size to keep the error close to the tolerance
		h *= 0.95 * pow(r, -0.2);
		// limit change of new step size
		h = clip(h, 0.1 * hTry, 5 * hTry);

	} while (r > 1 && h > minStep / c_light);
		//}while (1 == 2);
	/*
// rectilinear propagation if magnetic field is B=0
	if(step3d.getR2() != step3d.getR2()){
	  Vector3d dir = current.getDirection();
	  //current.setPosition(xi + dir * stepSize);
	  current.setPosition(xi + dir * fabs(step));
	  current.setDirection(dir);
	  candidate->setCurrentStep(fabs(step));
	  candidate->setNextStep(maxStep);
	  return;
	}
	*/
	current.setPosition(yOut.x);
	current.setDirection(yOut.u.getUnitVector());
	candidate->setCurrentStep(hTry * c_light);
	candidate->setNextStep(h * c_light);
}


void DiffusionModule::setMinimumStep(double min) {
	if (min < 0)
		throw std::runtime_error("PropagationCK: minStep < 0 ");
	if (min > maxStep)
		throw std::runtime_error("PropagationCK: minStep > maxStep");
	minStep = min;
}

void DiffusionModule::setMaximumStep(double max) {
	if (max < minStep)
		throw std::runtime_error("PropagationCK: maxStep < minStep");
	maxStep = max;
}

void DiffusionModule::setTolerance(double tol) {
	if ((tol > 1) or (tol < 0))
		throw std::runtime_error(
				"PropagationCK: target error not in range 0-1");
	tolerance = tol;
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

