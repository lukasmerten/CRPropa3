#include "crpropa/module/FieldLineIntegrator.h"

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



FieldLineIntegrator::FieldLineIntegrator(ref_ptr<MagneticField> field, double tolerance, 
				 double minStep, double maxStep, double direction) :
  minStep(0)
{
  setField(field);
  setMaximumStep(maxStep);
  setMinimumStep(minStep);
  setTolerance(tolerance);
  setDirection(direction);
  }


void FieldLineIntegrator::process(Candidate *candidate) const {

    // save the new previous particle state
	ParticleState &current = candidate->current;
	candidate->previous = current;
	
	double step = clip(candidate->getNextStep(), minStep, maxStep);
	Vector3d PosIn = current.getPosition();
	Vector3d PosOut, Err;


	double r=42.; //arbitrary number larger than one
	double newstep = step;

	while (r >1) {
		step = newstep;
		tryStep(PosIn, PosOut, Err, step / c_light);
		
		r = Err.getR() / tolerance;
		newstep = step * 0.95 * pow(r, -0.2);
		newstep = clip(newstep, 0.1 * step, 5 * step);
		newstep = clip(newstep, minStep, maxStep);


   		if (step == minStep)
			break;
	} 
	
	current.setPosition(PosOut);
	candidate->setCurrentStep(step);
	candidate->setNextStep(newstep);


}


void  FieldLineIntegrator::tryStep(const Vector3d &PosIn, Vector3d &PosOut, Vector3d &Err, double h) const {

	Vector3d k[] = {Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.)};
	PosOut = PosIn;
	Err = Vector3d(0.);

	//calculate the sum k_i * b_i
	for (size_t i = 0; i < 6; i++) {

		Vector3d y_n = PosIn;

		for (size_t j = 0; j < i; j++)
		  	y_n += k[j] * a[i * 6 + j] * h;
		
		// update k_i = direction of the regular magnetic mean field
		Vector3d BField(0.);
		try {
		  	BField = field->getField(y_n, 0.); // redshift evolution not taken into account
		} 
		catch (std::exception &e) {
		  	std::cerr << "DiffusionSDE: Exception in getField." << std::endl;
		  	std::cerr << e.what() << std::endl;
		}
		
		k[i] = BField.getUnitVector() * c_light * direction;

		PosOut += k[i] * b[i] * h;
		Err +=  (k[i] * (b[i] - bs[i])) * h / kpc;
		
	}
	//std::cout << "PosErr = " << PosErr <<"\n";
}



void FieldLineIntegrator::setMinimumStep(double min) {
	if (min < 0)
		throw std::runtime_error("DiffusionSDE: minStep < 0 ");
	if (min > maxStep)
		throw std::runtime_error("DiffusionSDE: minStep > maxStep");
	minStep = min;
}

void FieldLineIntegrator::setMaximumStep(double max) {
	if (max < minStep)
		throw std::runtime_error("DiffusionSDE: maxStep < minStep");
	maxStep = max;
}

void FieldLineIntegrator::setTolerance(double tol) {
	if ((tol > 1) or (tol < 0))
		throw std::runtime_error(
				"DiffusionSDE: tolerance error not in range 0-1");
	tolerance = tol;
}

void FieldLineIntegrator::setDirection(bool d) {
	if (d) {
		direction = 1.;
	}
	else {
		direction = -1;
	}
}
	

void FieldLineIntegrator::setField(ref_ptr<MagneticField> f) {
	field = f;
}

double FieldLineIntegrator::getMinimumStep() const {
	return minStep;
}

double FieldLineIntegrator::getMaximumStep() const {
	return maxStep;
}

double FieldLineIntegrator::getTolerance() const {
	return tolerance;
}

double FieldLineIntegrator::getDirection() const {
	return direction;
}



std::string FieldLineIntegrator::getDescription() const {
	std::stringstream s;
	s << "minStep: " << minStep / kpc  << " kpc, ";
	s << "maxStep: " << maxStep / kpc  << " kpc, ";
	s << "tolerance: " << tolerance << "\n";

	return s.str();
}
