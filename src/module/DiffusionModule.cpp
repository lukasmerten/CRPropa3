#include "crpropa/module/DiffusionModule.h"

#include <crpropa/Random.h>

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <stdexcept>


using namespace crpropa;
	
DiffusionModule::DiffusionModule(ref_ptr<MagneticField> field, std::string m, 
				 double minStep, double maxStep) :
  field(field), mode(m)
{
setMaximumStep(maxStep);
setMinimumStep(minStep);
  try {
    InputCheck(mode);
  }
  catch (std::invalid_argument& e) {
    std::cerr << "Error in DiffusionModule:" << std::endl;
    std::terminate();
  }
}

void DiffusionModule::process(Candidate *candidate) const {
	// save the new previous particle state
	ParticleState &current = candidate->current;
	candidate->previous = current;
	
	Vector3d xi = current.getPosition();
	double E = current.getEnergy();
	double C = current.getCharge();
	double rig = E / C;
	double DifCoeff = 6.1e24 * pow((rig / 4.0e9), 1./3.);
	double stepSize;
	
	double step = clip(candidate->getNextStep(), minStep, maxStep);

// rectilinear propagation for neutral particles
	if (C == 0) {
		Vector3d dir = current.getDirection();
		current.setPosition(xi + dir * step);
		candidate->setCurrentStep(fabs(step));
		candidate->setNextStep(maxStep);
		return;
	}

// Diffusion for charged particles 
	if(mode == singleString){ // single step distribution
	  stepSize = 2. / c_light * DifCoeff;
	  //step = stepSize;
	}else if(mode == logString){ // log-step distribution
	  //stepSize = 1. / c_light * DifCoeff;
	  //step = -1. * stepSize * log(1-Random::instance().rand());
	  stepSize =  -1. * 1. / c_light * DifCoeff * log(1-Random::instance().rand());
	}

   	if (Random::instance().rand() < .5) { // normal Diffusion
	//if (Random::instance().rand() < .0) { // along magnetic field lines: North direction
	//if (Random::instance().rand() < 1.) { //against magnetic field lines: South direction
	  step *= -1;
	}

// runge kutta step
    	Vector3d k1 = step * field->getField(xi).getUnitVector();
	Vector3d k2 = step * field->getField(xi + k1/2.0).getUnitVector();
	Vector3d k3 = step * field->getField(xi + k2/2.0).getUnitVector();
	Vector3d k4 = step * field->getField(xi + k3).getUnitVector();
    
	Vector3d step3d = (k1 + 2.0*(k2+k3) + k4)/6.; 

// rectilinear propagation if magnetic field is B=0
	if(step3d.getR2() != step3d.getR2()){
	  Vector3d dir = current.getDirection();
	  //current.setPosition(xi + dir * stepSize);
	  current.setPosition(xi + dir * fabs(step));
	  current.setDirection(dir);
	  candidate->setCurrentStep(fabs(step));
	  //candidate->setNextStep(maxStep);
	  candidate->setNextStep(fabs(step));
	  return;
	}

	current.setPosition(xi + step3d);
	//candidate->setCurrentStep(stepSize);
	candidate->setCurrentStep(fabs(step));
	candidate->setNextStep(stepSize);
	current.setDirection(step3d);
}

void DiffusionModule::InputCheck(std::string c) {
      check = c;
      singleString = "single";
      logString = "log";
      
      if (c == "single"){
	return;
      }
      if (c == "log"){
	return;
      }
      else{
	throw std::invalid_argument("Invalid second keyword: Use \"single\" or \"log\" instead");
      }
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

double DiffusionModule::getMinimumStep() const {
	return minStep;
}

double  DiffusionModule::getMaximumStep() const {
	return maxStep;
}


// Alternative Algorithm should be tested for efficiency
 /*
    // Runge Kutta 3/8
    Vector3d xi = candidate->current.getPosition();
    Vector3d k1 = step * field->getField(xi).getUnitVector();
    Vector3d k2 = step * field->getField(xi + k1/3.0).getUnitVector();
    Vector3d k3 = step * field->getField(xi - k1/3.0 + k2).getUnitVector();
    Vector3d k4 = step * field->getField(xi + k1 -k2 + k3).getUnitVector();
    
    // ToDo normalize
    Vector3d step3d = (k1 + 3.0*(k2+k3) + k4)/8.; 
    
    // Runge Kutta Fehlberg
    Vector3d xi = candidate->current.getPosition();
    Vector3d k1 = step * field->getField(xi).getUnitVector();
    Vector3d k2 = step * field->getField(xi + k1/4.0).getUnitVector();
    Vector3d k3 = step * field->getField(xi + k1*3./32. + k2*9./32.).getUnitVector();
    Vector3d k4 = step * field->getField(xi + k1*1932./2197. - k2*7200./2197. + k3*7296./2197.).getUnitVector();
    Vector3d k5 = step * field->getField(xi + k1*439./216. - k2*8. + k3*3680./513. - k4*845./4104.).getUnitVector();
    Vector3d k6 = step * field->getField(xi - k1*8./27. +k2*2. - k3*3544./2565. + k4*1859./4104. - k5*11./40.).getUnitVector();

    
    // 5-th order
    //Vector3d step3d = (k1*16./135. + k3*6656./12825. + k4*28561./56430. - k5*9./50. + k6*2./55.);
    // 4-th order
    Vector3d step3d = (k1*25./216. + k3*1408./2565. + k4*2197./4104. - k5*1./5.);
    */
