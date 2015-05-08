#include "crpropa/module/DiffusionModule.h"

#include <crpropa/Random.h>

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace crpropa;
	
DiffusionModule::DiffusionModule(ref_ptr<MagneticField> field) :
    field(field)
{
}

void DiffusionModule::process(Candidate *candidate) const {
	double rig = candidate->current.getEnergy() / candidate->current.getCharge();
	double DifCoeff = 6.1e24*pow((rig/4.0e9), 1./3.);
	double stepSize = 2. / c_light * DifCoeff; 


    double step = stepSize;
    if (Random::instance().rand() < .5) {
        step *= -1;
    }
    
    // runge kutta step
    
    Vector3d xi = candidate->current.getPosition();
    if (field->getField(xi).getR() == 0.){
      Vector3d step3d = candidate->current.getDirection() * stepSize;
      	candidate->previous = candidate->current;
	candidate->current.setPosition(xi + step3d);
	candidate->setCurrentStep(stepSize);
	candidate->setNextStep(stepSize);
	candidate->current.setDirection(step3d);

    }else{
    Vector3d k1 = step * field->getField(xi).getUnitVector();
    Vector3d k2 = step * field->getField(xi + k1/2.0).getUnitVector();
    Vector3d k3 = step * field->getField(xi + k2/2.0).getUnitVector();
    Vector3d k4 = step * field->getField(xi + k3).getUnitVector();
    
    // ToDo normalize
    Vector3d step3d = (k1 + 2.0*(k2+k3) + k4)/6.; 
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


	candidate->previous = candidate->current;
	candidate->current.setPosition(xi + step3d);
	candidate->setCurrentStep(stepSize);
	candidate->setNextStep(stepSize);
	candidate->current.setDirection(step3d);
    }
}
