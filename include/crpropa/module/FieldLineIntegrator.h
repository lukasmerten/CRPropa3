#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include <crpropa/Module.h>
#include <crpropa/magneticField/MagneticField.h>
#include <crpropa/Units.h>

namespace crpropa {

/**
 @class FieldLineIntegrator
 @brief Calculates the magnetic field line for a given start position for a given direction
 */


class FieldLineIntegrator : public Module{

private:
	    ref_ptr<MagneticField> field;
	    double minStep; // minStep/c_light is the minimum integration timestep
	    double maxStep; // maxStep/c_light is the maximum integration timestep
	    double tolerance; // tolerance is criterion for step adjustment. Step adjustment takes place when the tangential vector of the magnetic field line is calculated.
	    double direction;
	    
	    

public:
	    FieldLineIntegrator(ref_ptr<crpropa::MagneticField> field, double tolerance = 1e-4, 
	    		    double minStep=(10*pc), double maxStep=(1*kpc), double direction = true);

	    void process(crpropa::Candidate *candidate) const;
	   
	    void tryStep(const Vector3d &Pos, Vector3d &POut, Vector3d &Err, double h ) const;
	    
	    void setMinimumStep(double minStep);
	    void setMaximumStep(double maxStep);
	    void setTolerance(double tolerance);
	    void setDirection(bool direction);
	    void setField(ref_ptr<crpropa::MagneticField> field);

	    double getMinimumStep() const;
	    double getMaximumStep() const;
	    double getTolerance() const;
	    double getDirection() const;
	    std::string getDescription() const;

}; 

} //namespace crpropa
