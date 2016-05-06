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
 @class DiffusionModule
 @brief Propagates pseudo-particles as tracers of the phase space density.

 */


class DiffusionModule : public Module{

private:
    //std::vector<double> a, b, bs; /*< Cash-Karp coefficients */
	    ref_ptr<MagneticField> field;
	    double minStep;
	    double maxStep;
	    double tolerance;
	    double epsilon;
	    double alpha;
	    double scale;	   
	    

public:
	    DiffusionModule(ref_ptr<crpropa::MagneticField> field, double tolerance = 1e-4, 
			    double minStep = 10*parsec, double maxStep = 1 * kpc, double epsilon = 0.1);

	    void process(crpropa::Candidate *candidate) const;
	   
	  
	    void tryStep(const Vector3d &Pos, Vector3d &POut, Vector3d &PosErr, Vector3d &TVec,Vector3d &NVec,Vector3d &BVec, double z, double propStep ) const;

	    void setMinimumStep(double minStep);
	    void setMaximumStep(double maxStep);
	    void setTolerance(double tolerance);
	    void setEpsilon(double kappa);
	    void setAlpha(double alpha);
	    void setScale(double Scale);

	    double getMinimumStep() const;
	    double getMaximumStep() const;
	    double getTolerance() const;
	    double getEpsilon() const;
	    double getAlpha() const;
	    double getScale() const;
	    std::string getDescription() const;

}; 

} //namespace crpropa
