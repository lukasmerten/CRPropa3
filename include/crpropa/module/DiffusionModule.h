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
 @class MaximumTrajectoryLength
 @brief Deactivates the candidate beyond a maximum trajectory length

 This modules deactivates the candidate at a given maximum trajectory length.
 In that case the property ("Deactivated", module::description) is set.
 It also limits the candidates next step size to ensure the maximum trajectory length is no exceeded.
 */


class DiffusionModule : public Module{

private:
    //std::vector<double> a, b, bs; /*< Cash-Karp coefficients */
	    ref_ptr<MagneticField> field;
	    double minStep;
	    double maxStep;
	    double tolerance;
	    double kappaN;
	    double kappaB;
	    double alpha;
	    double scale;
	   
	    

public:
	    DiffusionModule(ref_ptr<crpropa::MagneticField> field, double tolerance = 1e-4, 
			    double minStep = 50*parsec, double maxStep = 10 * kpc, double kappa = 0.);

	    void process(crpropa::Candidate *candidate) const;
	   
	  
	    void tryStep(const Vector3d &Pos, Vector3d &POut, Vector3d &PosErr, Vector3d &TVec,Vector3d &NVec,Vector3d &BVec, double z, double propStep ) const;

	    void setMinimumStep(double minStep);
	    void setMaximumStep(double maxStep);
	    void setTolerance(double tolerance);
	    void setKappaN(double kappa);
	    void setKappaB(double kappa);
	    void setAlpha(double alpha);
	    void setScale(double Scale);

	    double getMinimumStep() const;
	    double getMaximumStep() const;
	    double getTolerance() const;
	    double getKappaN() const;
	    double getKappaB() const;
	    double getAlpha() const;
	    double getScale() const;
	    std::string getDescription() const;

}; 

} //namespace crpropa
