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
	    ref_ptr<MagneticField> field;
	    ref_ptr<MagneticField> turbField;
	    MagneticField NullField;
	    double minStep;
	    double maxStep;
	    double tolerance;
	    double epsilon;
	    double alpha;
	    double scale;
	    bool isTurbulent;
	    

public:
	    DiffusionModule(ref_ptr<crpropa::MagneticField> field, double tolerance = 1e-4, 
	    		    double minStep=(10*pc), double maxStep=(1*kpc), double epsilon=0.1);

	    DiffusionModule(ref_ptr<crpropa::MagneticField> field, ref_ptr<crpropa::MagneticField> turbField, double tolerance = 1e-4, 
			    double minStep = 10*parsec, double maxStep = 1*kpc, double epsilon = 0.1);

	    void process(crpropa::Candidate *candidate) const;
	   
	    void tryStep(const Vector3d &Pos, Vector3d &POut, Vector3d &PosErr, Vector3d &PosTest, Vector3d &TVec,Vector3d &NVec,Vector3d &BVec, double z, double propStep ) const;
	    
	    void calculateBTensor(double rig, double BTen[], Vector3d pos, Vector3d dir, double z) const;

	    void setMinimumStep(double minStep);
	    void setMaximumStep(double maxStep);
	    void setTolerance(double tolerance);
	    void setEpsilon(double kappa);
	    void setAlpha(double alpha);
	    void setScale(double Scale);
	    void setTurbulent(bool isTurb);
	    void setField(ref_ptr<crpropa::MagneticField> field);
	    void setTurbulentField(ref_ptr<crpropa::MagneticField> field);

	    double getMinimumStep() const;
	    double getMaximumStep() const;
	    double getTolerance() const;
	    double getEpsilon() const;
	    double getAlpha() const;
	    double getScale() const;
	    bool getTurbulent() const;
	    std::string getDescription() const;

}; 

} //namespace crpropa
