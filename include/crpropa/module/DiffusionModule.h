#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include <crpropa/Module.h>
#include <crpropa/magneticField/MagneticField.h>
#include <crpropa/Units.h>

namespace crpropa {
class DiffusionModule : public Module{

private:
	    std::vector<double> a, b, bs; /*< Cash-Karp coefficients */
	    ref_ptr<MagneticField> field;
	    double minStep;
	    double maxStep;
	    double tolerance;
	    double kappa;
	    std::vector<double> BTensor;
	    std::vector<double> eta;
	    

public:
	    DiffusionModule(ref_ptr<crpropa::MagneticField> field, double tolerance = 1e-4, 
			    double minStep = 50*parsec, double maxStep = 10 * kpc, double kappa = 0.);

	    void process(crpropa::Candidate *candidate) const;
	   
	    //void calBTensor(double rigidity, std::vector<double> BT) const;
	    void tryStep(const Vector3d &Pos, Vector3d &POut, Vector3d &DirErr, Vector3d &TVec,Vector3d &NVec,Vector3d &BVec,double &TStep, double t, double z ) const;

	    void setMinimumStep(double minStep);
	    void setMaximumStep(double maxStep);
	    void setTolerance(double tolerance);
	    void setKappa(double kappa);
	    //void setBTensor(std::vector<double> V) const;

	    double getMinimumStep() const;
	    double getMaximumStep() const;
	    double getTolerance() const;
	    double getKappa() const;
	    //std::vector<double> getBTensor() const;
	    Vector3d TVec, NVec, BVec, PosIn, PosOut, DirOut, DirErr;
	    double TStep, NStep, BStep;

}; 
} //namespace crpropa
