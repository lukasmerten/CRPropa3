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
	    ref_ptr<MagneticField> field;
	    std::string mode;
	    std::string singleString;
	    std::string logString;
	    std::string check;

	public:
	    DiffusionModule(ref_ptr<crpropa::MagneticField> field, std::string mode, 
			    double minStep = 0.1*parsec, double maxStep = 10 * kpc);

	    void process(crpropa::Candidate *candidate) const;
	    void InputCheck(std::string check);

	private:
	    double minStep;
	    double maxStep;

	public:
	void setMinimumStep(double minStep);
	void setMaximumStep(double maxStep);

	double getMinimumStep() const;
	double getMaximumStep() const;


}; 
} //namespace crpropa
