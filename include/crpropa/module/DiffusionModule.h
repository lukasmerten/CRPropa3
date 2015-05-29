#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include <crpropa/Module.h>
#include <crpropa/magneticField/MagneticField.h>
#include <crpropa/Units.h>

class DiffusionModule : public crpropa::Module{
	private:
	    crpropa::ref_ptr<crpropa::MagneticField> field;
	    std::string mode;
	    std::string singleString;
	    std::string logString;
	    std::string check;
	public:
	    DiffusionModule(crpropa::ref_ptr<crpropa::MagneticField> field, std::string mode);

	    void process(crpropa::Candidate *candidate) const;
	    void InputCheck(std::string check);
}; 
