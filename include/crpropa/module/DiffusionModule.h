#pragma once
#include <iostream>
#include <vector>
#include <cmath>

#include <crpropa/Module.h>
#include <crpropa/magneticField/MagneticField.h>

class DiffusionModule : public crpropa::Module{
	private:
	    crpropa::ref_ptr<crpropa::MagneticField> field;
	public:
		DiffusionModule(crpropa::ref_ptr<crpropa::MagneticField> field);

        void process(crpropa::Candidate *candidate) const;
}; 
