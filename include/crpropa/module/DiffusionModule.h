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
public:
        class Y {
	public:
		Vector3d x, u; /*< phase-point: position and direction */

		Y() {
		}

		Y(const Vector3d &x, const Vector3d &u) :
				x(x), u(u) {
		}

		Y(double f) :
				x(Vector3d(f, f, f)), u(Vector3d(f, f, f)) {
		}

		Y operator *(double f) const {
			return Y(x * f, u * f);
		}

		Y &operator +=(const Y &y) {
			x += y.x;
			u += y.u;
			return *this;
		}
	};

private:
	    std::vector<double> a, b, bs; /*< Cash-Karp coefficients */
	    ref_ptr<MagneticField> field;
	    double minStep;
	    double maxStep;
	    double tolerance;
	    double kappa;

public:
	    DiffusionModule(ref_ptr<crpropa::MagneticField> field, double tolerance = 1e-4, 
			    double minStep = 50*parsec, double maxStep = 10 * kpc, double kappa = 0.);

	    void process(crpropa::Candidate *candidate) const;

	    Y dYdt(const Y &y, ParticleState &p, double z) const;

	    void tryStep(const Y &y, Y &out, Y &error, double t,
			ParticleState &p, double z) const;
	   
	    void setMinimumStep(double minStep);
	    void setMaximumStep(double maxStep);
	    void setTolerance(double tolerance);
	    void setKappa(double kappa);

	    double getMinimumStep() const;
	    double getMaximumStep() const;
	    double getTolerance() const;
	    double getKappa() const;


}; 
} //namespace crpropa
