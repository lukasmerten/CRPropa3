#ifndef CRPROPA_JF12FIELD_IMPROVED_H
#define CRPROPA_JF12FIELD_IMPROVED_H

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Grid.h"
#include <crpropa/Units.h>

namespace crpropa {

/**
 @class JF12Field_improved
 @brief JF12Field_improved galactic magnetic field model

 Implements the JF2012 magnetic field model, consisting of a large-scale regular
 and random (striated) field and a small-scale random (turbulent) field.
 See:
 Jansson 2012a, ApJ. 757, A New Model of the Galactic Magnetic Field
 Jansson 2012b, arXiv:1210.7820, The Galactic Magnetic Field

 All three components may individually turned on and off.
 Currently only best fit values of the field paramaters are implemented and
 cannot be changed.

 The field is defined in the usual galactocentric coordinate system with the
 Galactic center at the origin, the x-axis pointing in the opposite direction of
 the Sun, and the z-axis pointing towards Galactic north.
 
 A solenoidal ring-spiral transition and a smooth X-field at z = 0
 according to Kleimann et al. were added.
 
 */
 
class JF12Field_improved: public MagneticField {
private:

	//allow for enabling and disabling of the individual components
	bool useRegular;
	bool useStriated;
	bool useTurbulent;
	
	bool useDisk;
	bool useToroidalHalo;
	bool useX;
	
	//transition widths of the spiral field
	double Delta;
	
	//height paramter of the modified X-field
	double hX;
	
	// disk spiral arms
	double rArms[8]; // radii where each arm crosses the negative x-axis
	
	// angles at which the dividing spirals intersect the r1-ring, 
	//longer array due to cyclicity
	double phi0Arms[11];
	
	//phi integration start angle and corresponding index
	double phi0;
	int idx0;
	
	//phi integration H(phi) = phiCoeff[j] + bDisk[j] * phi
	double phiCoeff[10]; //for phi0Arms[j-1] < phi < phi0Arms[j], hence only 10
	double corr; // correction term enforcing H(phi0) = 0
	
	//inner and outer boundaries disk field
	double r1; 
	double r2;
	
	//transitions region boundaries disk field
	double r1s;
	double r2s;
	
	double pitch;          // pitch angle
	double sinPitch, cosPitch, tanPitch, cotPitch, tan90MinusPitch;

	// Regular field ----------------------------------------------------------
	// disk
	double bDisk[11];       // field strengths of arms at r=5 kpc
	double bRing;          // ring field strength 3<r<5 kpc
	double hDisk, wDisk;   // disk/halo transistion and width
	// toroidal halo
	double bNorth, bSouth; // northern, southern halo field strength
	double rNorth, rSouth; // northern, southern transistion radius
	double wHalo, z0;      // transistion width and vertical scale height
	// poloidal halo
	double bX;             // field strength at origin
	double thetaX0;        // constant elevation angle at r > rXc, z = 0
	double sinThetaX0, cosThetaX0, tanThetaX0, cotThetaX0;
	double rXc;            // radius of varying elevation angle region
	double rX;             // exponential scale height

	// Striated field ---------------------------------------------------------
	double sqrtbeta;       // relative strength of striated field
	ref_ptr<ScalarGrid> striatedGrid;

	// Turbulent field --------------------------------------------------------
	ref_ptr<VectorGrid> turbulentGrid;
	// disk
	double bDiskTurb[8]; // field strengths in arms at r=5 kpc
	double bDiskTurb5;   // field strength at r<5kpc
	double zDiskTurb;	 // Gaussian scale height of disk
	// halo
	double bHaloTurb; // halo field strength
	double rHaloTurb; // exponential scale length
	double zHaloTurb; // Gaussian scale height

public:

	JF12Field_improved(bool usedisk=true, bool usetoroidalhalo=true, bool usex=true, double delta=(3*kpc), double hx=(500*pc));

	// Create and set a random realization for the striated field
	void randomStriated(int seed = 0);

#ifdef CRPROPA_HAVE_FFTW3F
	// Create a random realization for the turbulent field
	void randomTurbulent(int seed = 0);
#endif

	/**
	 * Set a striated grid and activate the striated field component
	 * @param grid	scalar grid containing random +1/-1 values, 100 parsec grid spacing
	 */
	void setStriatedGrid(ref_ptr<ScalarGrid> grid);

	/**
	 * Set a turbulent grid and activate the turbulent field component
	 * @param grid	vector grid containing a random field of Brms = 1
	 */
	void setTurbulentGrid(ref_ptr<VectorGrid> grid);

	ref_ptr<ScalarGrid> getStriatedGrid();
	ref_ptr<VectorGrid> getTurbulentGrid();

	void setUseRegular(bool use);
	void setUseStriated(bool use);
	void setUseTurbulent(bool use);

	bool isUsingRegular();
	bool isUsingStriated();
	bool isUsingTurbulent();

	// Regular field component
	Vector3d getRegularField(const Vector3d& pos) const;

	// Regular and striated field component
	Vector3d getStriatedField(const Vector3d& pos) const;

	// Brms of the turbulent field
	double getTurbulentStrength(const Vector3d& pos) const;

	// Turbulent field component
	Vector3d getTurbulentField(const Vector3d& pos) const;

	// All set field components
	Vector3d getField(const Vector3d& pos) const;
	
	Vector3d getXField(const double r, const double z, const double sinPhi, const double cosPhi) const;
	Vector3d getDiskField(const double r, const double z, const double phi, const double sinPhi, const double cosPhi) const;
	Vector3d getToroidalHaloField(const double r, const double z, const double sinPhi, const double cosPhi) const;
	
	//spiral field functions
	double logisticFunction(double x, double x0, double w) const;
	
	//transition polynomial p_delta(r)
	double p(const double r) const;
	
	//transition polynomial derivative
	double q(const double r) const;
	
	//evaluate H-Integral
	double PhiIntegralH(const double r, const double phi) const;
	
	//return correct field strength b_j in spiral arm
	double getSpiralStrength(const double r, const double phi) const;
	
};

} // namespace crpropa

#endif // CRPROPA_JF12Field_improved_H
