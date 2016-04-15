#include "crpropa/module/DiffusionModule.h"

#include <crpropa/Random.h>

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <stdexcept>


using namespace crpropa;

const double cash_karp_a[] = { 0., 0., 0., 0., 0., 0., 1. / 5., 0., 0., 0., 0.,
		0., 3. / 40., 9. / 40., 0., 0., 0., 0., 3. / 10., -9. / 10., 6. / 5.,
		0., 0., 0., -11. / 54., 5. / 2., -70. / 27., 35. / 27., 0., 0., 1631.
				/ 55296., 175. / 512., 575. / 13824., 44275. / 110592., 253.
				/ 4096., 0. };

const double cash_karp_b[] = { 37. / 378., 0, 250. / 621., 125. / 594., 0., 512.
		/ 1771. };

const double cash_karp_bs[] = { 2825. / 27648., 0., 18575. / 48384., 13525.
		/ 55296., 277. / 14336., 1. / 4. };


DiffusionModule::DiffusionModule(ref_ptr<MagneticField> field, double tolerance, 
				 double minStep, double maxStep, double kappa) :
  field(field)
{
  setMaximumStep(maxStep);
  setMinimumStep(minStep);
  setTolerance(tolerance);
  setKappa(kappa);
  //setBTensor(std::vector<double> (9));
  
// load CK-coefficients
  a.assign(cash_karp_a, cash_karp_a + 36);
  b.assign(cash_karp_b, cash_karp_b + 6);
  bs.assign(cash_karp_bs, cash_karp_bs + 6);
}



void DiffusionModule::process(Candidate *candidate) const {
	// save the new previous particle state
	ParticleState &current = candidate->current;
	candidate->previous = current;
	
	double step = clip(candidate->getNextStep(), minStep, maxStep);

	// rectilinear propagation for neutral particles
	if (current.getCharge() == 0) {
		Vector3d dir = current.getDirection();
		current.setPosition(current.getPosition() + dir * step);
		candidate->setCurrentStep(step);
		candidate->setNextStep(maxStep);
		return;
	}

	Vector3d PosIn = current.getPosition();
	Vector3d PosOut = Vector3d(0.);
	Vector3d DirOut = Vector3d(0.);
	Vector3d DirErr = Vector3d(0.);

	double h = step / c_light;
	double hTry, r;

	double z = candidate->getRedshift();
	double rig = current.getEnergy() / current.getCharge();

	//calBTensor(rig, BTensor);
	double DifCoeff = 6.1e24 * pow((rig / 4.0e9), 1./3.);
	std::vector<double> BTensor (9);
	BTensor.at(0) = pow( 2  * DifCoeff, 0.5);
	BTensor.at(4) = pow(2 * kappa * DifCoeff, 0.5);
	BTensor.at(8) = pow(2 * kappa * DifCoeff, 0.5);
	
	//~ std::cout << "Btensor= " << BTensor[0] <<"\n";
	std::vector<double> eta (3);
	for(size_t i=0; i < eta.size(); i++) {
	  eta.at(i) =  Random::instance().randNorm();
	  //~ std::cout <<"eta" << eta.at(i) <<"\n";
	}

	double TStep = BTensor[0] * eta[0];
	double NStep = BTensor[4] * eta[1];
	double BStep = BTensor[8] * eta[2];
	//~ std::cout <<"TStep " << TStep << "\n";
	//std::cout <<"eta" << eta.at(0) <<"\n";
    //~ std::cout << "Btensor= " << BTensor[0] <<"\n";

	Vector3d TVec(0.);
	Vector3d NVec(0.);
	Vector3d BVec(0.);
	 

	do {
	  hTry = h;
	  tryStep(PosIn, PosOut, DirErr, TVec, NVec, BVec, TStep, z, hTry);
	  r = DirErr.getR() / hTry / tolerance;
	  h *= 0.95 * pow(r, -0.2);
	  h = clip(h, 0.1 * hTry, 5 * hTry);
	} while (r > 1 && h > minStep /c_light);
	//std::cout <<"TStep = " << TStep <<"\n";
	//std::cout <<"hTry = " << hTry <<"\n";
	//std::cout <<"PosIn = " << PosIn <<"\n";
	//std::cout <<"TVec = " << TVec <<"\n";
	//std::cout <<"NVec = " << NVec <<"\n";
	//std::cout <<"BVec = " << BVec <<"\n";
	Vector3d PO = PosIn + (TVec * std::abs(TStep) + NVec * NStep + BVec * BStep) * pow(hTry, 0.5);
	
	//std::cout <<"TVec = " << TVec <<"\n";
	//std::cout <<"NVec = " << NVec <<"\n";
	//std::cout <<"BVec = " << BVec <<"\n";
	//std::cout <<"PosOut = " << PO <<"\n";
	DirOut = (PO -PosIn).getUnitVector();
	current.setPosition(PO);
	current.setDirection(DirOut);
	candidate->setCurrentStep(hTry * c_light);
	candidate->setNextStep(h * c_light);
}




/*
void DiffusionModule::calBTensor(double r, std::vector<double> BT) const {
  
  //double DifCoeff = 6.1e24 * pow((r / 4.0e9), 1./3.);
  // std::vector<double> BTen = std::vector<double> (9) ;
//std::vector<double> BTen = getBTensor();
  
  BTen.at(0) = pow( 2  * DifCoeff, 0.5);
  BTen.at(4) = pow(2 * kappa * DifCoeff, 0.5);
  BTen.at(8) = pow(2 * kappa * DifCoeff, 0.5);
  
  //setBTensor(BTen);
  setBTensor(std::vector<double> (9));
}
*/


void DiffusionModule::tryStep(const Vector3d &PosIn, Vector3d &POut, Vector3d &DirErr, Vector3d &TVec,Vector3d &NVec,Vector3d &BVec,double &TStep, double z, double h) const {
  	std::vector<Vector3d> k;
	k.reserve(6);
	std::vector<Vector3d> Pos;
	Pos.reserve(6);


	//PosOut = Vector3d(0., 0., 0.);
	//DirErr = Vector3d(0.);
	//~ std::cout <<"TStep " << TStep << "\n";
    //~ std::cout << "Btensor= " << BTensor[0] <<"\n";

	for (size_t i = 0; i < 6; i++) {
		Vector3d y_n = PosIn;
		for (size_t j = 0; j < i; j++)
		  y_n += k[j] * a[i * 6 + j] * TStep * pow(h, 0.5);
		Pos[i] = y_n;
		
		// update k_i
		Vector3d BField(0.);
		try {
		  BField = field->getField(y_n, z);
		} catch (std::exception &e) {
		  std::cerr << "PropagationCK: Exception in getField." << std::endl;
		  std::cerr << e.what() << std::endl;
		}
		
		k[i] = BField.getUnitVector();
		//~ std::cout <<"k[" << i << "]=\t" << k[i] << '\n';
		//~ std::cout <<"TStep" << TStep << "\n";
		//~ std::cout <<"sqrt(h)" << pow(h, 0.5)  << "\n";
		POut += k[i] * b[i] * TStep * pow(h, 0.5);
		//std::cout <<"Pos[i]= " << Pos[i] << '\n';
		DirErr +=  (k[i] * (b[i] - bs[i]) * TStep * pow(h, 0.5));
		//Pos[i+1] = POut;
		//std::cout << "Pos[i-1], Pos[i]" << Pos[i-1] << "\t" <<Pos[i]  <<"\n";
		if (i > 0)
		  //std::cout << ( Pos[i-1] - Pos[i] ).getR() <<"\n";
		  NVec += ( k[i] - k[i-1] ) / ( Pos[i-1] - Pos[i] ).getR();


	}
	TVec = POut.getUnitVector();
	//std::cout << NVec.getR()<<"\n";
	
	if (NVec.getR() == 0.){
	  NVec = TVec.cross( Random::instance().randVector() );
	    }
	
	NVec = NVec.getUnitVector();
	BVec = (TVec.cross(NVec)).getUnitVector();
	
}


void DiffusionModule::setMinimumStep(double min) {
	if (min < 0)
		throw std::runtime_error("DiffusionModule: minStep < 0 ");
	if (min > maxStep)
		throw std::runtime_error("DiffusionModule: minStep > maxStep");
	minStep = min;
}

void DiffusionModule::setMaximumStep(double max) {
	if (max < minStep)
		throw std::runtime_error("DiffusionModule: maxStep < minStep");
	maxStep = max;
}

void DiffusionModule::setTolerance(double tol) {
	if ((tol > 1) or (tol < 0))
		throw std::runtime_error(
				"DiffusionModule: tolerance error not in range 0-1");
	tolerance = tol;
}

void DiffusionModule::setKappa(double kap) {
	if ((kap > 1) or (kap < 0))
		throw std::runtime_error(
				"DiffusionModule: ratio error not in range 0-1");
	kappa = kap;
}
/*
void DiffusionModule::setBTensor(std::vector<double> V) const {
  if (V.size() != 9)
    throw std::runtime_error(
			     "DiffusionModule: BTensor needs nine elements!");
  BTensor = V;
}
*/

double DiffusionModule::getMinimumStep() const {
	return minStep;
}

double DiffusionModule::getMaximumStep() const {
	return maxStep;
}

double DiffusionModule::getTolerance() const {
	return tolerance;
}

double DiffusionModule::getKappa() const {
	return kappa;
}
/*
std::vector<double> DiffusionModule::getBTensor() const {
	return BTensor;
}
*/
