#include "crpropa/module/OutputTXT.h"
#include "crpropa/Units.h"

#include <stdio.h>

namespace crpropa {

  TrajectoryOutput::TrajectoryOutput(std::string name): SpatialScale(Mpc), EnergyScale(EeV) {
	setDescription("Trajectory output");
	fout.open(name.c_str());
	fout << "# D\tID\tE\tX\tY\tZ\tPx\tPy\tPz\n"
	     << "#\n"
	     << "# D           Trajectory length\n"
	     << "# ID          Particle type (PDG MC numbering scheme)\n"
	     << "# E           Energy [EeV]\n"
	     << "# X, Y, Z     Position [Mpc]\n"
	     << "# Px, Py, Pz  Heading (unit vector of momentum)\n"
	     << "#\n";
}

  TrajectoryOutput::TrajectoryOutput(std::string name, double s, double E): SpatialScale(s), EnergyScale(E) {
	setDescription("Trajectory output");
	fout.open(name.c_str());
	fout << "# D\tID\tE\tX\tY\tZ\tPx\tPy\tPz\n"
	     << "#\n"
	     << "# D           Trajectory length\n"
	     << "# ID          Particle type (PDG MC numbering scheme)\n"
	     << "# E           Energy ["<<EnergyScale/TeV << "*TeV]\n"
	     << "# X, Y, Z     Position ["<< SpatialScale/kpc <<"*kpc]\n"
	     << "# Px, Py, Pz  Heading (unit vector of momentum)\n"
	     << "#\n";
}
 
TrajectoryOutput::~TrajectoryOutput() {
	fout.close();
}

void TrajectoryOutput::process(Candidate *c) const {
	char buffer[1024];
	size_t p = 0;

	p += sprintf(buffer + p, "%8.3f\t", c->getTrajectoryLength() / SpatialScale);
	p += sprintf(buffer + p, "%10i\t", c->current.getId());
	p += sprintf(buffer + p, "%8.4f\t", c->current.getEnergy() / EnergyScale);
	Vector3d pos = c->current.getPosition() / SpatialScale;
	p += sprintf(buffer + p, "%8.4f\t%8.4f\t%8.4f\t", pos.x, pos.y, pos.z);
	const Vector3d &dir = c->current.getDirection();
	p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\n", dir.x, dir.y, dir.z);

#pragma omp critical
	fout.write(buffer, p);
}

void TrajectoryOutput::endRun() {
	fout.flush();
}

ConditionalOutput::ConditionalOutput(std::string fname, std::string cond) :
  condition(cond), SpatialScale(Mpc), EnergyScale(EeV) {
	setDescription(
			"Conditional output, condition: " + cond + ", filename: " + fname);
	fout.open(fname.c_str());
	fout << "# D\tID\tID0\tE\tE0\tX\tY\tZ\tX0\tY0\tZ0\tPx\tPy\tPz\tP0x\tP0y\tP0z\tz\n"
	     << "#\n"
	     << "# D           Trajectory length [Mpc]\n"
	     << "# ID          Particle type (PDG MC numbering scheme)\n"
	     << "# E           Energy [EeV]\n"
	     << "# X, Y, Z     Position [Mpc]\n"
	     << "# Px, Py, Pz  Heading (unit vector of momentum)\n"
	     << "# z           Current redshift\n"
	     << "# Initial state: ID0, E0, ...\n"
	     << "#\n";
}

  ConditionalOutput::ConditionalOutput(std::string fname, double s, double E, std::string cond) :
    condition(cond), SpatialScale(s), EnergyScale(E) {
	setDescription(
			"Conditional output, condition: " + cond + ", filename: " + fname);
	fout.open(fname.c_str());
	fout << "# D\tID\tID0\tE\tE0\tX\tY\tZ\tX0\tY0\tZ0\tPx\tPy\tPz\tP0x\tP0y\tP0z\tz\n"
	     << "#\n"
	     << "# D           Trajectory length ["<< SpatialScale/kpc <<"*kpc]\n"
	     << "# ID          Particle type (PDG MC numbering scheme)\n"
	     << "# E           Energy ["<<EnergyScale/TeV << "*TeV]\n"
	     << "# X, Y, Z     Position ["<< SpatialScale/kpc <<"*kpc]\n"
	     << "# Px, Py, Pz  Heading (unit vector of momentum)\n"
	     << "# z           Current redshift\n"
	     << "# Initial state: ID0, E0, ...\n"
	     << "#\n";
}
  
ConditionalOutput::~ConditionalOutput() {
	fout.close();
}

void ConditionalOutput::process(Candidate *c) const {
	if (not (c->hasProperty(condition)))
		return;

	c->removeProperty(condition);

	char buffer[1024];
	size_t p = 0;

	p += sprintf(buffer + p, "%8.3f\t", c->getTrajectoryLength() / SpatialScale);
	p += sprintf(buffer + p, "%10i\t", c->current.getId());
	p += sprintf(buffer + p, "%10i\t", c->source.getId());
	p += sprintf(buffer + p, "%8.4f\t", c->current.getEnergy() / EnergyScale);
	p += sprintf(buffer + p, "%8.4f\t", c->source.getEnergy() / EnergyScale);
	Vector3d pos = c->current.getPosition() / SpatialScale;
	p += sprintf(buffer + p, "%9.4f\t%9.4f\t%9.4f\t", pos.x, pos.y, pos.z);
	Vector3d ipos = c->source.getPosition() / SpatialScale;
	p += sprintf(buffer + p, "%9.4f\t%9.4f\t%9.4f\t", ipos.x, ipos.y, ipos.z);
	Vector3d dir = c->current.getDirection();
	p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\t", dir.x, dir.y, dir.z);
	Vector3d idir = c->source.getDirection();
	p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\t", idir.x, idir.y, idir.z);
	p += sprintf(buffer + p, "%1.3f\n", c->getRedshift());

#pragma omp critical
	fout.write(buffer, p);
}

void ConditionalOutput::endRun() {
	fout.flush();
}

  TrajectoryOutput1D::TrajectoryOutput1D(std::string filename) : SpatialScale(Mpc), EnergyScale(EeV){
	setDescription("TrajectoryOutput, filename: " + filename);
	fout.open(filename.c_str());
	fout << "#X\tID\tE\n"
	     << "#\n"
	     << "# X  Position [Mpc]\n"
	     << "# ID Particle type\n"
	     << "# E  Energy [EeV]\n";
}

  TrajectoryOutput1D::TrajectoryOutput1D(std::string filename, double s, double E) : SpatialScale(s), EnergyScale(E) {
	setDescription("TrajectoryOutput, filename: " + filename);
	fout.open(filename.c_str());
	fout << "#X\tID\tE\n"
	     << "#\n"
	     << "# X  Position [" << SpatialScale/kpc <<"*kpc]\n"
	     << "# ID Particle type\n"
	     << "# E  Energy ["<<EnergyScale/TeV << "*TeV]\n";
}
  
TrajectoryOutput1D::~TrajectoryOutput1D() {
	fout.close();
}

void TrajectoryOutput1D::process(Candidate *c) const {
	char buffer[1024];
	size_t p = 0;
	p += sprintf(buffer + p, "%8.4f\t", c->current.getPosition().x / SpatialScale);
	p += sprintf(buffer + p, "%10i\t", c->current.getId());
	p += sprintf(buffer + p, "%8.4f\n", c->current.getEnergy() / EnergyScale);

#pragma omp critical
	fout.write(buffer, p);
}

void TrajectoryOutput1D::endRun() {
	fout.flush();
}

  EventOutput1D::EventOutput1D(std::string filename) : SpatialScale(Mpc), EnergyScale(EeV) {
	setDescription("Conditional output, filename: " + filename);
	fout.open(filename.c_str());
	fout << "#ID\tE\tD\tID0\tE0\n"
	     << "#\n"
	     << "# ID  Particle type\n"
	     << "# E   Energy [EeV]\n"
	     << "# D   Comoving source distance [Mpc]\n"
	     << "# ID0 Initial particle type\n"
	     << "# E0  Initial energy [EeV]\n";
}

  EventOutput1D::EventOutput1D(std::string filename, double s, double E) : SpatialScale(s), EnergyScale(E) {
	setDescription("Conditional output, filename: " + filename);
	fout.open(filename.c_str());
	fout << "#ID\tE\tD\tID0\tE0\n"
	     << "#\n"
	     << "# ID  Particle type\n"
	     << "# E   Energy [EeV]\n"
	     << "# D   Comoving source distance ["<< SpatialScale/kpc <<"*kpc]\n"
	     << "# ID0 Initial particle type\n"
	     << "# E0  Initial energy ["<<EnergyScale/TeV << "*TeV]\n";
}

EventOutput1D::~EventOutput1D() {
	fout.close();
}

void EventOutput1D::process(Candidate *c) const {
	char buffer[1024];
	size_t p = 0;

	p += sprintf(buffer + p, "%10i\t", c->current.getId());
	p += sprintf(buffer + p, "%8.4f\t", c->current.getEnergy() / EnergyScale);
	p += sprintf(buffer + p, "%9.4f\t", c->source.getPosition().x / SpatialScale);
	p += sprintf(buffer + p, "%10i\t", c->source.getId());
	p += sprintf(buffer + p, "%8.4f\n", c->source.getEnergy() / EnergyScale);

#pragma omp critical
	fout.write(buffer, p);
}

void EventOutput1D::endRun() {
	fout.flush();
}

} // namespace crpropa
