#include "crpropa/module/PhotonOutput1D.h"
#include "crpropa/Units.h"

#include <iostream>
#include <sstream>
#include <cstdio>
#include <stdexcept>

#include "kiss/string.h"
#include "kiss/logger.h"

#ifdef CRPROPA_HAVE_ZLIB
#include <ozstream.hpp>
#endif

using namespace std;

namespace crpropa {

PhotonOutput1D::PhotonOutput1D() : out(&std::cout) {
	KISS_LOG_WARNING << "PhotonOutput1D is deprecated and will be removed in the future. Replace with TextOutput or HDF5Output with features ObserverNucleusVeto + ObserverDetectAll";
}

PhotonOutput1D::PhotonOutput1D(std::ostream &out) : out(&out) {
	KISS_LOG_WARNING << "PhotonOutput1D is deprecated and will be removed in the future. Replace with TextOutput or HDF5Output with features ObserverNucleusVeto + ObserverDetectAll";
}

PhotonOutput1D::PhotonOutput1D(const std::string &filename) : outfile(
	filename.c_str(), std::ios::binary), out(&outfile), filename(filename) {
	KISS_LOG_WARNING << "PhotonOutput1D is deprecated and will be removed in the future. Replace with TextOutput or HDF5Output with features ObserverNucleusVeto + ObserverDetectAll";
	if (kiss::ends_with(filename, ".gz"))
		gzip();

	*out << "#ID\tE\tD\tpID\tpE\tiID\tiE\tiD\n";
	*out << "#\n";
	*out << "# ID          Id of particle (photon, electron, positron)\n";
	*out << "# E           Energy [EeV]\n";
	*out << "# D           Comoving distance to origin [Mpc]\n";
	*out << "# pID         Id of parent particle\n";
	*out << "# pE          Energy [EeV] of parent particle\n";
	*out << "# iID         Id of source particle\n";
	*out << "# iE          Energy [EeV] of source particle\n";
	*out << "# iD          Comoving distance [Mpc] to source\n";
	*out << "#\n";
}

void PhotonOutput1D::process(Candidate *candidate) const {
	int pid = candidate->current.getId();
	if ((pid != 22) and (abs(pid) != 11))
		return;

	char buffer[1024];
	size_t p = 0;

	p += std::sprintf(buffer + p, "%4i\t", pid);
	p += std::sprintf(buffer + p, "%g\t", candidate->current.getEnergy() / EeV);
	p += std::sprintf(buffer + p, "%8.4f\t", candidate->current.getPosition().getR() / Mpc);

	p += std::sprintf(buffer + p, "%10i\t", candidate->created.getId());
	p += std::sprintf(buffer + p, "%8.4f\t", candidate->created.getEnergy() / EeV);

	p += std::sprintf(buffer + p, "%10i\t", candidate->source.getId());
	p += std::sprintf(buffer + p, "%8.4f\t", candidate->source.getEnergy() / EeV);
	p += std::sprintf(buffer + p, "%8.4f\n", candidate->source.getPosition().getR() / Mpc);

#pragma omp critical(FileOutput)
	{
		out->write(buffer, p);
	}

	candidate->setActive(false);
}

void PhotonOutput1D::close() {
	#ifdef CRPROPA_HAVE_ZLIB
		zstream::ogzstream *zs = dynamic_cast<zstream::ogzstream *>(out);
		if (zs) {
			zs->close();
			delete out;
			out = 0;
		}
	#endif
	outfile.flush();
}

string PhotonOutput1D::getDescription() const {
	std::stringstream s;
	s << "PhotonOutput1D: Output file = " << filename;
	return s.str();
}

PhotonOutput1D::~PhotonOutput1D() {
	close();
}

void PhotonOutput1D::gzip() {
	#ifdef CRPROPA_HAVE_ZLIB
		out = new zstream::ogzstream(*out);
	#else
		throw std::runtime_error("CRPropa was build without Zlib compression!");
	#endif
}

} // namespace crpropa
