#include "crpropa/module/Tools.h"
#include "crpropa/Clock.h"

#include <iostream>
#include <sstream>

using namespace std;

namespace crpropa {

PerformanceModule::~PerformanceModule() {
	double total = 0;
	for (size_t i = 0; i < modules.size(); i++) {
		_module_info &m = modules[i];
		total += m.time;
	}
	cout << "Performance for " << calls << " calls:" << endl;
	for (size_t i = 0; i < modules.size(); i++) {
		_module_info &m = modules[i];
		cout << " - " << floor((1000 * m.time / total) + 0.5) / 10 << "% -> "
				<< m.module->getDescription() << ": " << (m.time / calls)
				<< endl;
	}
}

void PerformanceModule::add(Module *module) {
	_module_info info;
	info.module = module;
	info.time = 0;
	modules.push_back(info);
}

void PerformanceModule::process(Candidate *candidate) const {
	vector<double> times(modules.size());
	for (size_t i = 0; i < modules.size(); i++) {
		_module_info &m = modules[i];
		double start = Clock::getInstance().getMillisecond();
		m.module->process(candidate);
		double end = Clock::getInstance().getMillisecond();
		times[i] = end - start;
	}

#pragma omp critical(PerformanceModule)
	{
		for (size_t i = 0; i < modules.size(); i++) {
			_module_info &m = modules[i];
			m.time += times[i];
		}
		calls++;
	}
}

string PerformanceModule::getDescription() const {
	stringstream sstr;
	sstr << "PerformanceModule (";
	for (size_t i = 0; i < modules.size(); i++) {
		_module_info &m = modules[i];
		if (i > 0)
			sstr << ", ";
		sstr << m.module->getDescription();
	}
	sstr << ")";
	return sstr.str();
}

// ----------------------------------------------------------------------------
ParticleFilter::ParticleFilter() {

}
ParticleFilter::ParticleFilter(const std::set<int> &ids) : ids(ids) {

}
void ParticleFilter::addId(int id) {
	ids.insert(id);
}
void ParticleFilter::removeId(int id) {
	ids.erase(id);
}

std::set<int> &ParticleFilter::getIds() {
	return ids;
}

void ParticleFilter::process(Candidate* candidate) const {
	if (ids.find(candidate->current.getId()) == ids.end())
		reject(candidate);
	else
		accept(candidate);
}

string ParticleFilter::getDescription() const {
	stringstream sstr;
	sstr << "ParticleFilter: ";
	for (std::set<int>::const_iterator i = ids.begin(); i != ids.end(); i++) {
		sstr << *i << ", ";
	}
	sstr << ")";
	return sstr.str();
}

// ----------------------------------------------------------------------------
EmissionMapFiller::EmissionMapFiller(EmissionMap *emissionMap) : emissionMap(emissionMap) {

}

void EmissionMapFiller::setEmissionMap(EmissionMap *emissionMap) {
	this->emissionMap = emissionMap;
}

void EmissionMapFiller::process(Candidate* candidate) const {
	if (emissionMap) {
		#pragma omp critical(EmissionMap)
		{
			emissionMap->fillMap(candidate->source);
		}
	}
}

string EmissionMapFiller::getDescription() const {
	return "EmissionMapFiller";
}

} // namespace crpropa
