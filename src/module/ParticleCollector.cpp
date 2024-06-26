#include "crpropa/module/ParticleCollector.h"
#include "crpropa/module/TextOutput.h"
#include "crpropa/Units.h"

namespace crpropa {

ParticleCollector::ParticleCollector() : nBuffer(10e6), clone(false), recursive(false)  {
        container.reserve(nBuffer); // for 1e6 candidates ~ 500MB of RAM
}

ParticleCollector::ParticleCollector(const std::size_t nBuffer) : clone(false), recursive(false)  {
	container.reserve(nBuffer);
}

ParticleCollector::ParticleCollector(const std::size_t nBuffer, const bool clone) : recursive(false) {
	container.reserve(nBuffer);
}

ParticleCollector::ParticleCollector(const std::size_t nBuffer, const bool clone, const bool recursive) {
	container.reserve(nBuffer);
}

void ParticleCollector::process(Candidate *c) const {
#pragma omp critical(ModifyContainer)
        {
		if(clone)
		       	container.push_back(c->clone(recursive));
		else
			container.push_back(c);
        }
}

void ParticleCollector::process(ref_ptr<Candidate> c) const {
	ParticleCollector::process((Candidate*) c);
}

void ParticleCollector::reprocess(Module *action) const {
	for (ParticleCollector::iterator itr = container.begin(); itr != container.end(); ++itr){
		if (clone)
			action->process((*(itr->get())).clone(false));
		else
        	        action->process(itr->get());
	}
}

void ParticleCollector::dump(const std::string &filename) const {
	TextOutput output(filename.c_str(), Output::Everything);
	reprocess(&output);
	output.close();
}

void ParticleCollector::load(const std::string &filename){
	TextOutput::load(filename.c_str(), this);
}

ParticleCollector::~ParticleCollector() {
        clearContainer();
}

std::size_t ParticleCollector::size() const {
        return container.size();
}

ref_ptr<Candidate> ParticleCollector::operator[](const std::size_t i) const {
	return container[i];
}

void ParticleCollector::clearContainer() {
        container.clear();
}

std::vector<ref_ptr<Candidate> >& ParticleCollector::getContainer() const {
        return container;
}

void ParticleCollector::setClone(bool b) {
        clone = b;
}

bool ParticleCollector::getClone() const {
        return clone;
}

std::string ParticleCollector::getDescription() const {
        return "ParticleCollector";
}

ParticleCollector::iterator ParticleCollector::begin() {
	return container.begin();
}

ParticleCollector::const_iterator ParticleCollector::begin() const {
	return container.begin();
}

ParticleCollector::iterator ParticleCollector::end() {
	return container.end();
}

ParticleCollector::const_iterator ParticleCollector::end() const {
	return container.end();
}

void ParticleCollector::getTrajectory(ModuleList* mlist, std::size_t i, Module *output) const {
	ref_ptr<Candidate> c_tmp = container[i]->clone();

	c_tmp->restart();

	mlist->add(output);
	mlist->run(c_tmp);
	mlist->remove(mlist->size()-1);
}

void ParticleCollector::getTrajectory(ref_ptr<ModuleList> mlist, std::size_t i, ref_ptr<Module> output) const {
	ParticleCollector::getTrajectory((ModuleList*) mlist, i, (Module*) output);
}

} // namespace crpropa
