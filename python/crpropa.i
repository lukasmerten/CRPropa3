/* CRPropa3 SWIG interface (for Python) */

/* Content:
 *
 * 1. SWIG settings and workarounds
 * 2. SWIG and CRPropa headers
 * 3. Pretty print for Python
 * 4. Magnetic Lens and Particle Maps Container
 *
 */


/* 1. SWIG settings and workarounds */
%include "1_swig.i"

/* 2: SWIG and CRPropa headers */

%include stl.i
%include std_set.i
%include std_multiset.i
%include std_map.i
%include std_pair.i
%include std_multimap.i
%include std_vector.i
%include std_string.i
%include std_list.i
%include stdint.i
%include std_container.i
%include "exception.i"

#ifdef CRPROPA_HAVE_QUIMBY
%import (module="quimby") quimby.i
#endif

%{
#include "crpropa/module/NuclearDecay.h"
#include "crpropa/module/ElectronPairProduction.h"
#include "crpropa/module/PhotoDisintegration.h"
#include "crpropa/module/PhotoPionProduction.h"
#include "crpropa/module/Redshift.h"
#include "crpropa/module/BreakCondition.h"
#include "crpropa/module/Boundary.h"
#include "crpropa/module/Observer.h"
#include "crpropa/module/OutputTXT.h"
#include "crpropa/module/OutputShell.h"
#include "crpropa/module/OutputROOT.h"
#include "crpropa/module/OutputCRPropa2.h"
#include "crpropa/module/PhotonDINT.h"
#include "crpropa/module/PhotonDINT1D.h"
#include "crpropa/module/PhotonEleCa.h"
#include "crpropa/module/PhotonOutput1D.h"
#include "crpropa/module/SimplePropagation.h"
#include "crpropa/module/PropagationCK.h"
#include "crpropa/module/Tools.h"
#include "crpropa/module/DiffusionModule.h"

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/magneticField/MagneticFieldGrid.h"
#include "crpropa/magneticField/QuimbyMagneticField.h"
#include "crpropa/magneticField/AMRMagneticField.h"
#include "crpropa/magneticField/JF12Field.h"
#include "crpropa/magneticField/TurbulentMagneticField.h"

#include "crpropa/Referenced.h"
#include "crpropa/Candidate.h"
#include "crpropa/ParticleState.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Module.h"
#include "crpropa/ModuleList.h"
#include "crpropa/Random.h"
#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Source.h"
#include "crpropa/Common.h"
#include "crpropa/Cosmology.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/PhotonPropagation.h"
#include "crpropa/Grid.h"
#include "crpropa/GridTools.h"

#include "crpropa/Version.h"
%}

%{
#include <iostream>
#include <iomanip>
%}

%ignore operator<<;
%ignore operator>>;
%ignore *::operator=;
%ignore operator crpropa::Source*;
%ignore operator crpropa::SourceList*;
%ignore operator crpropa::SourceInterface*;
%ignore operator crpropa::SourceFeature*;
%ignore operator crpropa::Candidate*;
%ignore operator crpropa::Module*;
%ignore operator crpropa::ModuleList*;
%ignore operator crpropa::MagneticField*;

%feature("ref")   crpropa::Referenced "$this->addReference();"
%feature("unref") crpropa::Referenced "$this->removeReference();"


%include "crpropa/Vector3.h"


%include "crpropa/Referenced.h"
%include "crpropa/Units.h"
%include "crpropa/Common.h"
%include "crpropa/Cosmology.h"
%include "crpropa/PhotonBackground.h"
%include "crpropa/PhotonPropagation.h"
%include "crpropa/Random.h"
%include "crpropa/ParticleState.h"
%include "crpropa/ParticleID.h"
%include "crpropa/ParticleMass.h"

%template(CandidateVector) std::vector< crpropa::ref_ptr<crpropa::Candidate> >;
%template(CandidateRefPtr) crpropa::ref_ptr<crpropa::Candidate>;
%include "crpropa/Candidate.h"

%template(ModuleRefPtr) crpropa::ref_ptr<crpropa::Module>;
%template(stdModuleList) std::list< crpropa::ref_ptr<crpropa::Module> >;
%feature("director") crpropa::Module;
%include "crpropa/Module.h"

%implicitconv crpropa::ref_ptr<crpropa::MagneticField>;
%template(MagneticFieldRefPtr) crpropa::ref_ptr<crpropa::MagneticField>;
%include "crpropa/magneticField/MagneticField.h"

%include "crpropa/Grid.h"
%include "crpropa/GridTools.h"

%implicitconv crpropa::ref_ptr<crpropa::Grid<crpropa::Vector3<float> > >;
%template(VectorGridRefPtr) crpropa::ref_ptr<crpropa::Grid<crpropa::Vector3<float> > >;
%template(VectorGrid) crpropa::Grid<crpropa::Vector3<float> >;

%implicitconv crpropa::ref_ptr<crpropa::Grid<float> >;
%template(ScalarGridRefPtr) crpropa::ref_ptr<crpropa::Grid<float> >;
%template(ScalarGrid) crpropa::Grid<float>;

%include "crpropa/magneticField/MagneticFieldGrid.h"
%include "crpropa/magneticField/QuimbyMagneticField.h"
%include "crpropa/magneticField/AMRMagneticField.h"
%include "crpropa/magneticField/JF12Field.h"
%include "crpropa/magneticField/TurbulentMagneticField.h"

%include "crpropa/module/BreakCondition.h"
%include "crpropa/module/Boundary.h"
%include "crpropa/module/Observer.h"
%include "crpropa/module/SimplePropagation.h"
%include "crpropa/module/PropagationCK.h"
%include "crpropa/module/OutputTXT.h"
%include "crpropa/module/OutputShell.h"
%include "crpropa/module/OutputROOT.h"
%include "crpropa/module/OutputCRPropa2.h"
%include "crpropa/module/PhotonDINT.h"
%include "crpropa/module/PhotonDINT1D.h"
%include "crpropa/module/PhotonEleCa.h"
%include "crpropa/module/PhotonOutput1D.h"
%include "crpropa/module/ElectronPairProduction.h"
%include "crpropa/module/NuclearDecay.h"
%include "crpropa/module/PhotoPionProduction.h"
%include "crpropa/module/PhotoDisintegration.h"
%include "crpropa/module/Redshift.h"
%include "crpropa/module/Tools.h"
%include "crpropa/module/DiffusionModule.h"

%template(SourceInterfaceRefPtr) crpropa::ref_ptr<crpropa::SourceInterface>;
%feature("director") crpropa::SourceInterface;
%template(SourceFeatureRefPtr) crpropa::ref_ptr<crpropa::SourceFeature>;
%feature("director") crpropa::SourceFeature;
%include "crpropa/Source.h"

%template(ModuleListRefPtr) crpropa::ref_ptr<crpropa::ModuleList>;
%include "crpropa/ModuleList.h"

%include "crpropa/Version.h"
//%pythoncode %{
//    __version__ = g_GIT_DESC 
//%}
%include "2_headers.i"

/* 3. Pretty print for Python */

%define __REPR__( classname ) 

%pythoncode %{
globals()["classname"["classname".find('::')+2:]].__repr__ = globals()["classname"["classname".find('::')+2:]].getDescription
%}

%enddef

%define VECTOR3__REPR__( classname ) 

%template(Vector3d) classname<double>;
%template(Vector3f) classname<float>;

%pythoncode %{

def Vector3__repr__(self):
    return "Vector(%.3g, %.3g, %.3g)" % (self.x, self.y, self.z)
Vector3d.__repr__ = Vector3__repr__
Vector3f.__repr__ = Vector3__repr__

%}

%enddef


%include "3_repr.i"

/* 4. Magnetic Lens */
%include "4_lens.i"

#ifdef WITH_GALACTIC_LENSES

%ignore Pixelization::nPix();

%pythoncode %{

def Pixelization_nonStaticnPix(self, order=None):
  if order == None:
    return Pixelization_nPix(self.getOrder())
  else:
    return Pixelization_nPix(order)
Pixelization.nPix = Pixelization_nonStaticnPix

MagneticLens.transformModelVector = MagneticLens.transformModelVector_numpyArray

ParticleMapsContainer.getMap = ParticleMapsContainer.getMap_numpyArray
ParticleMapsContainer.getParticleIds = ParticleMapsContainer.getParticleIds_numpyArray
ParticleMapsContainer.getEnergies = ParticleMapsContainer.getEnergies_numpyArray
ParticleMapsContainer.getRandomParticles = ParticleMapsContainer.getRandomParticles_numpyArray

%}

#endif // WITH_GALACTIC_LENSES_

