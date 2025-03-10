#ifndef __REMOLLGENHYPEP_HH
#define __REMOLLGENHYPEP_HH

#include "remollVEventGen.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Exp.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"


class remollGenHypep : public remollVEventGen {
public:
  remollGenHypep();
  virtual ~remollGenHypep();

private:
  void SamplePhysics(remollVertex *, remollEvent *);
};

#endif//__REMOLLGENHYPEP_HH
