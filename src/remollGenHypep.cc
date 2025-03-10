#include "remollGenHypep.hh"

#include "G4Material.hh"
#include "G4PhysicalConstants.hh"

#include "remollEvent.hh"
#include "remollVertex.hh"
#include "remolltypes.hh"

remollGenHypep::remollGenHypep()
  : remollVEventGen("hyperon_ep") {
  //fThCoM_min =    0.0*deg;
  //fThCoM_max =   180.0*deg;
  //fQ2 = 0;
  //fW2 = 1.7;
  fApplyMultScatt = false;
}

remollGenHypep::~remollGenHypep(){
}

void remollGenHypep::SamplePhysics(remollVertex *vert, remollEvent *evt){

  //G4ParticleDefinition* electron = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  G4ParticleDefinition* proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  G4ParticleDefinition* kaon = G4ParticleTable::GetParticleTable()->FindParticle("kaon+");
  G4ParticleDefinition* hyperon = G4ParticleTable::GetParticleTable()->FindParticle("lambda");
  //G4ParticleDefinition* hyperon = G4ParticleTable::GetParticleTable()->FindParticle("sigma0");
  
  //G4double e_mass = electron->GetPDGMass();
  G4double proton_mass = proton->GetPDGMass();
  G4double kaon_mass = kaon->GetPDGMass();
  G4double hyperon_mass = hyperon->GetPDGMass();
  
  /*
  G4cout << "electron mass is " << electron->GetPDGMass() << G4endl;
  G4cout << "proton mass is " << proton->GetPDGMass() << G4endl;
  G4cout << "kaon mass is " << kaon->GetPDGMass() << G4endl;
  G4cout << "lambda mass is " << lambda->GetPDGMass() << G4endl;
  //GetPDGMass returns a mass in MeV
  */

  G4double beamE =11.0*GeV;

  G4double theta_e, phi_e, theta_K, phi_K;

  theta_e = G4UniformRand()*pi;
  phi_e = 2*pi*(G4UniformRand() - .5);
  theta_K = G4UniformRand()*pi;
  phi_K = 2*pi*(G4UniformRand() - .5);
  
  // Hand Picked Values for Debugging
  /*
  theta_e = 0;
  phi_e = 0;
  theta_K = 3.1415;
  phi_K = 0;
  */

  //G4cout << "beam E = " << beamE << G4endl;
  //G4cout << "theta_e = " << theta_e << G4endl; 
  
  G4double W2_min = std::pow((kaon_mass + hyperon_mass),2); //W2 minimum is the mass required to produce K and Y with 0 momentum
  //G4cout << "W2_min = " << W2_min << G4endl;
  G4double recoilE_max = (std::pow(proton_mass,2) + 2*proton_mass*beamE - W2_min)/(4.0*beamE*std::pow(std::sin(theta_e/2.0),2) + 2.0*proton_mass);
  //G4cout << "recoilE_max = " << recoilE_max << G4endl;
  
  //G4double recoilE = G4UniformRand()*(beamE-recoilE_min) + recoilE_min;
  G4double recoilE = G4UniformRand()*recoilE_max;
  //recoilE = 4.0*GeV; //hand picked for debugging
  //G4cout << beamE << " - " << recoilE << " = " << beamE - recoilE << G4endl; 
  
  G4LorentzVector vecProton_Lab(0.0,0.0,0.0,proton_mass);
  G4LorentzVector vecElectron_Lab(0,0,0,0);
  G4LorentzVector vecPhoton_Lab(0,0,beamE,beamE);
  vecElectron_Lab.setRThetaPhi(recoilE,theta_e,phi_e);
  vecElectron_Lab.setE(recoilE); //~massless electron
  //vecPhoton_Lab(0,0,beamE,beamE);
  vecPhoton_Lab = vecPhoton_Lab - vecElectron_Lab;
  
  //vecPhoton_Lab.setRThetaPhi(beamE - recoilE,-1*theta_e,phi_e);
  //vecPhoton_Lab.setE(beamE - recoilE);

  //G4cout << "vecElectron_Lab.pz() = " << vecElectron_Lab.pz() << G4endl;
  //G4cout << "vecElectron_Lab.px() = " << vecElectron_Lab.px() << G4endl;
  //G4cout << "vecElectron_Lab.py() = " << vecElectron_Lab.py() << G4endl;
  //G4cout << "vecElectron_Lab.e() = " << vecElectron_Lab.e() << G4endl;
  //G4cout << "vecElectron_Lab.m() = " << vecElectron_Lab.m() << G4endl;

  
  //G4cout << "vecPhoton_Lab.pz() = " << vecPhoton_Lab.pz() << G4endl;
  //G4cout << "vecPhoton_Lab.px() = " << vecPhoton_Lab.px() << G4endl;
  //G4cout << "vecPhoton_Lab.py() = " << vecPhoton_Lab.py() << G4endl;
  //G4cout << "vecPhoton_Lab.e() = " << vecPhoton_Lab.e() << G4endl;
  //G4cout << "vecPhoton_Lab.m() = " << vecPhoton_Lab.m() << G4endl;

  G4LorentzVector vecCMFrame = vecProton_Lab+vecPhoton_Lab;
  G4ThreeVector CMtoLabBoost = vecCMFrame.boostVector();
  //G4cout << "vecCMFrame.pz() = " << vecCMFrame.pz() << G4endl;
  //G4cout << "vecCMFrame.px() = " << vecCMFrame.px() << G4endl;
  //G4cout << "vecCMFrame.py() = " << vecCMFrame.py() << G4endl;
  //G4cout << "vecCMframe.e() = " << vecCMFrame.e() << G4endl;

  G4LorentzVector testBoost = vecCMFrame;
  testBoost.boost(-1*CMtoLabBoost);
  //G4cout << "testBoost.pz() = " << testBoost.pz() << G4endl;
  
  G4double Q2 = 4*beamE*recoilE*std::pow(std::sin(theta_e/2.0),2);
  //G4cout << Q2 << " ?= " << -1*vecPhoton_Lab.mag2() << G4endl;
  G4double W2 = vecCMFrame.m2();
  //G4cout << "W2 = " << W2 << G4endl;

  G4LorentzVector vecHyperon_CM, vecKaon_CM;
  G4double CM_momentum = std::sqrt(std::pow( W2 + std::pow(kaon_mass,2) - std::pow(hyperon_mass,2) ,2) - 4*W2*std::pow(kaon_mass,2))/(2*std::sqrt(W2));
  //G4cout << "CM_momentum = " << CM_momentum << G4endl; 
  vecKaon_CM.setRThetaPhi(CM_momentum, theta_K, phi_K);
  vecHyperon_CM.setVect(-1*vecKaon_CM.vect());
  //vecHyperon_CM.setRThetaPhi(CM_momentum, pi - theta_K, phi_K + pi);
  vecKaon_CM.setE(std::sqrt(std::pow(kaon_mass,2) + std::pow(CM_momentum,2)));
  vecHyperon_CM.setE(std::sqrt(std::pow(hyperon_mass,2) + std::pow(CM_momentum,2)));

  //G4cout << "vecHyperon_CM = (" << vecHyperon_CM.px() << "," << vecHyperon_CM.py() << "," << vecHyperon_CM.pz() << "," << vecHyperon_CM.e() << ")" << G4endl;
  //G4cout << "vecKaon_CM = (" << vecKaon_CM.px() << "," << vecKaon_CM.py() << "," << vecKaon_CM.pz() << "," << vecKaon_CM.e() << ")" << G4endl;
  
  //G4cout << "vecHyperon_CM.pz() = " << vecHyperon_CM.pz() << G4endl;
  //G4cout << "vecHyperon_CM.pz() = " << vecHyperon_CM.pz() << G4endl;
  //G4cout << "vecHyperon_CM.px() = " << vecHyperon_CM.px() << G4endl;
  //G4cout << "vecHyperon_CM.py() = " << vecHyperon_CM.py() << G4endl;
  //G4cout << "vecHyperon_CM.e() = " << vecHyperon_CM.e() << G4endl;
  //G4cout << "vecHyperon_CM.m() = " << vecHyperon_CM.m() << G4endl;
  
  G4LorentzVector vecHyperon_Lab, vecKaon_Lab;
  vecHyperon_Lab = vecHyperon_CM.boost(CMtoLabBoost);
  vecKaon_Lab = vecKaon_CM.boost(CMtoLabBoost);

  //G4cout << "vecHyperon_Lab.pz() = " << vecHyperon_Lab.pz() << G4endl;
  //G4cout << "vecHyperon_Lab.px() = " << vecHyperon_Lab.px() << G4endl;
  //G4cout << "vecHyperon_Lab.py() = " << vecHyperon_Lab.py() << G4endl;
  //G4cout << "vecHyperon_Lab.e() = " << vecHyperon_Lab.e() << G4endl;
  //G4cout << "vecHyperon_Lab.m() = " << vecHyperon_Lab.m() << G4endl;

  //G4cout << "vecKaon_Lab.pz() = " << vecKaon_Lab.pz() << G4endl;
  //G4cout << "vecKaon_Lab.px() = " << vecKaon_Lab.px() << G4endl;
  //G4cout << "vecKaon_Lab.py() = " << vecKaon_Lab.py() << G4endl;
  //G4cout << "vecKaon_Lab.e() = " << vecKaon_Lab.e() << G4endl;
  //G4cout << "vecKaon_Lab.m() = " << vecKaon_Lab.m() << G4endl;

  //G4cout << "total Lab px = " << vecHyperon_Lab.px() + vecKaon_Lab.px() + vecElectron_Lab.px() << G4endl;
  //G4cout << "total Lab py = " << vecHyperon_Lab.py() + vecKaon_Lab.py() + vecElectron_Lab.py() << G4endl;
  //G4cout << "total Lab pz = " << vecHyperon_Lab.pz() + vecKaon_Lab.pz() + vecElectron_Lab.pz() << G4endl;
  //G4cout << "total Lab E = " << vecHyperon_Lab.e() + vecKaon_Lab.e() + vecElectron_Lab.e() << G4endl;
  
  
  //G4String Hyperon = "lambda";
  G4ThreeVector Vertex_Prod(0,0,0);
  //G4ThreeVector HyperonMomentum_Lab(0,0,5);
  //G4ThreeVector ElecScatMomentum_Lab(0,0,1);
  //G4ThreeVector KaonMomentum_Lab(0,0,5);

  evt->fBeamE = beamE;
  evt->SetQ2(Q2);
  evt->SetW2(W2);
  evt->SetAsymmetry(0.0);
  evt->SetEffCrossSection(1.0);
  evt->SetRate(3);

  evt->ProduceNewParticle(Vertex_Prod, vecHyperon_Lab.vect(), hyperon->GetParticleName());
  evt->ProduceNewParticle(Vertex_Prod, vecElectron_Lab.vect(), "e-");
  evt->ProduceNewParticle(Vertex_Prod, vecKaon_Lab.vect(), kaon->GetParticleName());
  
  return;
}
