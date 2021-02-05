#ifndef analysis_AFP_h
#define analysis_AFP_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"

namespace MA5
{
class AFP : public AnalyzerBase
{
  INIT_ANALYSIS(AFP,"AFP")

 public:
  virtual bool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
  virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
  virtual bool Execute(SampleFormat& sample, const EventFormat& event);

 private:
  // Declaring particle containers
  std::ofstream myfile;
  std::vector<const MCParticleFormat*> leptons;
  std::vector<const MCParticleFormat*> lepton2;
  std::vector<const MCParticleFormat*> lepton1;

  MAbool is_lepton(const MCParticleFormat* part) const {
     if ( part==0 ) return false;
     if ( !PHYSICS->Id->IsFinalState(part) ) return false;
     if ( (std::abs(part->pdgid())!=13) && (std::abs(part->pdgid())!=11) ) return false;
     if (part->pt() < 20.) return false;
     if (part->abseta() > 2.5) return false;
     return true; }

  std::vector<const MCParticleFormat*> jets;

  MAbool is_jet(const MCParticleFormat* part) const {
     if ( part==0 ) return false;
     if ( !PHYSICS->Id->IsFinalState(part) ) return false;
     if ( (part->pdgid()!=-5)&&(part->pdgid()!=-4)&&(part->pdgid()!=-3)&&(part->pdgid()!=-2)&&(part->pdgid()!=-1)&&\
          (part->pdgid()!=5)&&(part->pdgid()!=4)&&(part->pdgid()!=3)&&(part->pdgid()!=2)&&(part->pdgid()!=1)&&\
          (part->pdgid()!=21) ) return false;
     if (part->pt() < 15.) return false;
     return true; }

  std::vector<const MCParticleFormat*> initialquarks;

  MAbool is_initialquark(const MCParticleFormat* part) const {
     if ( part==0 ) return false;
     if ( !PHYSICS->Id->IsInitialState(part) ) return false;
     if ( (part->pdgid()!=-1)&&(part->pdgid()!=-2)&&(part->pdgid()!=-3)&&(part->pdgid()!=-4)&&(part->pdgid()!=-5)&&\
          (part->pdgid()!=1)&&(part->pdgid()!=2)&&(part->pdgid()!=3)&&(part->pdgid()!=4)&&(part->pdgid()!=5)\
        ) return false;
     return true; }  

};
}

#endif
