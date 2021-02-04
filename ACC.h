#ifndef analysis_user_h
#define analysis_user_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"

namespace MA5
{
class ACC : public AnalyzerBase
{
  INIT_ANALYSIS(ACC,"MadAnalysis5job")

 public : 
  virtual MAbool Initialize(const MA5::Configuration& cfg,
                          const std::map<std::string,std::string>& parameters);
  virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
  virtual MAbool Execute(SampleFormat& sample, const EventFormat& event);

 private : 
  // Declaring particle containers
  std::ofstream ACCfile;
  std::vector<const MCParticleFormat*> leptons;
  MAbool is_lepton(const MCParticleFormat* part) const {
     if ( part==0 ) return false;
     if ( !PHYSICS->Id->IsFinalState(part) ) return false;
     if ( (std::abs(part->pdgid())!=13) && (std::abs(part->pdgid())!=11) ) return false;
     if (part->pt() < 25.) return false;
     if (part->abseta() > 2.5) return false;     
     return true; }

  std::vector<const MCParticleFormat*> electrons;
  bool is_electron(const MCParticleFormat* part) const {
     if ( part==0 ) return false;
     if ( !PHYSICS->Id->IsFinalState(part) ) return false;
     if ( (part->pdgid()!=-11)&&(part->pdgid()!=11) ) return false;
     if (part->pt() < 25.) return false;
     if (part->abseta() > 2.5) return false;     
     return true; }
  std::vector<const MCParticleFormat*> muons;
  bool is_muon(const MCParticleFormat* part) const {
     if ( part==0 ) return false;
     if ( !PHYSICS->Id->IsFinalState(part) ) return false;
     if ( (part->pdgid()!=-13)&&(part->pdgid()!=13) ) return false;
     if (part->pt() < 25.) return false;
     if (part->abseta() > 2.5) return false;     
     return true; }

};
}

#endif
