#include "SampleAnalyzer/User/Analyzer/AFP.h"
#include <cmath>
#include <iostream>
#include "SampleAnalyzer/Commons/Service/PDGService.h"

using namespace MA5;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
MAbool AFP::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  // Initializing PhysicsService for MC
  PHYSICS->mcConfig().Reset();  
	
  // definition of the multiparticle "invisible"
  PHYSICS->mcConfig().AddInvisibleId(12);

  // ===== Signal region ===== //
  Manager()->AddRegionSelection("SR1");

  // ===== Selections ===== //
  Manager()->AddCut("$N_l=2$");
  Manager()->AddCut("OSSF");

  // ===== Histograms ===== //
  Manager()->AddHisto("pT",   20,500.,1500.);

//  Manager()->AddHisto("AFP_Left",  20,500.,1500.);
//  Manager()->AddHisto("AFP_Right", 20,500.,1500.);

  return true;
}



// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
MAbool AFP::Execute(SampleFormat& sample, const EventFormat& event)
{
  MAfloat32 EvWeight = 1.0;
  if (weighted_events_ && event.mc()!=0) EvWeight *= event.mc()->weight();

  //EvWeight *= sample.mc()->xsection(); //* 1000. * lumi;
  if (sample.mc()!=0) sample.mc()->addWeightedEvents(EvWeight);
  Manager()->InitializeForNewEvent(EvWeight);

  // Clearing particle containers
  leptons.clear();
//  lepton1.clear();
//  lepton2.clear();

  // Filling particle containers
  for (MAuint32 i=0;i<event.mc()->particles().size();i++)
  {
    if (is_lepton(&(event.mc()->particles()[i])))
      leptons.push_back(&(event.mc()->particles()[i]));
  }

  // Sorting particles
  // Sorting particle collection according to PTordering
  SORTER->sort(leptons,PTordering);
  SORTER->sort(jets,PTordering);
  // for getting 2th particle
//  lepton2=SORTER->rankFilter(leptons,2,PTordering);
  // for getting 1th particle
//  lepton1=SORTER->rankFilter(leptons,1,PTordering);


  // DiLepton reconstruction
  MAuint32 hasOS   = 0;
  MAuint32 hasOSSF = 0;
  for (MAuint32 lid = 0; lid < leptons.size(); lid++)
  {
    for (MAuint32 sid = 0; (sid < leptons.size()) && (sid != lid); sid++)
    {
        int sumflvr = leptons[lid]->pdgid() + leptons[sid]->pdgid();
        // OS: 0, 2,  OSSF : 0
        // SSOF : 24, SSSF : 22, 26
        // SSOF : 24, OSOF : 2
        if (std::abs(sumflvr) == 0) {hasOSSF++; hasOS++;}
    }
  }

  if(!Manager()->ApplyCut(leptons.size()==2, "$N_l=2$")) return true;
  if(!Manager()->ApplyCut(hasOSSF>0, "OSSF")) return true;
 
 double pTlepton_Min = 500;
 double pTlepton_Max = 1600;
 double MassZprime = 3000;
 double FB = pTlepton_Min + 0.1*MassZprime;

 for (MAuint32 i=0; i<leptons.size(); i++)
 {   
 Manager()->FillHisto("pT", leptons[i]->pt());

 	if (leptons[i]->pt() > pTlepton_Min and leptons[i]->pt() < FB)
 	{
	        Manager()->FillHisto("AFP_Left", leptons[i]->pt());
	}
	if (leptons[i]->pt() > FB and leptons[i]->pt() < pTlepton_Max)
	{
		Manager()->FillHisto("AFP_Right", leptons[i]->pt());
	}
 }	


  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void AFP::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
}
