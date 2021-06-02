#include "SampleAnalyzer/User/Analyzer/ACC.h"
using namespace MA5;
using namespace std;

MAbool ACC::Initialize(const MA5::Configuration& cfg,
                      const std::map<std::string,std::string>& parameters)
{
  // Initializing PhysicsService for MC
  PHYSICS->mcConfig().Reset();

  // definition of the multiparticle "invisible"
  PHYSICS->mcConfig().AddInvisibleId(12);
  PHYSICS->mcConfig().AddInvisibleId(14);
  PHYSICS->mcConfig().AddInvisibleId(-12);
  PHYSICS->mcConfig().AddInvisibleId(-14);
//  PHYSICS->mcConfig().AddInvisibleId(9900012);
//  PHYSICS->mcConfig().AddInvisibleId(-9900012);

  // ===== Signal region ===== //
  Manager()->AddRegionSelection("SR_Lepton");

  // ===== Selections ===== //
  Manager()->AddCut("One Lepton");  
//  Manager()->AddCut("MtMiss > 2000");

  // ===== Histograms ===== //
  // ===== AbsEta ===== //
  Manager()->AddHisto("AbsEta(Lepton_plus)", 20, 0., 2.5);
  Manager()->AddHisto("AbsEta(Lepton_minus)", 20, 0., 2.5);
 
  Manager()->AddHisto("MET One Lepton", 25,0.,2500.);  
  Manager()->AddHisto("MtMiss One Lepton", 15,0.,3000.);

  // No problem during initialization

  ACCfile.open("ACC.csv");//, std::ios_base::app);
  ACCfile << "MET"    << ","
       	  << "MtMiss" << ","
          << "Eta" << ","
	  << "Lepton_PT" << ","
	  << "Lepton_ID" << ","
          << "BoostedLep_Eta" << ","
	  << "EvWeight" << "\n";
  
  return true;
}


MAbool ACC::Execute(SampleFormat& sample, const EventFormat& event)
{
  MAfloat32 __event_weight__ = 1.0;
  if (weighted_events_ && event.mc()!=0) __event_weight__ = event.mc()->weight();

  if (sample.mc()!=0) sample.mc()->addWeightedEvents(__event_weight__);
  Manager()->InitializeForNewEvent(__event_weight__);

  // Clearing particle containers
  leptons.clear();

  // Filling particle containers
  for (MAuint32 i=0;i<event.mc()->particles().size();i++)
  {
    if (is_lepton(&(event.mc()->particles()[i])))
      leptons.push_back(&(event.mc()->particles()[i]));
  }

  // Sorting particles
  SORTER->sort(leptons,PTordering);

  MAbool cutPosLep=false;
  MAbool cutNegLep=false;
  MAbool countLepton=false;
  if(leptons.size()==1)
  {
   MAbool cutPosLepton = ( leptons[0]->pdgid() < 0);
   MAbool cutNegLepton = ( leptons[0]->pdgid() > 0);
   cutPosLep = cutPosLepton;
   cutNegLep = cutNegLepton;
   countLepton=true;
  } 

  if(!Manager()->ApplyCut(countLepton, "One Lepton")) return true;

  // Useful quantities 
  // MET & MT
  double MET_double = PHYSICS->Transverse->EventMET(event.mc());

  MALorentzVector MET = event.mc()->MET().momentum();
  double MtMiss = leptons[0]->mt_met(MET);

 MALorentzVector lep_vec;

  for (MAuint32 i=0; i<leptons.size(); i++)
  {
      lep_vec.SetPxPyPzE(leptons[i]->px(),leptons[i]->py(),leptons[i]->pz(),leptons[i]->e());
  }

 MALorentzVector lep_true = lep_vec + MET;
 MABoost recoil_true;
 recoil_true.setBoostVector(-lep_true.Px()/lep_true.E(), -lep_true.Py()/lep_true.E(), -lep_true.Pz()/lep_true.E());
 MALorentzVector boosted_lep = lep_vec;
 recoil_true.boost(boosted_lep);

//  double MtMiss = 0;
//  if ( leptons.size() == 1 ) { MtMiss = sqrt(2*leptons[0]->pt()*MET*(1-cos(leptons[0]->dphi_0_pi(event.mc()->MET())))); }

//  cout << MET_double << " " << MtMiss << "\n";

//  if(!Manager()->ApplyCut(MtMiss > 2000,"MtMiss > 2000")) return true;


  
// Filling Eta
  // Computing the ACC components for AbsEta
  if(cutPosLep) Manager()->FillHisto("AbsEta(Lepton_plus)",  leptons[0]->eta());
  if(cutNegLep) Manager()->FillHisto("AbsEta(Lepton_minus)", leptons[0]->eta());

  if(countLepton)  Manager()->FillHisto("MET One Lepton", MET_double);
  if(countLepton)  Manager()->FillHisto("MtMiss One Lepton", MtMiss);

  
  double Eta = leptons[0]->eta();
  double Lepton_ID = leptons[0]->pdgid();
  double BoostedLep_Eta = boosted_lep.Eta();
  double Lep_PT = leptons[0]->pt(); 
  
  ACCfile  << MET_double << ","
           << MtMiss << ","
           << Eta << ","
	   << Lep_PT  << ","
	   << Lepton_ID << ","
	   << BoostedLep_Eta << ","
           << __event_weight__ << "\n";

  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void ACC::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
	ACCfile.close();
}

