#include "SampleAnalyzer/User/Analyzer/AFB.h"
#include <cmath>
#include <iostream>
#include "SampleAnalyzer/Commons/Service/PDGService.h"

using namespace MA5;
using namespace std;

MAbool AFB::Initialize(const MA5::Configuration& cfg,
                      const std::map<std::string,std::string>& parameters)
{
  // Initializing PhysicsService for MC
  PHYSICS->mcConfig().Reset();

  // definition of the multiparticle "invisible"
  PHYSICS->mcConfig().AddInvisibleId(12);

  // ===== Signal region ===== //
  Manager()->AddRegionSelection("Rap_No");
  Manager()->AddRegionSelection("Rap_00");
  Manager()->AddRegionSelection("Rap_01");
  Manager()->AddRegionSelection("Rap_02");  
  Manager()->AddRegionSelection("Rap_04");
  Manager()->AddRegionSelection("Rap_06");

  // ===== Selections ===== //
  Manager()->AddCut("$N_l=2$");
  Manager()->AddCut("OSSF");
  Manager()->AddCut("$|y_{ll}|>0.0$", "Rap_00");
  Manager()->AddCut("$|y_{ll}|>0.1$", "Rap_01");  
  Manager()->AddCut("$|y_{ll}|>0.2$", "Rap_02");  
  Manager()->AddCut("$|y_{ll}|>0.4$", "Rap_04");
  Manager()->AddCut("$|y_{ll}|>0.6$", "Rap_06");

  // ===== Histograms ===== //

//  Manager()->AddHisto("AFB_pos_Rap_No",   30,50.,200., "Rap_No");
//  Manager()->AddHisto("AFB_neg_Rap_No",   30,50.,200., "Rap_No");

//  Manager()->AddHisto("AFB_pos_Rap_00",   15,50.,200., "Rap_00");
//  Manager()->AddHisto("AFB_neg_Rap_00",   15,50.,200., "Rap_00");

//  Manager()->AddHisto("AFB_pos_Rap_01",  15,50.,200., "Rap_01");
//  Manager()->AddHisto("AFB_neg_Rap_01",  15,50.,200., "Rap_01");

//  Manager()->AddHisto("AFB_pos_Rap_02",  15,50.,200., "Rap_02");
//  Manager()->AddHisto("AFB_neg_Rap_02",  15,50.,200., "Rap_02");

//  Manager()->AddHisto("AFB_pos_Rap_04",  15,50.,200., "Rap_04");
//  Manager()->AddHisto("AFB_neg_Rap_04",  15,50.,200., "Rap_04");

//  Manager()->AddHisto("AFB_pos_Rap_06",  15,50.,200., "Rap_06");
//  Manager()->AddHisto("AFB_neg_Rap_06",  15,50.,200., "Rap_06");


  Manager()->AddHisto("AFB_pos_Rap_No",  13,2400.,3600., "Rap_No");
  Manager()->AddHisto("AFB_neg_Rap_No",  13,2400.,3600., "Rap_No");  

  Manager()->AddHisto("AFB_pos_Rap_00",  13,2400.,3600., "Rap_00");
  Manager()->AddHisto("AFB_neg_Rap_00",  13,2400.,3600., "Rap_00");

  Manager()->AddHisto("AFB_pos_Rap_01",  13,2400.,3600., "Rap_01");
  Manager()->AddHisto("AFB_neg_Rap_01",  13,2400.,3600., "Rap_01");  

  Manager()->AddHisto("AFB_pos_Rap_02",  13,2400.,3600., "Rap_02");
  Manager()->AddHisto("AFB_neg_Rap_02",  13,2400.,3600., "Rap_02");

  Manager()->AddHisto("AFB_pos_Rap_04",  6,2400.,3600., "Rap_04");
  Manager()->AddHisto("AFB_neg_Rap_04",  6,2400.,3600., "Rap_04");

  Manager()->AddHisto("AFB_pos_Rap_06",  6,2400.,3600., "Rap_06");
  Manager()->AddHisto("AFB_neg_Rap_06",  6,2400.,3600., "Rap_06");


  Manager()->AddHisto("Eta_lp_Rap_No",  16,-4.,4., "Rap_No");
  Manager()->AddHisto("Eta_lm_Rap_No",  16,-4.,4., "Rap_No");
  
  Manager()->AddHisto("Eta_lp_Rap_00",  16,-4.,4., "Rap_00");
  Manager()->AddHisto("Eta_lm_Rap_00",  16,-4.,4., "Rap_00");

  Manager()->AddHisto("Eta_lp_Rap_01",  16,-4.,4., "Rap_01");
  Manager()->AddHisto("Eta_lm_Rap_01",  16,-4.,4., "Rap_01");

  Manager()->AddHisto("Eta_lp_Rap_02",  16,-4.,4., "Rap_02");
  Manager()->AddHisto("Eta_lm_Rap_02",  16,-4.,4., "Rap_02");  

  Manager()->AddHisto("Eta_lp_Rap_04",  16,-4.,4., "Rap_04");
  Manager()->AddHisto("Eta_lm_Rap_04",  16,-4.,4., "Rap_04");  

  Manager()->AddHisto("Eta_lp_Rap_06",  16,-4.,4., "Rap_06");
  Manager()->AddHisto("Eta_lm_Rap_06",  16,-4.,4., "Rap_06");  

//  Manager()->AddHisto("Mll_Rap_No",   15,50.,200., "Rap_No");
//  Manager()->AddHisto("Mll_Rap_00",   15,50.,200., "Rap_00");
//  Manager()->AddHisto("Mll_Rap_01",   15,50.,200., "Rap_01");
//  Manager()->AddHisto("Mll_Rap_02",   15,50.,200., "Rap_02");
//  Manager()->AddHisto("Mll_Rap_04",   15,50.,200., "Rap_04");
//  Manager()->AddHisto("Mll_Rap_06",   15,50.,200., "Rap_06");

//  Manager()->AddHisto("Mll_Rap_No",   13,2400.,3600., "Rap_No");  
//  Manager()->AddHisto("Mll_Rap_00",   13,2400.,3600., "Rap_00");
//  Manager()->AddHisto("Mll_Rap_01",   13,2400.,3600., "Rap_01");  
//  Manager()->AddHisto("Mll_Rap_02",   13,2400.,3600., "Rap_02");  
//  Manager()->AddHisto("Mll_Rap_04",    6,2400.,3600., "Rap_04");
//  Manager()->AddHisto("Mll_Rap_06",    6,2400.,3600., "Rap_06");


  // No problem during initialization

  myfile.open("mydata.csv");//, std::ios_base::app);
  myfile << "Mll" << ","
         << "Eta_lp" << ","
	 << "Eta_lm" << ","
	 << "Eta_l1" << ","
	 << "Eta_l2" << ","
         << "EvWeight" << "\n";

  return true;
}

MAbool AFB::Execute(SampleFormat& sample, const EventFormat& event)
{
  MAfloat32 EvWeight = 1.0;
  if (weighted_events_ && event.mc()!=0) EvWeight *= event.mc()->weight();

  //EvWeight *= sample.mc()->xsection(); //* 1000. * lumi;
  if (sample.mc()!=0) sample.mc()->addWeightedEvents(EvWeight);
  Manager()->InitializeForNewEvent(EvWeight);


  // Clearing particle containers
  leptons.clear();
  jets.clear();
  initialquarks.clear();
//  lepton1.clear();
//  lepton2.clear();

  // Filling particle containers
  for (MAuint32 i=0;i<event.mc()->particles().size();i++)
  {	  
    if (is_lepton(&(event.mc()->particles()[i])))
      leptons.push_back(&(event.mc()->particles()[i]));
    if (is_jet(&(event.mc()->particles()[i])))
      jets.push_back(&(event.mc()->particles()[i]));
    if (is_initialquark(&(event.mc()->particles()[i])))
      initialquarks.push_back(&(event.mc()->particles()[i])); 
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
  
  // Useful quantities
  double mll    = (leptons[0]->momentum() + leptons[1]->momentum()).M();
  double eta_l1 = leptons[0]->eta();
  double eta_l2 = leptons[1]->eta();
//  double angle_lab  = 0;
  MAdouble64 angle_CM_true = 0;
  MAdouble64 angle_CM_rec  = 0;  
  MAdouble64 cos_angle_CM_star = 0;
  double etalp  = 0;
  double etalm  = 0;

  MALorentzVector lplepton;
  MALorentzVector lmlepton;

  for (MAuint32 i=0; i<leptons.size(); i++)
  {
  	if (leptons[i]->pdgid()<0)
	{
		lplepton.SetPxPyPzE(leptons[i]->px(),leptons[i]->py(),leptons[i]->pz(),leptons[i]->e());    
	}	
        if (leptons[i]->pdgid()>0)
        {
                lmlepton.SetPxPyPzE(leptons[i]->px(),leptons[i]->py(),leptons[i]->pz(),leptons[i]->e());
        } 	
  }

 MALorentzVector lplm_true = lplepton + lmlepton;
 MABoost recoil_true;
 recoil_true.setBoostVector(-lplm_true.Px()/lplm_true.E(), -lplm_true.Py()/lplm_true.E(), -lplm_true.Pz()/lplm_true.E());
 MALorentzVector boosted_lp = lplepton;
 MALorentzVector boosted_lm = lmlepton; 
// recoil_true.boost(boosted_lp);
 recoil_true.boost(boosted_lm);

 double mll_boosted = (boosted_lp + boosted_lm).M();
 double eta_boosted_lp = 0;
 double eta_boosted_lm = 0;
 double quarkcharge  = 0;
 double leptoncharge = 0;

 MALorentzVector pos_quarks;
 MALorentzVector neg_quarks;

 MAdouble64 yll_lab    = (leptons[0]->momentum() + leptons[1]->momentum()).Rapidity();
 MAdouble64 yll_CM     = (boosted_lp + boosted_lm).Rapidity();

 signed int signofyll_lab = yll_lab/std::fabs(yll_lab);
 signed int signofyll_CM  = yll_CM/std::fabs(yll_CM);

 if(!Manager()->ApplyCut(std::fabs(yll_lab>0.0), "$|y_{ll}|>0.0$")) return true; 
 if(!Manager()->ApplyCut(std::fabs(yll_lab>0.1), "$|y_{ll}|>0.1$")) return true;
 if(!Manager()->ApplyCut(std::fabs(yll_lab>0.2), "$|y_{ll}|>0.2$")) return true; 
 if(!Manager()->ApplyCut(std::fabs(yll_lab>0.4), "$|y_{ll}|>0.4$")) return true;
 if(!Manager()->ApplyCut(std::fabs(yll_lab>0.6), "$|y_{ll}|>0.6$")) return true;

 for (MAuint32 i=0; i<initialquarks.size(); i++)
 {
	quarkcharge = (PDG->GetCharge(initialquarks[i]->pdgid()))/3.; 
	if ( initialquarks[i]->pdgid() > 0 )
	{	
		pos_quarks.SetPxPyPzE(initialquarks[i]->px(),initialquarks[i]->py(),initialquarks[i]->pz(),initialquarks[i]->e());
	}
        if ( initialquarks[i]->pdgid() < 0 )
	{
                neg_quarks.SetPxPyPzE(initialquarks[i]->px(),initialquarks[i]->py(),initialquarks[i]->pz(),initialquarks[i]->e());
       	}	
 }

 for (MAuint32 i=0; i<leptons.size(); i++)
 {      
	leptoncharge =  (PDG->GetCharge(leptons[i]->pdgid()))/3.; 
        if ( leptoncharge > 0. )
        {
//      angle_lab = leptons[i]->angle(pos_quarks);
//	angle_CM_true  = boosted_lp.Angle(pos_quarks);
//	angle_CM_rec   = boosted_lp.Angle(boosted_lp + boosted_lm);
//	cos_angle_CM_star = cos(angle_CM_rec)*( yll_lab/std::fabs(yll_lab) );
        etalp = leptons[i]->eta();
//      eta_boosted_lp = boosted_lp.Eta();	

        Manager()->FillHisto("Eta_lp_Rap_No", etalp);
        Manager()->FillHisto("Eta_lp_Rap_00", etalp);
        Manager()->FillHisto("Eta_lp_Rap_01", etalp);	
        Manager()->FillHisto("Eta_lp_Rap_02", etalp);
        Manager()->FillHisto("Eta_lp_Rap_04", etalp);
        Manager()->FillHisto("Eta_lp_Rap_06", etalp);

//	Manager()->FillHisto("Eta_lp_Rap_No", eta_boosted_lp);
//      Manager()->FillHisto("Eta_lp_Rap_00", eta_boosted_lp);
//      Manager()->FillHisto("Eta_lp_Rap_02", eta_boosted_lp);	
//      Manager()->FillHisto("Eta_lp_Rap_04", eta_boosted_lp);
//      Manager()->FillHisto("Eta_lp_Rap_06", eta_boosted_lp);
	}

	if ( leptoncharge < 0. )
        {
//      angle_lab = leptons[i]->angle(pos_quarks);
        angle_CM_true  = boosted_lm.Angle(pos_quarks);
        angle_CM_rec   = boosted_lm.Angle(boosted_lp + boosted_lm);
//	angle_lab_rec   = boosted_lm.Angle(lplepton + lmlepton);
	cos_angle_CM_star = cos(angle_CM_rec)*( yll_CM/std::fabs(yll_CM) );
        etalm = leptons[i]->eta();
        eta_boosted_lm = boosted_lm.Eta();
 
//	double test_eta = -log(tan(angle_CM_true/2.)); 
//      cout << angle_CM_true << " " << test_eta  << " " << etalm << " " << eta_boosted_lm <<"\n";

//      double test_eta = -log(tan(angle_CM_rec/2.));
//      cout << cos(angle_CM_rec) << " " << test_eta  << " " << etalm << " " << eta_boosted_lm <<"\n";

//      double test_eta = -log(tan(cos_angle_CM_star/2.));
//      cout << cos(angle_CM_rec) << " " << test_eta  << " " << etalm << " " << eta_boosted_lm <<"\n";

        Manager()->FillHisto("Eta_lm_Rap_No", etalm);
        Manager()->FillHisto("Eta_lm_Rap_00", etalm);
        Manager()->FillHisto("Eta_lm_Rap_01", etalm);	
        Manager()->FillHisto("Eta_lm_Rap_02", etalm);
        Manager()->FillHisto("Eta_lm_Rap_04", etalm);
        Manager()->FillHisto("Eta_lm_Rap_06", etalm);

//      Manager()->FillHisto("Eta_lm_Rap_No", eta_boosted_lm);
//      Manager()->FillHisto("Eta_lm_Rap_00", eta_boosted_lm);
//      Manager()->FillHisto("Eta_lm_Rap_02", eta_boosted_lm);	
//      Manager()->FillHisto("Eta_lm_Rap_04", eta_boosted_lm);
//      Manager()->FillHisto("Eta_lm_Rap_06", eta_boosted_lm);
        }

        if ( leptoncharge < 0 and cos(angle_CM_true)>0 )
		
        {
        Manager()->FillHisto("AFB_pos_Rap_No", mll);
	}
        
	if ( leptoncharge < 0 and cos(angle_CM_true)<0 )
        {
        Manager()->FillHisto("AFB_neg_Rap_No", mll);
	}
  }	

//if ( (signofyll_lab>0 and etalm<0) or (signofyll_lab<0 and etalm>0) ) 
  if ( (signofyll_lab>0 and eta_boosted_lm>0) or (signofyll_lab<0 and eta_boosted_lm<0) )
  {
  Manager()->FillHisto("AFB_pos_Rap_00", mll);
  Manager()->FillHisto("AFB_pos_Rap_01", mll);
  Manager()->FillHisto("AFB_pos_Rap_02", mll);
  Manager()->FillHisto("AFB_pos_Rap_04", mll);
  Manager()->FillHisto("AFB_pos_Rap_06", mll);
  }

//if ( (signofyll_lab>0 and etalm>0) or (signofyll_lab<0 and etalm<0) )
  if ( (signofyll_lab>0 and eta_boosted_lm<0) or (signofyll_lab<0 and eta_boosted_lm>0) )	  
  {
  Manager()->FillHisto("AFB_neg_Rap_00", mll);
  Manager()->FillHisto("AFB_neg_Rap_01", mll);
  Manager()->FillHisto("AFB_neg_Rap_02", mll);
  Manager()->FillHisto("AFB_neg_Rap_04", mll);
  Manager()->FillHisto("AFB_neg_Rap_06", mll);
  }


//  Manager()->FillHisto("Mll_Rap_No",  mll);
//  Manager()->FillHisto("Mll_Rap_00",  mll);
//  Manager()->FillHisto("Mll_Rap_01",  mll);  
//  Manager()->FillHisto("Mll_Rap_02",  mll);  
//  Manager()->FillHisto("Mll_Rap_04",  mll);
//  Manager()->FillHisto("Mll_Rap_06",  mll);
    
  myfile  << mll << ","
	  << etalp << ","
	  << etalm << ","
	  << eta_l1 << ","
	  << eta_l2 << ","
          << EvWeight << "\n";

  return true;
}

void AFB::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
	myfile.close();
}
