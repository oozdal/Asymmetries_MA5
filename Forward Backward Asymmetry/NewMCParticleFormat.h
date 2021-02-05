#ifndef NEWMCPART_H
#define NEWMCPART_H

// STL headers
#include <algorithm>
#include <cmath>

#include "SampleAnalyzer/Commons/DataFormat/SampleFormat.h"
#include "SampleAnalyzer/Commons/Service/PDGService.h"


namespace MA5
{
  class NewMCParticleFormat
  {
    private:

        /// PDG numbering of the particle
        MAint32 pdgid_;
        MALorentzVector momentum_;
        MAbool charge_;       /// charge of the particle 0 = -1, 1 = +1

    public :
        /// Constructor without argument
        NewMCParticleFormat()
        {
            momentum_.SetPxPyPzE(0.,0.,0.,0.);
            pdgid_      = 0; 
        }

        /// Destructor
        virtual ~NewMCParticleFormat() {}

        /// Clear all information
        virtual void Reset()
        {
          momentum_.SetPxPyPzE(0.,0.,0.,0.);
          pdgid_      = 0; 
        }

        const MAint32& pdgid()  const {return pdgid_;}
        void setMomentum(const MALorentzVector& v)  {momentum_=v;}
        void setPdgid(MAint32 v)   {pdgid_=v;}

        /// Accessor to 4-vector momentum (read-only)
        const MALorentzVector& momentum() const {return momentum_;}

        /// Accessor to 4-vector momentum
        MALorentzVector& momentum() {return momentum_;}

//	MAint32& charge() {return PDG->GetCharge(pdgid_);}

        /// Accessor to the electric charge
  	virtual const MAint32 charge() const
  	{ if (charge_) return +1; else return -1; }

        /// Mutator related to the electric charge
        virtual void SetCharge(MAint32 charge)
        { if (charge>0) charge_=true; else charge_=false; }
	
	MAdouble64 px() {return momentum_.Px();}
        MAdouble64 py() {return momentum_.Py();}
        MAdouble64 pz() {return momentum_.Pz();}

  };
}
#endif
