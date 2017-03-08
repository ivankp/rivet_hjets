// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

#define NJETS 2

#define h_(NAME) Histo1DPtr _h_##NAME;
#define BOOK(NAME,nbins,min,max) _h_##NAME = bookHisto1D(#NAME,nbins,min,max);

namespace Rivet {
  /// @brief Add a short analysis description here
  class hjets : public Analysis {
  private:
    /// @name Histograms
    //@{
    CounterPtr _h_xsec;

    h_(jets_N_excl)

    h_(H_pT)
    h_(H_y)

    h_(HT)
    h_(Hjets_mass)

#if NJETS >= 1
    h_(j1_pT)
    h_(j1_y)
#elif NJETS >= 2
    h_(j2_pT)
    h_(j2_y)
#elif NJETS >= 3
    h_(j3_pT)
    h_(j3_y)
#endif
    //@}
  public:
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(hjets);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const FinalState fs;
      IdentifiedFinalState fs_higgs(PID::HIGGS);
      VetoedFinalState fs_rest(fs);
      fs_rest.addVetoOnThisFinalState(fs_higgs);
      addProjection(fs, "FS");
      addProjection(fs_higgs, "Higgs");
      addProjection(fs_rest, "Rest");
      addProjection(FastJets(fs_rest, FastJets::ANTIKT, 0.4), "Jets");

      // Book histograms
      _h_xsec = bookCounter(3, 1, 1);

      BOOK(jets_N_excl,NJETS+1,0,NJETS+1)

      BOOK(H_pT,100,0,1000)
      BOOK(H_y,90,-4.5,4.5)

      BOOK(HT,150,0,1500)
      BOOK(Hjets_mass,200,100,2100)

#if NJETS >= 1
      BOOK(j1_pT,70,0,700)
      BOOK(j1_y,90,-4.5,4.5)
#elif NJETS >= 2
      BOOK(j2_pT,70,0,700)
      BOOK(j2_y,90,-4.5,4.5)
#elif NJETS >= 3
      BOOK(j3_pT,70,0,700)
      BOOK(j3_y,90,-4.5,4.5)
#endif
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // https://rivet.hepforge.org/trac/wiki/WritingAnAnalysis
      const double weight = event.weight();

      const auto higgs = applyProjection<IdentifiedFinalState>(event,"Higgs")
        .particles();
      const auto jets  = applyProjection<FastJets>(event,"Jets")
        .jetsByPt(Cuts::pT > 30.*GeV && Cuts::abseta < 4.4);

      if (higgs.size()!=1) return;

      _h_jets_N_excl->fill(higgs.size());

      if (jets.size() < NJETS) return;

      const auto& hmom = higgs[0].momentum();
      const double H_pT = hmom.pT();
      _h_H_pT->fill(H_pT,weight);
      _h_H_y ->fill(hmom.rap(),weight);

      double HT = H_pT;
      auto Hjets = hmom;
      for (const auto& j : jets) {
        HT += j.pT();
        Hjets += j;
      }

      _h_HT->fill(HT);
      _h_Hjets_mass->fill(Hjets.mass());

#if NJETS >= 1
      _h_j1_pT->fill(jets[0].pT(),weight);
      _h_j1_y ->fill(jets[0].rap(),weight);
#elif NJETS >= 2
      _h_j2_pT->fill(jets[1].pT(),weight);
      _h_j2_y ->fill(jets[1].rap(),weight);
#elif NJETS >= 3
      _h_j3_pT->fill(jets[2].pT(),weight);
      _h_j3_y ->fill(jets[2].rap(),weight);
#endif
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_xsec, crossSection()/picobarn/sumOfWeights()); // norm to cross section
    }
    //@}
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(hjets);
}
