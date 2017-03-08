// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {
  /// @brief Add a short analysis description here
  class hjets : public Analysis {
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
      _h_H_pT = bookHisto1D("H_pT",100,0,2e3);
      _h_xsec = bookCounter(3, 1, 1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // https://rivet.hepforge.org/trac/wiki/WritingAnAnalysis
      const double weight = event.weight();

      const auto higgs = applyProjection<IdentifiedFinalState>(event,"Higgs")
        .particles();
      const auto jets  = applyProjection<FastJets>(event,"Jets")
        .jetsByPt(Cuts::pT > 30*GeV && Cuts::abseta < 4.4);

      MSG_DEBUG("Jet multiplicity = " << jets.size());

      if (higgs.size()==1) {
        const auto hmom = higgs[0].momentum();
        _h_H_pT->fill(hmom.pT(),weight);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize(_h_YYYY); // normalize to unity
      scale(_h_xsec, crossSection()/picobarn/sumOfWeights()); // norm to cross section
    }
    //@}

  private:
    /// @name Histograms
    //@{
    Histo1DPtr _h_H_pT;
    CounterPtr _h_xsec;
    //@}
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(hjets);
}
