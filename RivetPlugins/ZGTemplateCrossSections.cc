// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  /// @brief Measurement of pp->Z+gamma (Z->e+e- or Z->mu+mu-) inclusive (not differential yet) cross sections at sqrt(s) = 13 TeV
  // apply generator level selection
  /// @author Kezhu Guo () <>
  class ZGTemplateCrossSections : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ZGTemplateCrossSections);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // prompt final state electrons
      const PromptFinalState el_pfs = PromptFinalState(Cuts::abspid == PID::ELECTRON);
      declare(el_pfs, "PromptFinalStateElectrons");

      // prompt final state muons
      const PromptFinalState mu_pfs = PromptFinalState(Cuts::abspid == PID::MUON);
      declare(mu_pfs, "PromptFinalStateMuons");

      // prompt final state photons
      const PromptFinalState g_pfs = PromptFinalState(Cuts::abspid == PID::PHOTON);
      declare(g_pfs, "PromptFinalStatePhotons");

      book(_h_l0l1_M, 40, 50, 130);
      book(_h_llg_M, 40,50,200);
      book(_h_l0l1_dr, 30, 0, 6);
      book(_h_l0l1_pt, 40, 0, 100);
      book(_h_l0l1_eta, 40, -2.5, 2.5);
      // book(_h_l0l1_deta, );
      // book(_h_l0l1_dphi, );
      book(_h_l0_pt, 40, 0, 100);
      book(_h_l0_eta, 40,-2.5,2.5);
      book(_h_l1_pt, 40, 0, 100);
      book(_h_l1_eta, 40,-2.5,2.5);
      book(_h_g_pt, 40, 0, 100);
      book(_h_g_eta, 40,-2.5,2.5);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      double pTCut_lead = 30.;
      double pTCut_sub = 20.0; // same for both channel
      double etaCut_lead = 2.5;
      double etaCut_sub = etaCut_lead;

      double pTCut_g = 30.;
      double etaCut_g_03 = 2.5;
      double etaCut_g_20 = 1.566;
      double etaCut_g_01 = 1.4442;

      const PromptFinalState PFS_e = apply<PromptFinalState>(event, "PromptFinalStateElectrons");
      vector<Particle> vec_PFSByPt_e = PFS_e.particlesByPt([&](Particle const& e) { // descending order
        return e.pT() > pTCut_sub && e.abseta() < etaCut_sub;
      });

      const PromptFinalState PFS_m = apply<PromptFinalState>(event, "PromptFinalStateMuons");
      vector<Particle> vec_PFSByPt_m = PFS_m.particlesByPt([&](Particle const& m) {
        return m.pT() > pTCut_sub && m.abseta() < etaCut_sub;
      });

      const PromptFinalState PFS_g = apply<PromptFinalState>(event, "PromptFinalStatePhotons");
      vector<Particle> vec_PFSByPt_g = PFS_g.particlesByPt([&](Particle const& g) {
        return g.pT() > pTCut_g && ((g.abseta() < etaCut_g_03 && g.abseta() > etaCut_g_20) 
                                      || g.abseta() < etaCut_g_01);
      });

      for (int i = 0; i < (int) vec_PFSByPt_g.size() - 1; i++) {
        std::cout << vec_PFSByPt_g[i] << std::endl;
      }


      int i1_e = -1, i2_e = -1, i3_e = -1;
      double max_ll_pt_e = -1;
      double max_g_pt_e = -1;
      int sel_e = findPFSLeptonPair(vec_PFSByPt_e, vec_PFSByPt_g, 
                                    // leptonID, 
                                    pTCut_lead,
                                    i1_e, i2_e, i3_e, max_ll_pt_e, max_g_pt_e);

      int i1_m = -1, i2_m = -1, i3_m = -1;
      double max_ll_pt_m = -1;
      double max_g_pt_m = -1;
      int sel_m = findPFSLeptonPair(vec_PFSByPt_m, vec_PFSByPt_g, 
                                    // leptonID, 
                                    pTCut_lead,
                                    i1_m, i2_m, i3_m, max_ll_pt_m, max_g_pt_m);

      int gen_flvr = -1;
      if (sel_e == 1 && sel_m == 1) { // both channels available
        if (max_ll_pt_m > max_ll_pt_e) {
          gen_flvr = 2;
        } else if (max_ll_pt_m < max_ll_pt_e || max_g_pt_m < max_g_pt_e) {
          gen_flvr = 1;
        } else { // max_ll_pt_m == max_ll_pt_e && max_g_pt_m >= max_g_pt_e; thus favors muon here
          gen_flvr = 2;
        }
      } else if (sel_e == 1) {  // electron channel
          gen_flvr = 1;
      } else if (sel_m == 1) {  // muon channel
          gen_flvr = 2;
      }

      if (gen_flvr == 1) {
        fillHistogram(vec_PFSByPt_e, vec_PFSByPt_g, i1_e, i2_e, i3_e);
      } else if (gen_flvr == 2) {
        fillHistogram(vec_PFSByPt_m, vec_PFSByPt_g, i1_m, i2_m, i3_m);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (auto const& x :
          {_h_l0l1_M, _h_llg_M, _h_l0l1_dr, _h_l0l1_pt, 
            _h_l0l1_eta, //_h_l0l1_deta, _h_l0l1_dphi, 
            _h_l0_pt, _h_l0_eta, _h_l1_pt, _h_l1_eta, 
            _h_g_pt, _h_g_eta}) {
        scale(x, crossSection() / picobarn / sumOfWeights());
      }
    }

    //@}


  private:


    /// @name Histograms
    //@{
    Histo1DPtr _h_l0l1_M;
    Histo1DPtr _h_llg_M;
    Histo1DPtr _h_l0l1_dr;
    Histo1DPtr _h_l0l1_os;
    Histo1DPtr _h_l0l1_pt;
    Histo1DPtr _h_l0l1_eta;
    Histo1DPtr _h_l0l1_deta;
    Histo1DPtr _h_l0l1_dphi;
    Histo1DPtr _h_l0_pt;
    Histo1DPtr _h_l0_eta;
    Histo1DPtr _h_l1_pt;
    Histo1DPtr _h_l1_eta;
    Histo1DPtr _h_g_pt;
    Histo1DPtr _h_g_eta;
    //@}


    int findPFSLeptonPair(vector<Particle>& vec_PFSByPt_l, vector<Particle>& vec_PFSByPt_g, 
                          // int pdgID, 
                          double l0_pt_cut,
                          int& i1_l, int& i2_l, int& i3_l, double& max_ll_pt, double& max_g_pt) {
      int l_lepton = (int) vec_PFSByPt_l.size();
      if (l_lepton < 2 || vec_PFSByPt_g.size() < 1) {
        return 0;   // no pair selected
      }
      if (l_lepton == 2 && vec_PFSByPt_l[0].charge() + vec_PFSByPt_l[1].charge() != 0) {
        return 0;
      }

      // already in descending order

      // why not get two sets (pid and -pid), sort, get first element of each; because there are also other conditions
      for (int i = 0; i < l_lepton - 1; i++) {
        auto & l0 = vec_PFSByPt_l[i];
        for (int j = i + 1; j < l_lepton; j++) {
          auto & l1 = vec_PFSByPt_l[j];
          // std::cout << "(" << l0.vector().pt() << ", " 
          //           << l1.vector().pt() << ")" << "\n";
          if (l0.charge() + l1.charge() != 0) {continue;}
          if (l0.pT() <= l0_pt_cut) {continue;} // already in descending order
          if (deltaR(l0.mom(), l1.mom()) <= 0.4) {continue;}

          auto l0l1 = l0.mom() + l1.mom();
          auto l0l1_M_ = l0l1.mass();
          if (l0l1_M_ <= 50) {continue;}

          for (int k = 0; k < (int) vec_PFSByPt_g.size(); k++) {
            // std::cout << "testing: (" << i << j << k << ")" <<"\n";
            auto const& g = vec_PFSByPt_g[k];
            // if ((l0.vector() + l1.vector() + g.vector()).M() <= 182) {continue;}
            if (deltaR(l0.mom(), g.mom()) > 0.4 && deltaR(l1.mom(), g.mom()) > 0.4) {
              if (l0l1.pt() > max_ll_pt || (l0l1.pt() == max_ll_pt && g.pt() > max_g_pt)) {
                i1_l = i, i2_l = j, i3_l = k;
                // std::cout << "(" << i << j << k << ")" <<"\n";
                max_ll_pt = l0l1.pT();
                max_g_pt = g.pT();
              }
            }
          } // end k
        } // end j
      } // end i

      if (i1_l == -1 && i2_l == -1) {
        return 0;
      }
      return 1;

    }


    void fillHistogram (vector<Particle>& selected_leptons, vector<Particle>& selected_photons, 
                        int i1_l, int i2_l, int i3_l) {

      auto const& l0 = selected_leptons[i1_l];
      auto const& l1 = selected_leptons[i2_l];
      auto const& g = selected_photons[i3_l];
      auto l0l1 = l0.mom() + l1.mom();
      // auto gen_l0l1_M_ = l0l1.mass();
      // auto gen_llg_M_ = (l0.mom() + l1.mom() + g.mom()).mass();
      // auto gen_l0l1_dr_ = deltaR(&l0, &l1);
      // auto gen_l0l1_os_ = l0.charge() != l1.charge();
      // auto gen_l0l1_pt_ = l0l1.pT();
      // auto gen_l0l1_eta_ = l0l1.eta();
      // auto gen_l0l1_deta_ = abs(l0.eta() - l1.eta());
      double dphi = abs(l0.phi()-l1.phi());
      if (dphi > PI) {dphi = 2*PI-dphi;}
      // auto gen_l0l1_dphi_ = dphi;
      // auto gen_l0_pt_ = l0.pT();
      // auto gen_l0_eta_ = l0.eta();
      // auto gen_l1_pt_ = l1.pT();
      // auto gen_l1_eta_ = l1.eta();
      // auto gen_g_pt_ = g.pT();
      // auto gen_g_eta_ = g.eta();
     
      _h_l0l1_M->fill(l0l1.mass()/GeV);
      _h_llg_M->fill((l0.mom() + l1.mom() + g.mom()).mass()/GeV);
      _h_l0l1_dr->fill(deltaR(l0.mom(), l1.mom()));
      _h_l0l1_pt->fill(l0l1.pT()/GeV);
      _h_l0l1_eta->fill(l0l1.eta());
      // _h_l0l1_deta->fill(abs(l0.eta() - l1.eta()));
      // _h_l0l1_dphi->fill(dphi);
      _h_l0_pt->fill(l0.pT()/GeV);
      _h_l0_eta->fill(l0.eta());
      _h_l1_pt->fill(l1.pT()/GeV);
      _h_l1_eta->fill(l1.eta());
      _h_g_pt->fill(g.pT()/GeV);
      _h_g_eta->fill(g.eta());

    }

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ZGTemplateCrossSections);


}
