#ifndef LJMet_Com_interface_BaseEventSelector_h
#define LJMet_Com_interface_BaseEventSelector_h

/*
 Interface class for FWLite PAT analyzer-selectors
 Specific selectors must implement the () operator
 
 Author: Gena Kukartsev, 2010, 2012
         Orduna@2014
 */

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <utility> // std::pair

#include "FWCore/Framework/interface/Event.h"
#include "LJMet/Com/interface/LjmetEventContent.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "PhysicsTools/SelectorUtils/interface/EventSelector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "TLorentzVector.h"
#include "LJMet/Com/interface/BTagSFUtil.h"
#include "LJMet/Com/interface/BtagHardcodedConditions.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "TMath.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodBDT.h"

//#include "TROOT.h"
//#include "TVector3.h"

struct MVAElectronVars {
    Float_t see, spp, OneMinusE1x5E5x5, R9, etawidth, phiwidth, HoE, PreShowerOverRaw, kfhits, kfchi2, gsfchi2, fbrem, convVtxFitProbability, EoP, eleEoPout, IoEmIoP, deta, dphi, detacalo, gsfhits, expectedMissingInnerHits, pt, isBarrel, isEndcap, SCeta, eClass, pfRelIso, expectedInnerHits, vtxconv, mcEventWeight, mcCBmatchingCategory;
};

class BaseEventSelector : public EventSelector {
    //
    // Base class for all event selector plugins
    //
    
    friend class LjmetFactory;
    
public:
    BaseEventSelector();
    virtual ~BaseEventSelector() { };
    virtual void BeginJob(std::map<std::string, edm::ParameterSet const > par);
    virtual bool operator()( edm::EventBase const & event, pat::strbitset & ret) = 0;
    virtual void EndJob() {

      delete ResJetPar;
      delete L3JetPar;
      delete L2JetPar;
      delete L1JetPar;

      delete ResJetParAK8;
      delete L3JetParAK8;
      delete L2JetParAK8;
      delete L1JetParAK8;

      delete jecUnc;

      delete JetCorrector;
      delete JetCorrectorAK8;

    }
    virtual void AnalyzeEvent( edm::EventBase const & event, LjmetEventContent & ec ) { }
    virtual void JECbyIOV(edm::EventBase const & event);
    std::string GetName() { return mName; }
    /// Evaluates a signed perp components of v1 relative to v2. The sign is defined by Phi
    double GetPerp(TVector3 & v1, TVector3 & v2);
    bool AreMatched(const reco::Candidate & c1, const reco::Candidate & c2, double DR, double DPtRel) { return ((reco::deltaR(c1, c2) < DR) && (abs(c2.pt() - c1.pt())/ c2.pt() < DPtRel)); }
    
    std::vector<edm::Ptr<pat::Jet>> const & GetAllJets() const { return mvAllJets; }
    std::vector<edm::Ptr<pat::Jet>> const & GetSelectedJets() const { return mvSelJets; }
    std::vector<pat::Jet> const & GetSelectedCleanedJets() const { return mvSelJetsCleaned; }
    std::vector<pat::Jet> const & GetSelectedCorrJets() const { return mvSelCorrJets; }
    std::vector<TLorentzVector> const & GetSelectedCorrJets_jesup() const { return mvCorrJets_jesup; }
    std::vector<TLorentzVector> const & GetSelectedCorrJets_jesdn() const { return mvCorrJets_jesdn; }
    std::vector<TLorentzVector> const & GetSelectedCorrJets_jerup() const { return mvCorrJets_jerup; }
    std::vector<TLorentzVector> const & GetSelectedCorrJets_jerdn() const { return mvCorrJets_jerdn; }
    std::vector<edm::Ptr<pat::Jet>> const & GetLooseJets() const { return mvSelJets; }
    std::vector<edm::Ptr<pat::Jet>> const & GetSelectedBtagJets() const { return mvSelBtagJets; }
    std::vector<std::pair<TLorentzVector, bool>> const & GetCorrJetsWithBTags() const { return mvCorrJetsWithBTags; }
    std::vector<edm::Ptr<pat::Muon>> const & GetAllMuons() const { return mvAllMuons; }
    std::vector<edm::Ptr<pat::Muon>> const & GetSelectedMuons() const { return mvSelMuons; }
    std::vector<edm::Ptr<pat::Muon>> const & GetLooseMuons() const { return mvLooseMuons; }
    std::vector<edm::Ptr<pat::Electron>> const & GetAllElectrons() const { return mvAllElectrons; }
    std::vector<edm::Ptr<pat::Electron>> const & GetSelectedElectrons() const { return mvSelElectrons; }
    std::vector<edm::Ptr<pat::Electron>> const & GetLooseElectrons() const { return mvLooseElectrons; }
    edm::Ptr<pat::MET> const & GetMet() const { return mpMet; }
    edm::Ptr<reco::PFMET> const & GetType1CorrMet() const { return mpType1CorrMet; }
    TLorentzVector const & GetCorrectedMet() const { return correctedMET_p4; }
    std::vector<unsigned int> const & GetSelectedTriggers() const { return mvSelTriggers; }
    std::map<std::string, unsigned int> const & GetSelectedTriggersEl() const { return mvSelTriggersEl; }
    std::map<std::string, unsigned int> const & GetSelectedTriggersMu() const { return mvSelTriggersMu; }
    std::map<std::string, unsigned int> const & GetSelectedMCTriggersEl() const { return mvSelMCTriggersEl; }
    std::map<std::string, unsigned int> const & GetSelectedMCTriggersMu() const { return mvSelMCTriggersMu; }
    std::vector<edm::Ptr<reco::Vertex>> const & GetSelectedPVs() const { return mvSelPVs; }
    double const & GetTestValue() const { return mTestValue; }
    void SetMc(bool isMc) { mbIsMc = isMc; }
    bool IsMc() { return mbIsMc; }
    bool const & GetIsTau() { return mbIsTau; }
    
    // LJMET event content setters
    void Init( void );
    void SetEventContent(LjmetEventContent * pEc) { mpEc = pEc; }
    /// Declare a new histogram to be created for the module
    void SetHistogram(std::string name, int nbins, double low, double high) { mpEc->SetHistogram(mName, name, nbins, low, high); }
    void SetHistValue(std::string name, double value) { mpEc->SetHistValue(mName, name, value); }
    void SetTestValue(double & test) { mTestValue = test; }
    
    void SetCorrectedMet(TLorentzVector & met) { correctedMET_p4 = met; }
    void SetCorrJetsWithBTags(std::vector<std::pair<TLorentzVector, bool>> & jets) { mvCorrJetsWithBTags = jets; }
    
    bool isJetTagged(const pat::Jet &jet, edm::EventBase const & event, bool applySF = true, int shiftflag = 0, bool subjetflag = false);
    TLorentzVector correctJetForMet(const pat::Jet & jet, edm::EventBase const & event, unsigned int syst = 0);
    TLorentzVector correctJet(const pat::Jet & jet, edm::EventBase const & event, bool doAK8Corr = false, bool forceCorr = false, unsigned int syst = 0);
    pat::Jet correctJetReturnPatJet(const pat::Jet & jet, edm::EventBase const & event, bool doAK8Corr = false, bool forceCorr = false, unsigned int syst = 0);
    TLorentzVector correctMet(const pat::MET & met, edm::EventBase const & event, unsigned int syst = 0, bool useHF = true);
    TLorentzVector correctMet(const pat::MET & met, edm::EventBase const & event, std::vector<pat::Jet> jets, unsigned int syst = 0, bool useHF = true);
    TLorentzVector correctMet(const pat::MET & met, edm::EventBase const & event, std::vector<edm::Ptr<pat::Jet> > jets, unsigned int syst = 0, bool useHF = true);
    double mvaValue(const pat::Electron & electron, edm::EventBase const & event);
    double mvaValue_alt(const pat::Electron & electron, edm::EventBase const & event);
    
protected:
    std::vector<edm::Ptr<pat::Jet>> mvAllJets;
    std::vector<edm::Ptr<pat::Jet>> mvSelJets;
    std::vector<pat::Jet> mvSelJetsCleaned;
    std::vector<pat::Jet> mvSelCorrJets;
    std::vector<TLorentzVector> mvCorrJets_jesup;
    std::vector<TLorentzVector> mvCorrJets_jesdn;
    std::vector<TLorentzVector> mvCorrJets_jerup;
    std::vector<TLorentzVector> mvCorrJets_jerdn;
    std::vector<edm::Ptr<pat::Jet>> mvLooseJets;
    std::vector<std::pair<TLorentzVector, bool>> mvCorrJetsWithBTags;
    std::vector<edm::Ptr<pat::Jet>> mvSelBtagJets;
    std::vector<edm::Ptr<pat::Muon>> mvAllMuons;
    std::vector<edm::Ptr<pat::Muon>> mvSelMuons;
    std::vector<edm::Ptr<pat::Muon>> mvLooseMuons;
    std::vector<edm::Ptr<pat::Electron>> mvAllElectrons;
    std::vector<edm::Ptr<pat::Electron>> mvSelElectrons;
    std::vector<edm::Ptr<pat::Electron>> mvLooseElectrons;
    edm::Ptr<pat::MET> mpMet;
    edm::Ptr<reco::PFMET> mpType1CorrMet;
    TLorentzVector correctedMET_p4;
    std::vector<unsigned int> mvSelTriggers;
    std::map<std::string, unsigned int> mvSelTriggersEl;
    std::map<std::string, unsigned int> mvSelTriggersMu;
    std::map<std::string, unsigned int> mvSelMCTriggersEl;
    std::map<std::string, unsigned int> mvSelMCTriggersMu;
    std::vector<edm::Ptr<reco::Vertex>> mvSelPVs;
    double mTestValue;
    
    // containers for config parameter values
    std::map<std::string, bool> mbPar;
    std::map<std::string, int> miPar;
    std::map<std::string, double> mdPar;
    std::map<std::string, std::vector<double>> mvdPar;
    std::map<std::string, std::string> msPar;
    std::map<std::string, edm::InputTag> mtPar;
    std::map<std::string, std::vector<std::string>> mvsPar;

    std::string mName;
    std::string mLegend;
    bool mbIsMc;
    bool mbIsTau;
    
private:
    int mNCorrJets;
    int mNBtagSfCorrJets;
    double bTagCut;
    BTagSFUtil mBtagSfUtil;
    BtagHardcodedConditions mBtagCond;
    JME::JetResolution resolution;
    JME::JetResolution resolutionAK8;
    JME::JetResolutionScaleFactor resolution_SF;
    JetCorrectionUncertainty *jecUnc;
    JetCorrectorParameters *L3JetPar;
    JetCorrectorParameters *L2JetPar;
    JetCorrectorParameters *L1JetPar;    				   
    JetCorrectorParameters *L3JetParAK8;
    JetCorrectorParameters *L2JetParAK8;
    JetCorrectorParameters *L1JetParAK8;
    JetCorrectorParameters *ResJetPar; 
    JetCorrectorParameters *ResJetParAK8; 
    FactorizedJetCorrector *JetCorrector;
    FactorizedJetCorrector *JetCorrectorAK8;
    LjmetEventContent * mpEc;
    MVAElectronVars allMVAVars;
    TMVA::Reader tmpTMVAReader_EB;
    TMVA::Reader tmpTMVAReader_EE;
    MVAElectronVars allMVAVars_alt;
    TMVA::Reader tmpTMVAReader_EB_alt;
    TMVA::Reader tmpTMVAReader_EE_alt;
    
    TRandom3 JERrand;
    
    /// Private init method to be called by LjmetFactory when registering the selector
    void init() { mLegend = "[" + mName + "]: "; std::cout << mLegend << "registering " << mName << std::endl; }
    void setName(std::string name) { mName = name; }
    /// Do what any event selector must do before event gets checked
    void BeginEvent(edm::EventBase const & event, LjmetEventContent & ec) { mNCorrJets = 0; mNBtagSfCorrJets = 0; }
    /// Do what any event selector must do after event processing is done, but before event content gets saved to file
    void EndEvent(edm::EventBase const & event, LjmetEventContent & ec) { SetHistValue("nBtagSfCorrections", mNBtagSfCorrJets); }
};

#endif
