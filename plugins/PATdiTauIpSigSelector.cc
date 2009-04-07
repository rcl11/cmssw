#include "TauAnalysis/CandidateTools/plugins/PATdiTauIpSigSelector.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

template <typename T1, typename T2, typename TExtr1>
PATdiTauIpSigSelector<T1,T2,TExtr1>::PATdiTauIpSigSelector(const edm::ParameterSet& cfg) 
{
  vertexSrc_ = cfg.getParameter<edm::InputTag>("vertexSource");

  if ( cfg.exists("IpMin") ) {
    ipMin_ = cfg.getParameter<double>("IpMin");
    applyIpMin_ = true;
  } else {
    applyIpMin_ = false;
  }

  if ( cfg.exists("IpMax") ) {
    ipMax_ = cfg.getParameter<double>("IpMax");
    applyIpMax_ = true;
  } else {
    applyIpMax_ = false;
  }
}

template <typename T1, typename T2, typename TExtr>
PATdiTauIpSigSelector<T1,T2,TExtr>::~PATdiTauIpSigSelector() 
{
//--- nothing to be done yet...
}

template <typename T1, typename T2, typename TExtr>
void PATdiTauIpSigSelector<T1,T2,TExtr>::select(const edm::Handle<collection>& patDiTauCollection,
					   edm::Event& evt, const edm::EventSetup& es) 
{
  selected_.clear();
  
  edm::Handle<reco::VertexCollection> primaryEventVertexCollection;
  evt.getByLabel(vertexSrc_, primaryEventVertexCollection);
  if ( primaryEventVertexCollection->size() < 1 ) {
    return;
  } else if ( primaryEventVertexCollection->size() > 1 ) {
    edm::LogError ("PATdiTauIpSigSelector::select") << " Cannot have more than one primary event vertex --> skipping !!";
    return;
  } 

  const reco::Vertex& thePrimaryEventVertex = (*primaryEventVertexCollection->begin());

  for ( typename collection::const_iterator patDiTau = patDiTauCollection->begin(); 
	patDiTau != patDiTauCollection->end(); ++patDiTau ) {
    
    bool ipIsValid = false;
    double ip = ipSigExtractor_(*patDiTau, thePrimaryEventVertex, ipIsValid);
    if ( ipIsValid && (!applyIpMin_ || ip > ipMin_) &&
	 (!applyIpMax_ || ip < ipMax_) ) {
      selected_.push_back(&(*patDiTau));
    }
  }
}

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

struct PATElectronIpSigExtractor
{
  double operator()(const pat::Electron& electron, const reco::Vertex& primaryEventVertex, bool& isValid) 
  {
    if ( electron.gsfTrack().isNonnull() ) {
      isValid = true;
      return (electron.gsfTrack()->dxy(primaryEventVertex.position())/electron.gsfTrack()->dxyError());
    } else {
      isValid = false;
      return 0.;
    }
  }
};

struct PATMuonIpSigExtractor
{
  double operator()(const pat::Muon& muon, const reco::Vertex& primaryEventVertex, bool& isValid) 
  {
    if ( muon.track().isNonnull() ) {
      isValid = true;
      return (muon.track()->dxy(primaryEventVertex.position())/muon.track()->dxyError());
    } else {
      isValid = false;
      return 0.;
    }
  }
};

struct PATTauIpSigExtractor
{
  double operator()(const pat::Tau& tau, const reco::Vertex& primaryEventVertex, bool& isValid) 
  {
    if ( tau.leadTrack().isNonnull() ) {
      isValid = true;
      return (tau.leadTrack()->dxy(primaryEventVertex.position())/tau.leadTrack()->dxyError());
    } else {
      isValid = false;
      return 0.;
    }
  }
};

template <typename T1, typename T2, typename Extr1, typename Extr2>
struct PATdiTauIpSigExtractor
{
  double operator()(const CompositePtrCandidateT1T2MEt<T1,T2>& diTau, const reco::Vertex& primaryEventVertex, bool& isValid)
  {
    Extr1 extr1_;
    Extr2 extr2_;
    bool valid1 = false;
    bool valid2 = false;
    double sig1 = extr1_(*(diTau.leg1()),primaryEventVertex,valid1);
    double sig2 = extr2_(*(diTau.leg2()),primaryEventVertex,valid2);
    if (valid1&&valid2) {
      isValid = true;
      return (sqrt(sig1*sig1+sig2*sig2));
    } else {
      isValid = false;
      return 999;
    }
  }
};

#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"

typedef PATdiTauIpSigExtractor<pat::Electron, pat::Muon, PATElectronIpSigExtractor, PATMuonIpSigExtractor> PATElecMuIpSigExtractor;
typedef PATdiTauIpSigExtractor<pat::Electron, pat::Tau, PATElectronIpSigExtractor, PATTauIpSigExtractor> PATElecTauIpSigExtractor;
typedef PATdiTauIpSigExtractor<pat::Muon, pat::Tau, PATMuonIpSigExtractor, PATTauIpSigExtractor> PATMuTauIpSigExtractor;

typedef ObjectSelector<PATdiTauIpSigSelector<pat::Electron, pat::Muon, PATElecMuIpSigExtractor> > PATElecMuIpSigSelector;
typedef ObjectSelector<PATdiTauIpSigSelector<pat::Electron, pat::Tau, PATElecTauIpSigExtractor> > PATElecTauIpSigSelector;
typedef ObjectSelector<PATdiTauIpSigSelector<pat::Muon, pat::Tau, PATMuTauIpSigExtractor> > PATMuTauIpSigSelector;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_ANOTHER_FWK_MODULE(PATElecMuIpSigSelector);
DEFINE_ANOTHER_FWK_MODULE(PATElecTauIpSigSelector);
DEFINE_ANOTHER_FWK_MODULE(PATMuTauIpSigSelector);
