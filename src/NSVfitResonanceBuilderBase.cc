#include "TauAnalysis/CandidateTools/interface/NSVfitResonanceBuilderBase.h"

NSVfitResonanceBuilderBase::NSVfitResonanceBuilderBase(const edm::ParameterSet& cfg)
  : NSVfitBuilderBase(cfg),
    prodResonanceLabel_(cfg.getParameter<std::string>("prodResonanceLabel")),
    numDaughterBuilders_(0)
{
  edm::ParameterSet cfg_daughters = cfg.getParameter<edm::ParameterSet>("daughters");
  typedef std::vector<std::string> vstring;
  vstring daughterNames = cfg_daughters.getParameterNamesForType<edm::ParameterSet>();
  for ( vstring::const_iterator daughterName = daughterNames.begin();
	daughterName != daughterNames.end(); ++daughterName ) {
    edm::ParameterSet cfg_daughter = cfg_daughters.getParameter<edm::ParameterSet>(*daughterName);
    edm::ParameterSet cfg_builder = cfg_daughter.getParameter<edm::ParameterSet>("builder");
    cfg_builder.addParameter<std::string>("prodParticleLabel", *daughterName);
    std::string pluginType = cfg_builder.getParameter<std::string>("pluginType");
    NSVfitSingleParticleBuilderBase* daughterBuilder =
      NSVfitSingleParticleBuilderPluginFactory::get()->create(pluginType, cfg_builder);
    daughterBuilders_.push_back(daughterBuilder);
    ++numDaughterBuilders_;
  }
}

NSVfitResonanceBuilderBase::~NSVfitResonanceBuilderBase()
{
  for ( std::vector<NSVfitSingleParticleBuilderBase*>::iterator it = daughterBuilders_.begin();
	it != daughterBuilders_.end(); ++it ) {
    delete (*it);
  }
}

void NSVfitResonanceBuilderBase::beginJob(NSVfitAlgorithmBase* algorithm)
{
  for ( std::vector<NSVfitSingleParticleBuilderBase*>::iterator daughterBuilder = daughterBuilders_.begin();
	daughterBuilder != daughterBuilders_.end(); ++daughterBuilder ) {
    (*daughterBuilder)->beginJob(algorithm);
  }
}

void NSVfitResonanceBuilderBase::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
  NSVfitBuilderBase::beginEvent(evt, es);
  for ( std::vector<NSVfitSingleParticleBuilderBase*>::iterator daughterBuilder = daughterBuilders_.begin();
	daughterBuilder != daughterBuilders_.end(); ++daughterBuilder ) {
    (*daughterBuilder)->beginEvent(evt, es);
  }
}

NSVfitResonanceHypothesis* NSVfitResonanceBuilderBase::build(
    const inputParticleMap& inputParticles) const
{
  NSVfitResonanceHypothesis* resonanceHypothesis = new NSVfitResonanceHypothesis();

  reco::Candidate::LorentzVector p4(0,0,0,0);

  for ( std::vector<NSVfitSingleParticleBuilderBase*>::const_iterator daughterBuilder = daughterBuilders_.begin();
	daughterBuilder != daughterBuilders_.end(); ++daughterBuilder ) {
    NSVfitSingleParticleHypothesisBase* daughterHypothesis =
      (*daughterBuilder)->build(inputParticles);
    daughterHypothesis->setMother(resonanceHypothesis);

    p4 += daughterHypothesis->p4();

    resonanceHypothesis->daughters_.push_back(daughterHypothesis);
  }

  resonanceHypothesis->p4_ = p4;

  resonanceHypothesis->name_ = prodResonanceLabel_;

  resonanceHypothesis->barcode_ = barcodeCounter_;
  ++barcodeCounter_;

  return resonanceHypothesis;
}

void NSVfitResonanceBuilderBase::applyFitParameter(NSVfitResonanceHypothesis* resonanceHypothesis, double* params) const
{
  std::cout << "<NSVfitResonanceBuilderBase::applyFitParameter>:" << std::endl;
  for ( unsigned iDaughterBuilder = 0; iDaughterBuilder < numDaughterBuilders_; ++iDaughterBuilder ) {
    daughterBuilders_[iDaughterBuilder]->applyFitParameter(&resonanceHypothesis->daughters_[iDaughterBuilder], params);
  }

  reco::Candidate::LorentzVector dp4(0,0,0,0);

  edm::OwnVector<NSVfitSingleParticleHypothesisBase> daughters = resonanceHypothesis->daughters_;
  for ( edm::OwnVector<NSVfitSingleParticleHypothesisBase>::const_iterator daughter = daughters.begin();
	daughter != daughters.end(); ++daughter ) {
    dp4 += daughter->dp4_fitted();
  }

  resonanceHypothesis->dp4_ = dp4;
}

void NSVfitResonanceBuilderBase::print(std::ostream& stream) const
{
  stream << "<NSVfitResonanceBuilderBase::print>:" << std::endl;
  stream << " pluginName = " << pluginName_ << std::endl;
  stream << " pluginType = " << pluginType_ << std::endl;
  stream << " prodResonanceLabel = " << prodResonanceLabel_ << std::endl;
  for ( std::vector<NSVfitSingleParticleBuilderBase*>::const_iterator daughterBuilder = daughterBuilders_.begin();
	daughterBuilder != daughterBuilders_.end(); ++daughterBuilder ) {
    (*daughterBuilder)->print(stream);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

EDM_REGISTER_PLUGINFACTORY(NSVfitResonanceBuilderPluginFactory, "NSVfitResonanceBuilderPluginFactory");


