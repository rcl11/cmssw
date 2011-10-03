#ifndef TauAnalysis_CandidateTools_interface_NSVfitEventAnalyzer_h
#define TauAnalysis_CandidateTools_interface_NSVfitEventAnalyzer_h

#include <map>
#include <string>

#include "TH1.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

class PFMEtSignInterface;

/**
   \class NSVfitEventAnalyzer NSVfitEventAnalyzer.h "TauAnalysis/CandidateTools/interface/NSVfitEventAnalyzer.h"
   \brief Basic edm and fwlite friendly analyzer class to do basic testing of NSVfit

   This is an example for keeping classes that can be used both within FWLite and within the full 
   framework. The class is derived from the BasicAnalyzer base class, which is an interface for 
   the two wrapper classes EDAnalyzerWrapper and FWLiteAnalyzerWrapper. The latter provides basic 
   configuration file reading and event looping equivalent to the FWLiteHistograms executable of 
   this package. You can see the FWLiteAnalyzerWrapper class at work in the FWLiteWithBasicAnalyzer
   executable of this package.
*/

class NSVfitEventAnalyzer : public edm::EDAnalyzer {

 public:
  /// default constructor
  NSVfitEventAnalyzer(const edm::ParameterSet& cfg);
  /// default destructor
  virtual ~NSVfitEventAnalyzer();
  /// everything that needs to be done before the event loop
  void beginJob(){};
  /// everything that needs to be done after the event loop
  void endJob(){};
  /// everything that needs to be done during the event loop
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

 private:
  /// input tag for MET
  edm::InputTag met_;
  /// input tag for electrons
  edm::InputTag leps1_;
  /// input tag for muons
  edm::InputTag leps2_;

  /// MET significance interface
  PFMEtSignInterface* metSign_;
  /// histograms
  std::map<std::string, TH1*> hists_;
};

#endif
