#ifndef TauAnalysis_CandidateTools_NSVfitSingleTauBuilderBase_h
#define TauAnalysis_CandidateTools_NSVfitSingleTauBuilderBase_h

/** \class NSVfitTauDecayBuilderBase
 *
 * Base-class for building objects that come from tau decays.
 *
 * \author Evan K. Friis, Christian Veelken, UC Davis
 *
 * \version $Revision: 1.4 $
 *
 * $Id: NSVfitTauDecayBuilderBase.h,v 1.4 2011/03/29 14:53:26 veelken Exp $
 *
 */

#include "DataFormats/Candidate/interface/Candidate.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitSingleParticleBuilderBase.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitTrackService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitTauDecayHypothesis.h"

class NSVfitSingleParticleHypothesisBase;
class NSVfitAlgorithmBase;

class NSVfitTauDecayBuilderBase : public NSVfitSingleParticleBuilderBase 
{
  public:
    NSVfitTauDecayBuilderBase(const edm::ParameterSet& cfg)
      : NSVfitSingleParticleBuilderBase(cfg),
        algorithm_(0),
        idxFitParameter_nuInvMass_(-1)
    {}
    virtual ~NSVfitTauDecayBuilderBase() {}

    // Setup the parameters of the fit.
    virtual void beginJob(NSVfitAlgorithmBase*);

    // Build the tau decay hypothesis from the fit parameters
    virtual void applyFitParameter(NSVfitSingleParticleHypothesisBase*, double*) const;

    /* Abstract functions overridden by the different decay type builders */
    // Overridden to allocate the specific decay type.
    virtual bool nuSystemIsMassless() const = 0;
    // The decay mode
    virtual int getDecayMode(const reco::Candidate*) const = 0;
    // Get the track(s) associated to a given Candidate
    virtual std::vector<reco::TrackBaseRef> extractTracks(const reco::Candidate*) const = 0;

    virtual void print(std::ostream&) const;

  protected:
    // Initialize data-members common to tau --> e/mu and tau --> had decays
    void initialize(NSVfitTauDecayHypothesis*, const reco::Candidate*) const;

    NSVfitAlgorithmBase* algorithm_;

    edm::Service<NSVfitTrackService> trackService_;

    int idxFitParameter_visEnFracX_;
    int idxFitParameter_phi_lab_;
    int idxFitParameter_nuInvMass_; // used for leptonic decays only.
    int idxFitParameter_deltaR_;

    int idxFitParameter_pvShiftX_;
    int idxFitParameter_pvShiftY_;
    int idxFitParameter_pvShiftZ_;
};

void applyOptionalFitParameter(double*, int, double&);

#endif /* end of include guard: TauAnalysis_CandidateTools_NSVfitSingleTauBuilderBase_h */
