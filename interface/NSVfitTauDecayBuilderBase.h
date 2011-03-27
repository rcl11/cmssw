#ifndef TauAnalysis_CandidateTools_NSVfitSingleTauBuilderBase_h
#define TauAnalysis_CandidateTools_NSVfitSingleTauBuilderBase_h

/** \class NSVfitTauDecayBuilderBase
 *
 * Base-class for building objects that come from tau decays.
 *
 * \author Evan K. Friis, Christian Veelken, UC Davis
 *
 * \version $Revision:$
 *
 * $Id:$
 *
 */

#include "TauAnalysis/CandidateTools/interface/NSVfitSingleParticleBuilderBase.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitTauDecayHypothesis.h"

class NSVfitSingleParticleHypothesisBase;
class NSVfitAlgorithmBase;

class NSVfitTauDecayBuilderBase : public NSVfitSingleParticleBuilderBase {
  public:
    NSVfitTauDecayBuilderBase(const edm::ParameterSet& cfg):
      NSVfitSingleParticleBuilderBase(cfg),algorithm_(NULL){}
    virtual ~NSVfitTauDecayBuilderBase() {}
    // Build the basic single particle hypothesis corresponding to this tau
    NSVfitSingleParticleHypothesisBase* build(
        const inputParticleMap& inputParticles) const;

    // Build the tau decay hypothesis from the fit parameters
    void applyFitParameter(NSVfitSingleParticleHypothesisBase*, double*) const;
    // Setup the parameters of the fit.
    void beginJob(NSVfitAlgorithmBase*);

    /* Abstract functions overridden by the different decay type builders */
    // Overridden to allocate the specific decay type.
    virtual NSVfitTauDecayHypothesis* buildSpecific(
        const edm::Ptr<reco::Candidate> particle,
        const std::string& label, int barcode) const = 0;
    virtual bool nuSystemIsMassless() const=0;
    // The decay mode
    virtual int getDecayMode(const reco::Candidate*) const = 0;
    // Get the track(s) associated to a given Candidate
    virtual std::vector<reco::TrackBaseRef> extractTracks(
        const reco::Candidate*) const=0;
    virtual void beginJobSpecific(NSVfitAlgorithmBase* algorithm) = 0;
    // Apply any extra fit parameters.
    virtual void applyFitParameterSpecific(
        NSVfitTauDecayHypothesis* hypotheis, double *param) const = 0;

    virtual void print(std::ostream&) const;

  protected:
    NSVfitAlgorithmBase* algorithm_;
  private:
    int idxFitParameter_visEnFracX_;
    int idxFitParameter_phi_lab_;
    int idxFitParameter_nuInvMass_; // only used for leptonic decays.
    int idxFitParameter_deltaR_;
};
#endif /* end of include guard: TauAnalysis_CandidateTools_NSVfitSingleTauBuilderBase_h */
