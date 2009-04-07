#ifndef TauAnalysis_RecoTools_PATdiTauIpSigSelector_h
#define TauAnalysis_RecoTools_PATdiTauIpSigSelector_h

/** \class PATdiTauIpSigSelector
 *
 * Select di-tau piars based on track 
 * transverse impact parameter significance
 * with respect to primary event vertex
 * 
 * \author G. Cerati
 *
 * \version $Revision: 1.1 $
 *
 * $Id: PATdiTauIpSigSelector.h,v 1.2 2009/03/23 14:18:59 cerati Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"

#include <vector>

template <typename T1, typename T2, typename TExtr>
class PATdiTauIpSigSelector
{
  public:
    
    typedef std::vector<CompositePtrCandidateT1T2MEt<T1,T2> > collection;
  
    explicit PATdiTauIpSigSelector(const edm::ParameterSet&);
    ~PATdiTauIpSigSelector();

    typename std::vector<const CompositePtrCandidateT1T2MEt<T1,T2> *>::const_iterator begin() const { return selected_.begin(); }
    typename std::vector<const CompositePtrCandidateT1T2MEt<T1,T2> *>::const_iterator end() const { return selected_.end(); }

    void select(const edm::Handle<collection>&, edm::Event&, const edm::EventSetup&);
    
    size_t size() const { return selected_.size(); }

  private:
    std::vector<const CompositePtrCandidateT1T2MEt<T1,T2> *> selected_;
    edm::InputTag vertexSrc_;
    double ipMin_;
    bool applyIpMin_;
    double ipMax_;
    bool applyIpMax_;

    TExtr ipSigExtractor_;
};

#endif
