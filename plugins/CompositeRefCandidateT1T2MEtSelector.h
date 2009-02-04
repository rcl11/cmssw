//
// $Id: CompositeRefCandidateT1T2MEtSelector.h,v 1.3 2008/06/19 13:22:12 gpetrucc Exp $
//

#ifndef TauAnalysis_CandidateTools_CompositeRefCandidateT1T2MEtSelector_h
#define TauAnalysis_CandidateTools_CompositeRefCandidateT1T2MEtSelector_h

#include "DataFormats/Common/interface/RefVector.h"

#include "PhysicsTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "PhysicsTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"
#include "PhysicsTools/UtilAlgos/interface/SingleElementCollectionSelector.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositeRefCandidateT1T2MEt.h"

#include <vector>

typedef SingleObjectSelector<
            std::vector<DiCandidatePair>,
            StringCutObjectSelector<DiCandidatePair>
        > DiCandidatePairSelector;
typedef SingleObjectSelector<
            std::vector<PATElecTauPair>,
            StringCutObjectSelector<PATElecTauPair>
        > PATElecTauPairSelector;
typedef SingleObjectSelector<
            std::vector<PATMuTauPair>,
            StringCutObjectSelector<PATMuTauPair>
        > PATMuTauPairSelector;
typedef SingleObjectSelector<
            std::vector<PATDiTauPair>,
            StringCutObjectSelector<PATDiTauPair>
        > PATDiTauPairSelector;
typedef SingleObjectSelector<
            std::vector<PATElecMuPair>,
            StringCutObjectSelector<PATElecMuPair>
        > PATElecMuPairSelector;

typedef SingleObjectSelector<
            std::vector<DiCandidatePair>,
            StringCutObjectSelector<DiCandidatePair>,
            edm::RefVector<std::vector<DiCandidatePair> >
        > DiCandidatePairRefSelector;
typedef SingleObjectSelector<
            std::vector<PATElecTauPair>,
            StringCutObjectSelector<PATElecTauPair>,
            edm::RefVector<std::vector<PATElecTauPair> >
        > PATElecTauPairRefSelector;
typedef SingleObjectSelector<
            std::vector<PATMuTauPair>,
            StringCutObjectSelector<PATMuTauPair>,
            edm::RefVector<std::vector<PATMuTauPair> >
        > PATMuTauPairRefSelector;
typedef SingleObjectSelector<
            std::vector<PATDiTauPair>,
            StringCutObjectSelector<PATDiTauPair>,
            edm::RefVector<std::vector<PATDiTauPair> >
        > PATDiTauPairRefSelector;
typedef SingleObjectSelector<
            std::vector<PATElecMuPair>,
            StringCutObjectSelector<PATElecMuPair>,
            edm::RefVector<std::vector<PATElecMuPair> >
        > PATElecMuPairRefSelector;

#endif
