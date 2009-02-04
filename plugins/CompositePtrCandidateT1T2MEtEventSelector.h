//
// $Id: CompositePtrCandidateT1T2MEtEventSelector.h,v 1.1 2009/02/04 15:41:35 veelken Exp $
//

#ifndef TauAnalysis_CandidateTools_CompositePtrCandidateT1T2MEtEventSelector_h
#define TauAnalysis_CandidateTools_CompositePtrCandidateT1T2MEtEventSelector_h

#include "PhysicsTools/UtilAlgos/interface/AnySelector.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectCountEventSelector.h"
#include "PhysicsTools/UtilAlgos/interface/MinNumberSelector.h"
#include "PhysicsTools/PatUtils/interface/MaxNumberSelector.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"

#include <vector>

typedef ObjectCountEventSelector<std::vector<DiCandidatePair>, AnySelector, MinNumberSelector> DiCandidatePairMinEventSelector;
typedef ObjectCountEventSelector<std::vector<PATElecTauPair>, AnySelector, MinNumberSelector> PATElecTauPairMinEventSelector;
typedef ObjectCountEventSelector<std::vector<PATMuTauPair>, AnySelector, MinNumberSelector> PATMuTauPairMinEventSelector;
typedef ObjectCountEventSelector<std::vector<PATDiTauPair>, AnySelector, MinNumberSelector> PATDiTauPairMinEventSelector;
typedef ObjectCountEventSelector<std::vector<PATElecMuPair>, AnySelector, MinNumberSelector> PATElecMuPairMinEventSelector;

typedef ObjectCountEventSelector<std::vector<DiCandidatePair>, AnySelector, MaxNumberSelector> DiCandidatePairMaxEventSelector;
typedef ObjectCountEventSelector<std::vector<PATElecTauPair>, AnySelector, MaxNumberSelector> PATElecTauPairMaxEventSelector;
typedef ObjectCountEventSelector<std::vector<PATMuTauPair>, AnySelector, MaxNumberSelector> PATMuTauPairMaxEventSelector;
typedef ObjectCountEventSelector<std::vector<PATDiTauPair>, AnySelector, MaxNumberSelector> PATDiTauPairMaxEventSelector;
typedef ObjectCountEventSelector<std::vector<PATElecMuPair>, AnySelector, MaxNumberSelector> PATElecMuPairMaxEventSelector;

#endif
