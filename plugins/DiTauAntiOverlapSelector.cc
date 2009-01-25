#include "TauAnalysis/CandidateTools/plugins/DiTauAntiOverlapSelector.h"

#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"

#include "FWCore/Framework/interface/MakerMacros.h"

typedef ObjectSelector<DiTauAntiOverlapSelectorImp> DiTauAntiOverlapSelector;
DEFINE_ANOTHER_FWK_MODULE(DiTauAntiOverlapSelector);


