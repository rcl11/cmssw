
// system include files
#include <memory>

#include "TauAnalysis/CandidateTools/plugins/DiTauProducer.h"

#include "TauAnalysis/CandidateTools/interface/FetchCollection.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

using namespace std;
using namespace edm;
using namespace reco;


DiTauProducer::DiTauProducer(const edm::ParameterSet& iConfig)
  : 
  inputTagHadronicTaus_ ( iConfig.getParameter<InputTag>("hadronicTaus") ),
  inputTagLeptonicTaus_ ( iConfig.getParameter<InputTag>("leptonicTaus") ),
  inputTagMETs_ ( iConfig.getParameter<InputTag>("METs") ),
  verbose_ ( iConfig.getUntrackedParameter<bool>("verbose", true) ),
  useLeadingTaus_ ( iConfig.getParameter<bool>("useLeadingTaus") ),
  metMode_( iConfig.getParameter<int>("metMode"))
{ 

  switch (metMode_) {
  case COLLINEAR_APPROX: 
  case NO_MET:
  case TRANSVERSE_RECO:
    break; // fine
  default:
    // throw exception
    assert(0);
  }

  produces<DiTauCollection >("");
}


DiTauProducer::~DiTauProducer() {}



void
DiTauProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {


  typedef CandidateView::const_iterator IC;

  if(verbose_) {
    cout<<"DiTauProducer "<<endl;
    cout<<"\t"<<inputTagHadronicTaus_<<endl;
    cout<<"\t"<<inputTagLeptonicTaus_<<endl;
    cout<<"\t"<<inputTagMETs_<<endl;
    cout<<"\tuse leading taus: "<<useLeadingTaus_<<endl;
    cout<<"\tMET mode        : "<<metMode_<<endl;
  }

  Handle<CandidateView> hadTaus;
  pf::fetchCollection(hadTaus, 
		      inputTagHadronicTaus_, 
		      iEvent );

  Handle<CandidateView> lepTaus;
  pf::fetchCollection(lepTaus, 
		      inputTagLeptonicTaus_, 
		      iEvent );


  assert( ! hadTaus->empty() );
  assert( ! lepTaus->empty() );

  // in case the input tag for MET is empty, 
  // we can still compute the ditau:
  // - just sum up the 2 legs: Z->ee and mumu
  // add a flag controlling the usage of MET:
  // 0: Collinear approx, 1:no MET, 2: transverse


  Handle<CandidateView> METs;

  CandidatePtr  metPtr;
  if( metMode_ != NO_MET ) {
    // the only case where MET is not necessary 

    if(inputTagMETs_ == InputTag() ) {
      // MET not specified, and asking for the MET to be used!
      // need to replace this by an exception
      assert(0);
    }
  
    pf::fetchCollection(METs, 
			inputTagMETs_, 
			iEvent );

    assert( METs->size() == 1 );

    metPtr = METs->ptrAt(0);
  }

  // here, check whether the same collection is put on 
  // both legs.
  bool sameCollection = false;
   
  if( hadTaus.id() == lepTaus.id() )
    sameCollection = true;

  if(verbose_) {

    if( METs.isValid() ) 
      cout<<"\tMET size = "<<METs->size()<<", ";
    else
      cout<<"\tno MET specified, ";
    cout<<"HTAU size = "<<hadTaus->size()<<", "
	<<"LTAU size = "<<lepTaus->size()<<", ";
    if(sameCollection) 
      cout<<"same collection"<<endl;
    else 
      cout<<"different collections"<<endl;
  }


  unsigned int hTauIdx = 0;
  unsigned int lTauIdx = 0;


  CandidatePtr lTauPtr;
  CandidatePtr hTauPtr;

  std::auto_ptr<DiTauCollection> pOutDiTau(new DiTauCollection() );

  const CandidateView& myHtaus = *hadTaus;
  const CandidateView& myLtaus = *lepTaus;

  if(useLeadingTaus_) {
    // take the leading leptonic tau and the leading hadronic tau.
    // only one DiTau per event

    // this mode should be the default. 

    if( sameCollection ) {
      // cannot put the same collection on both legs and ask to 
      // take the leading particle on both legs!
      // replace the assert by an exception.
      assert(0);
    }

    int idxTmp = 0; double pTtmp = 0;
    for(CVIT iT = myLtaus.begin(); iT != myLtaus.end(); ++iT, idxTmp++ ) {
      if( (*iT).pt() > pTtmp){
	pTtmp = (*iT).pt();
	lTauIdx = idxTmp;
	lTauPtr = myLtaus.ptrAt(idxTmp);
      }
    }
    idxTmp = 0;  pTtmp = 0;
    for(CVIT iT = myHtaus.begin(); iT != myHtaus.end(); ++iT, idxTmp++ ) {
      if( (*iT).pt() > pTtmp){
	pTtmp = (*iT).pt();
	hTauIdx = idxTmp;
	hTauPtr = myHtaus.ptrAt(idxTmp);
      }
    }

    fillDiTau( pOutDiTau, lTauPtr, hTauPtr, metPtr );
  }
  else {
    // make a DiTau resonance for each combination of leptonic and 
    // hadronic tau.

    // here, we get all possible combinations.

    // in the case where the same object (same 3-momentum...) is 
    // found on both legs, the combination is discarded. 
    
    // when the same collection is put on both legs, 
    // we get twice the same ditau:
    // 1-1 discarded
    // 1-2 ok
    // 2-1 ok, = 1-2
    // 2-2 discarded 
    // is this the intended behaviour? 
    // we could check whether a given combination has already been found. 
    
    // possible to check whether the same collection is on both legs 
    // using the productID of the collection. 
    // in this case: 
    // if i1 = i2 : discard
    // if i1>i2 discard
    //    this is equivalent to discard i1 >= i2


    unsigned lIdx=0;
    for(CVIT iL = myLtaus.begin(); 
	iL != myLtaus.end(); ++iL, ++lIdx ) {
     
      unsigned hIdx = 0;
      for(CVIT iH = myHtaus.begin(); iH != myHtaus.end(); ++iH, ++hIdx ) {

	// need to check that the same collection is on both legs
	// in case it is, do the following:
	if( sameCollection && lIdx>=hIdx) continue; 

	lTauPtr = myLtaus.ptrAt(lIdx);
	hTauPtr = myHtaus.ptrAt(hIdx);
	if(verbose_) 
	  cout<<"\t* combination: "<<lIdx<<" "<<hIdx<<endl;
	fillDiTau( pOutDiTau, lTauPtr, hTauPtr, metPtr );
      }
    }
  }

  iEvent.put(pOutDiTau, "");

}


void DiTauProducer::fillDiTau(std::auto_ptr<DiTauCollection>& diTauCollection,
			      const CandidatePtr& lTauIt, 
			      const CandidatePtr& hTauIt, 
			      const CandidatePtr& metIt) const {
  
  
  math::XYZTLorentzVector hP4;
  double hadrA = -1;
  double leptA = -1;
  
  const Candidate& lTau = *lTauIt;
  const Candidate& hTau = *hTauIt;

  checkCombination( lTau, hTau );

  bool valid = false;

  switch (metMode_) {
  case COLLINEAR_APPROX: 
    {
      if( metIt.isNull() ) 
	assert(0);    
      valid = computeDiTau4P(hTau.p4(),
			     lTau.p4(),
			     metIt->p4(),
			     hadrA,leptA,
			     hP4);
      break;
    }
  case NO_MET:
    // just sum the 4 mom of the 2 taus 
    hP4 = lTau.p4();
    hP4 += hTau.p4();
    break; 
  case TRANSVERSE_RECO:
    // just sum the transverse component of the 4 mom of the 2 taus and MET
    // should it be more sophisticated??
    // mass of this diTau is 0, should be m or mT ??
    if( metIt.isNull() ) 
      assert(0);    
    hP4 = lTau.p4();
    hP4 += hTau.p4();
//     hP4 = math::XYZTLorentzVector( hP4.px()+metIt->px(),
//  				   hP4.py()+metIt->py(),
//  				   0,
// 				   hP4.pt()+metIt->pt()
//  				   );
    // more sophisticated solution:
    // P=pT_vis+MET, E=sqrt(mT_vis)+sqrt(mT_invis)
    // jacobian pick of the mass distribution should be eqal m(diTau)
    // is it correct??
    hP4 = math::XYZTLorentzVector( hP4.px()+metIt->px(),
				   hP4.py()+metIt->py(),
				   0,
				   sqrt( hP4.pt()*hP4.pt()+
					 hP4.mass()*hP4.mass() ) +
				   sqrt( metIt->pt()*metIt->pt()+
					 hP4.mass()*hP4.mass() )
				   );
    break; 
  default:
    // throw exception
    assert(0);
  }

  int charge = lTau.charge()+hTau.charge();
  
  DiTauCandidate diTau = DiTauCandidate( charge, hP4, math::XYZPoint(0,0,0), 
					 hadrA, leptA, 
					 cosPhiT( lTau.p4(), hTau.p4() ),
					 valid );  
  diTau.setHadrTau( hTauIt );
  diTau.setLeptTau( lTauIt );
  diTau.setMet( metIt );

  if(verbose_) {
    cout<<"\tProducing DiTau: "<<endl;
    cout<<diTau.print()<<endl;
  }

  diTauCollection->push_back( diTau );
}



//compute diTau mass with collinear approximation
bool
DiTauProducer::computeDiTau4P( const math::XYZTLorentzVector& tau1,  
			       const math::XYZTLorentzVector& tau2,
			       const math::XYZTLorentzVector& met, 
			       double& a1, double &a2,
			       math::XYZTLorentzVector& diTau ) const {
  
  //a1(a2) - fraction of momentum carried by visible products of tau  
  a1 = ( tau1.x()*tau2.y() - tau1.y()*tau2.x() ) /
    ( ( tau1.x()*tau2.y() - tau1.y()*tau2.x() ) 
      + met.x()*tau2.y() -  met.y()*tau2.x() );
  a2 = ( tau1.x()*tau2.y() - tau1.y()*tau2.x() ) /
    ( ( tau1.x()*tau2.y() - tau1.y()*tau2.x() ) 
      - met.x()*tau1.y() + met.y()*tau1.x() );
  
  if(verbose_) {
    cout<<"\t\ta1="<<a1<<", a2="<<a2<<endl;
  }
    
  if( a1 > 1 || a1 <= 0 || a2 > 1 || a2 <= 0 ) {
    diTau = math::XYZTLorentzVector(0,0,0,0);
    if(verbose_)
      cout<<"\t\ta1 or a2 out of range "<<endl;
    return false;
  }
  else {
    diTau = tau1/a1 + tau2/a2;
    return true;
  }
}


bool DiTauProducer::checkCombination(const Candidate& tau1, 
				     const Candidate& tau2) const {

  // check if the 2 taus are in fact the same particle, 
  // by looking at the difference between their 3-momentum
  if( fabs( ( tau1.p4().Vect() - tau2.p4().Vect() ).Mag2() )<1e-9 ) {
    // use the message logger to send an error or a warning?
    cout<<"\t\tDiTauProducer::checkCombination "<<endl;
    cout<<"\t\tsame particles on both legs? "<<endl;
    
    return false;
  }
  else return true; 

}

// compute cosine of angle phi between tau1 and tau2 (in transversal plane)
double DiTauProducer::cosPhiT( const math::XYZTLorentzVector& tau1,  
			       const math::XYZTLorentzVector& tau2 ) const {
  double cos = ( tau1.px()*tau2.px() + tau1.py()*tau2.py() ) /
    ( tau1.pt()*tau2.pt() );
  if(verbose_) {
    cout<<"\tcos(phi)= "<<cos<<endl;
  }
  return cos;
}

void 
DiTauProducer::beginJob(const edm::EventSetup&) {
}

void 
DiTauProducer::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(DiTauProducer);
