#ifndef L1Trigger_DTTrigger_DTTrigTest_cc
#define L1Trigger_DTTrigger_DTTrigTest_cc

#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/Event.h>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "L1Trigger/DTTrigger/interface/DTTrig.h"
#include "L1Trigger/DTTriggerServerPhi/interface/DTChambPhSegm.h"
#include "SimDataFormats/Track/interface/EmbdSimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/EmbdSimVertexContainer.h"
#include <CLHEP/Vector/LorentzVector.h>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <math.h>

using namespace std;
using namespace edm;

class DTTrigTest: public EDAnalyzer{
public:
  DTTrigTest(const ParameterSet& pset){
    MyTrig = new DTTrig();
    string outputfile = pset.getUntrackedParameter<string>("outputFileName");
    string rootext = ".root";
    f = new TFile((outputfile+rootext).c_str(),"RECREATE");
    theTree = new TTree("h1","GMT",0);
    //MyTrig->config()->setParam("Debugging level","fullTRACO");
    //cout << MyTrig->config()->trigSetupGeom() << endl;
    cout << "constructor executed!!!" << endl;
  }

  ~DTTrigTest(){ 
    delete MyTrig;
    delete f;
    cout << "destructor executed!!!" << endl;
  }

  void endJob(){
    cout << "Writing Tree and Closing File" << endl;
    theTree->Write();
    delete theTree;
    f->Close();
  }

  void beginJob(const EventSetup & iEventSetup){
    
    MyTrig->createTUs(iEventSetup);
    cout << "****TU's Created" << endl;
    
    // BOOKING of the tree's varables
    // GENERAL block branches
    theTree->Branch("Run",&runn,"Run/I");
    theTree->Branch("Event",&eventn,"Event/I");
    theTree->Branch("Weight",&weight,"Weight/F");  
    // GEANT block branches
    theTree->Branch("Ngen",&ngen,"Ngen/I");
    theTree->Branch("Pxgen",pxgen,"Pxgen[Ngen]/F");
    theTree->Branch("Pygen",pygen,"Pygen[Ngen]/F");
    theTree->Branch("Pzgen",pzgen,"Pzgen[Ngen]/F");
    theTree->Branch("Ptgen",ptgen,"Ptgen[Ngen]/F");
    theTree->Branch("Etagen",etagen,"Etagen[Ngen]/F");
    theTree->Branch("Phigen",phigen,"Phigen[Ngen]/F");
    theTree->Branch("Chagen",chagen,"Chagen[Ngen]/I");
    theTree->Branch("Vxgen",vxgen,"Vxgen[Ngen]/F");
    theTree->Branch("Vygen",vygen,"Vygen[Ngen]/F");
    theTree->Branch("Vzgen",vzgen,"Vzgen[Ngen]/F");
    // L1MuDTBtiChipS block
    theTree->Branch("Nbti",&nbti,"Nbti/I");
    theTree->Branch("bwh",bwh,"bwh[Nbti]/I"); 
    theTree->Branch("bstat",bstat,"bstat[Nbti]/I");    
    theTree->Branch("bsect",bsect,"bsect[Nbti]/I");  
    theTree->Branch("bsl",bsl,"bsl[Nbti]/I");
    theTree->Branch("bnum",bnum,"bnum[Nbti]/I");
    theTree->Branch("bbx",bbx,"bbx[Nbti]/I");
    theTree->Branch("bcod",bcod,"bcod[Nbti]/I");
    theTree->Branch("bk",bk,"bk[Nbti]/I");
    theTree->Branch("bx",bx,"bx[Nbti]/I");
    theTree->Branch("bposx",bposx,"bposx[Nbti]/F");
    theTree->Branch("bposy",bposy,"bposy[Nbti]/F");
    theTree->Branch("bposz",bposz,"bposz[Nbti]/F");
    theTree->Branch("bdirx",bdirx,"bdirx[Nbti]/F");
    theTree->Branch("bdiry",bdiry,"bdiry[Nbti]/F");
    theTree->Branch("bdirz",bdirz,"bdirz[Nbti]/F");
    // L1MuDTTracoChipS block
    theTree->Branch("Ntraco",&ntraco,"Ntraco/I");
    theTree->Branch("twh",twh,"twh[Ntraco]/I"); 
    theTree->Branch("tstat",tstat,"tstat[Ntraco]/I");    
    theTree->Branch("tsect",tsect,"tsect[Ntraco]/I");  
    theTree->Branch("tnum",tnum,"tnum[Ntraco]/I"); 
    theTree->Branch("tbx",tbx,"tbx[Ntraco]/I");
    theTree->Branch("tcod",tcod,"tcod[Ntraco]/I");
    theTree->Branch("tk",tk,"tk[Ntraco]/I");
    theTree->Branch("tx",tx,"tx[Ntraco]/I");
    theTree->Branch("tposx",tposx,"tposx[Ntraco]/F");
    theTree->Branch("tposy",tposy,"tposy[Ntraco]/F");
    theTree->Branch("tposz",tposz,"tposz[Ntraco]/F");
    theTree->Branch("tdirx",tdirx,"tdirx[Ntraco]/F");
    theTree->Branch("tdiry",tdiry,"tdiry[Ntraco]/F");
    theTree->Branch("tdirz",tdirz,"tdirz[Ntraco]/F");
    // TSPHI block
    theTree->Branch("Ntsphi",&ntsphi,"Ntsphi/I");
    theTree->Branch("swh",swh,"swh[Ntsphi]/I"); 
    theTree->Branch("sstat",sstat,"sstat[Ntsphi]/I");    
    theTree->Branch("ssect",ssect,"ssect[Ntsphi]/I");  
    theTree->Branch("sbx",sbx,"sbx[Ntsphi]/I");
    theTree->Branch("scod",scod,"scod[Ntsphi]/I");
    theTree->Branch("sphi",sphi,"sphi[Ntsphi]/I");
    theTree->Branch("sphib",sphib,"sphib[Ntsphi]/I");
    theTree->Branch("sposx",sposx,"sposx[Ntsphi]/F");
    theTree->Branch("sposy",sposy,"sposy[Ntsphi]/F");
    theTree->Branch("sposz",sposz,"sposz[Ntsphi]/F");
    theTree->Branch("sdirx",sdirx,"sdirx[Ntsphi]/F");
    theTree->Branch("sdiry",sdiry,"sdiry[Ntsphi]/F");
    theTree->Branch("sdirz",sdirz,"sdirz[Ntsphi]/F");
    // TSTHETA block
    theTree->Branch("Ntstheta",&ntstheta,"Ntstheta/I");
    theTree->Branch("thwh",thwh,"thwh[Ntstheta]/I"); 
    theTree->Branch("thstat",thstat,"thstat[Ntstheta]/I");    
    theTree->Branch("thsect",thsect,"thsect[Ntstheta]/I");  
    theTree->Branch("thbx",thbx,"thbx[Ntstheta]/I");
    theTree->Branch("thcode",thcode,"thcode[Ntstheta][7]/I");
    theTree->Branch("thpos",thpos,"thpos[Ntstheta][7]/I");
    theTree->Branch("thqual",thqual,"thqual[Ntstheta][7]/I");
  }
  
  void analyze(const Event & iEvent, const EventSetup& iEventSetup){

    const int MAXGEN  = 10;
    const float ptcut  = 1.0;
    const float etacut = 2.4;

    MyTrig->triggerReco(iEvent,iEventSetup);
    cout << "****Trigger algorithm executed for run " << iEvent.id().run() <<" event " << iEvent.id().event() << endl;
    
    // GENERAL Block
    runn   = iEvent.id().run();
    eventn = iEvent.id().event();
    weight = 1; // FIXME what to do with this varable?
    
    // GEANT Block
    Handle<vector<EmbdSimTrack> > MyTracks;
    Handle<vector<EmbdSimVertex> > MyVertexes;
    iEvent.getByLabel("SimG4Object",MyTracks);
    iEvent.getByLabel("SimG4Object",MyVertexes);
    vector<EmbdSimTrack>::const_iterator itrack;
    ngen=0;
    cout  << "Tracks found in the detector" << MyTracks->size() <<endl;
    for (itrack=MyTracks->begin(); itrack!=MyTracks->end(); itrack++){
      if ( abs(itrack->type())==13){
	float pt  = itrack->momentum().perp();
	float eta = itrack->momentum().pseudoRapidity();
	if ( pt>ptcut && fabs(eta)<etacut ){
	  HepLorentzVector momentum = itrack->momentum();
	  float phi = momentum.phi();
	  int charge = -1; // static_cast<int> (itrack->charge()); charge() still to be implemented
	  if ( phi<0 ) phi = 2*M_PI + phi;
	  int vtxindex = itrack->vertIndex();
	    float gvx=0,gvy=0,gvz=0;
	  if (vtxindex >-1){
	    gvx=MyVertexes->at(vtxindex).position().x();
	    gvy=MyVertexes->at(vtxindex).position().y();
	    gvz=MyVertexes->at(vtxindex).position().z();
	  }
	  if ( ngen < MAXGEN ) {
	    pxgen[ngen]=momentum.x();
	    pygen[ngen]=momentum.y();
	    pzgen[ngen]=momentum.z();
	    ptgen[ngen]=pt;
	    etagen[ngen]=eta;
	    phigen[ngen]=phi;
	    chagen[ngen]=charge;
	    vxgen[ngen]=gvx;
	    vygen[ngen]=gvy;
	    vzgen[ngen]=gvz;
   	    ngen++;
	  }
	}
      }
    }
    
    // L1 Local Trigger Block
    // BTI
    vector<DTBtiTrigData> btitrigs = MyTrig->BtiTrigs();
    vector<DTBtiTrigData>::const_iterator pbti;
    int ibti = 0;
    cout << btitrigs.size() << " BTI triggers found" << endl;
    for ( pbti = btitrigs.begin(); pbti != btitrigs.end(); pbti++ ) {
      if ( ibti < 100 ) {
	bwh[ibti]=pbti->wheel();
	//cout << "BTI Wheel :" << pbti->wheel() <<"    " << bwh[ibti]<< endl;
	bstat[ibti]=pbti->station();
	bsect[ibti]=pbti->sector();
	bsl[ibti]=pbti->btiSL();
	bnum[ibti]=pbti->btiNumber();
	bbx[ibti]=pbti->step();
	bcod[ibti]=pbti->code();
	bk[ibti]=pbti->K();
	bx[ibti]=pbti->X();
	GlobalPoint pos = MyTrig->CMSPosition(&(*pbti));
	GlobalVector dir = MyTrig->CMSDirection(&(*pbti));
	bposx[ibti] = pos.x();
	bposy[ibti] = pos.y();
	bposz[ibti] = pos.z();
	bdirx[ibti] = dir.x();
	bdiry[ibti] = dir.y();
	bdirz[ibti] = dir.z();
	ibti++;
      }
    } 
    nbti = ibti;
    //cout << nbti << endl;
    
    //TRACO
    vector<DTTracoTrigData> tracotrigs = MyTrig->TracoTrigs();
    vector<DTTracoTrigData>::const_iterator ptc;
    int itraco = 0;
    cout << tracotrigs.size() << " TRACO triggers found" << endl;
    for (ptc=tracotrigs.begin(); ptc!=tracotrigs.end(); ptc++) {
      if (itraco<80) {
	twh[itraco]=ptc->wheel();
	tstat[itraco]=ptc->station();
	tsect[itraco]=ptc->sector();
	tnum[itraco]=ptc->tracoNumber();
	tbx[itraco]=ptc->step();
	tcod[itraco]=ptc->code();
	tk[itraco]=ptc->K();
	tx[itraco]=ptc->X();
	GlobalPoint pos = MyTrig->CMSPosition(&(*ptc));
	GlobalVector dir = MyTrig->CMSDirection(&(*ptc));
	tposx[itraco] = pos.x();
	tposy[itraco] = pos.y();
	tposz[itraco] = pos.z();
	tdirx[itraco] = dir.x();
	tdiry[itraco] = dir.y();
	tdirz[itraco] = dir.z();
	itraco++;
      }
    }
    ntraco = itraco;
    //cout << ntraco<<endl;
    
    //TSPHI
    vector<DTChambPhSegm> tsphtrigs = MyTrig->TSPhTrigs();
    vector<DTChambPhSegm>::const_iterator ptsph;
    int itsphi = 0; 
    cout << tsphtrigs.size() << " TSPhi triggers found" << endl;
    for (ptsph=tsphtrigs.begin(); ptsph!=tsphtrigs.end(); ptsph++) {
      if (itsphi<40 ) {
	const DTChambPhSegm& seg = (*ptsph);
	swh[itsphi] = ptsph->wheel();
	sstat[itsphi] = ptsph->station();
	ssect[itsphi] = ptsph->sector();
	sbx[itsphi] = ptsph->step();      
	scod[itsphi] = ptsph->oldCode();
	sphi[itsphi] = ptsph->phi();
	sphib[itsphi] = ptsph->phiB();
	GlobalPoint pos = MyTrig->CMSPosition(&seg); 
	GlobalVector dir = MyTrig->CMSDirection(&seg);
	sposx[itsphi] = pos.x();
	sposy[itsphi] = pos.y();
	sposz[itsphi] = pos.z();
	sdirx[itsphi] = dir.x();
	sdiry[itsphi] = dir.y();
	sdirz[itsphi] = dir.z();
	itsphi++;
      }
    }
    ntsphi = itsphi;
    //cout << ntsphi << endl;

//TSPHI
    vector<DTChambThSegm> tsthtrigs = MyTrig->TSThTrigs();
    vector<DTChambThSegm>::const_iterator ptsth;
    int itstheta = 0; 
    cout << tsthtrigs.size() << " TSTheta triggers found" << endl;
    for (ptsth=tsthtrigs.begin(); ptsth!=tsthtrigs.end(); ptsth++) {
      if (itstheta<40 ) {
	thwh[itstheta] = ptsth->ChamberId().wheel();
	thstat[itstheta] = ptsth->ChamberId().station();
	thsect[itstheta] = ptsth->ChamberId().sector();
	thbx[itstheta] = ptsth->step();
	for(int i=0;i<7;i++) {
	  thcode[itstheta][i] = ptsth->code(i);
	  thpos[itstheta][i] = ptsth->position(i);
	  thqual[itstheta][i] = ptsth->quality(i);
	}
	itstheta++;
      }
    }
    ntstheta = itstheta;
    
    //Fill the tree
    theTree->Fill();
}
  
private:

  //trigger istance
  DTTrig* MyTrig;

  // tree
  TTree* theTree;
  // TFile
  TFile *f;

  //GENERAL block
  int             runn;
  int             eventn;
  float           weight;

  //GEANT block
  int             ngen;
  float           pxgen[10];
  float           pygen[10];
  float           pzgen[10];
  float           ptgen[10];
  float           etagen[10];
  float           phigen[10];
  int             chagen[10];
  float           vxgen[10];
  float           vygen[10];
  float           vzgen[10];
  
  // BTI
  int nbti;
  int bwh[100];
  int bstat[100];
  int bsect[100];
  int bsl[100];
  int bnum[100];
  int bbx[100];
  int bcod[100];
  int bk[100];
  int bx[100];
  float bposx[100];
  float bposy[100];
  float bposz[100];
  float bdirx[100];
  float bdiry[100];
  float bdirz[100];
  
  // TRACO
  int ntraco;
  int twh[80];
  int tstat[80];
  int tsect[80];
  int tnum[80];
  int tbx[80];
  int tcod[80];
  int tk[80];
  int tx[80];
  float tposx[100];
  float tposy[100];
  float tposz[100];
  float tdirx[100];
  float tdiry[100];
  float tdirz[100];
  
  // TSPHI
  int ntsphi;
  int swh[40];
  int sstat[40]; 
  int ssect[40];
  int sbx[40];
  int scod[40];
  int sphi[40];
  int sphib[40];
  float sposx[100];
  float sposy[100];
  float sposz[100];
  float sdirx[100];
  float sdiry[100];
  float sdirz[100]; 

  // TSTHETA
  int ntstheta;
  int thwh[40];
  int thstat[40]; 
  int thsect[40];
  int thbx[40];
  int thcode[40][7];
  int thpos[40][7];
  int thqual[40][7];

};
 
#endif

