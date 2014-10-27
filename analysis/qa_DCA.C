TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;
const char *run_config = "PrimTrk.ClosePrimVtx";

//================================================
void qa_DCA()
{
  gStyle->SetOptStat(0);

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.%s.root",run_config),"read");
  
  gDCA();
}

//================================================
void gDCA(const Int_t save = 1)
{
  TH2F *hTrkDca = (TH2F*)f->Get(Form("hTrkDca_qa_%s",trigName[kTrigType]));
  c = draw2D(hTrkDca,Form("Au+Au %s: gDCA vs p_{T} of %s tracks%s",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_DCA/%s.%s_gDCA_vs_pt_%s.png",run_config,trk_name[trk_index],trigName[kTrigType]));

  TH1F *hDca = (TH1F*)hTrkDca->ProjectionY(Form("hDca_qa_%s",trigName[kTrigType]));
  c = draw1D(hDca,Form("Au+Au %s: gDCA of %s tracks%s",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]),kTRUE,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_DCA/%s.%s_gDCA_%s.png",run_config,trk_name[trk_index],trigName[kTrigType]));
}


void qa()
{
  TString mGeomTag = "y2012a";
  TString ts = Form("$STAR/StarVMC/Geometry/macros/loadStarGeometry.C(\"%s\",1)",mGeomTag.Data());
  Int_t ierr=0;
  gROOT->Macro(ts.Data(),&ierr);
  assert(!ierr);
  assert(gGeoManager);

  // intialize backleg/module geometry
  TGeoVolume *mMtdGeom = gGeoManager->FindVolumeFast("MUTD");
  const char *elementName = mMtdGeom->GetName();
  if(elementName){
    TGeoIterator next(mMtdGeom);
    next.SetTopName("/HALL_1/CAVE_1/MUTD_1");
    TGeoNode   *node = 0;
    TGeoVolume *blVol = 0;
    TGeoVolume *modVol = 0;
    while ( (node=(TGeoNode*)next()) )
      {
	TString name = node->GetName();
	TString path;
	next.GetPath(path);
	if(!gGeoManager->CheckPath(path.Data())){
	  LOG_WARN<<"Path "<<path.Data()<<" is not found"<<endm;
	  continue;
	}
	gGeoManager->cd(path.Data());
	TGeoVolume *detVol = gGeoManager->GetCurrentVolume();    

	//fill GeoBLs 
	if(IsMTTG(detVol)){
	  blVol = (TGeoVolume *)detVol;
	  Double_t op[3];
	  Double_t local[3] = {0,0,0};
	  gGeoManager->LocalToMaster(local,op);
	  TGeoHMatrix *mat = gGeoManager->GetCurrentMatrix();
	}
	//fill GeoModules
	if(IsMTRA(detVol)){
	  modVol = (TGeoVolume *)detVol;
	  Double_t op[3];
	  Double_t local[3] = {0,0,0};
	  gGeoManager->LocalToMaster(local,op);
	  TGeoHMatrix *mat = gGeoManager->GetCurrentMatrix();
	}
      }
  }
}
