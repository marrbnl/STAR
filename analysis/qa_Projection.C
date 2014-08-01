TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;

const char *run_config = "NoEloss.Global.HLT";

//================================================
void qa_Projection()
{						
  gStyle->SetOptStat(0);

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/output/jpsi.AuAu200.Run14.Match.%s.root",run_config),"read");

  //trackYZ();
  //EnergyLoss();
  trackProjection();
}

//================================================
void trackYZ(const Int_t save = 0)
{
  TList *list = new TList;

  THnSparseF *hTrkProjYZ = (THnSparseF*)f->Get(Form("hTrkProjYZ_qa_%s",trigName[kTrigType]));

  const char *title[2] = {"global z","local y"};
  const char *name[2] = {"z","y"};
  TH2F *hBL[2];
  TH1F *hAll[2];

  hTrkProjYZ->GetAxis(5)->SetRange(2,2);
  for(Int_t i=0; i<2; i++)
    {
      hBL[i] = (TH2F*)hTrkProjYZ->Projection(2-i,3);
      hBL[i]->SetName(Form("Track_%s_vs_BL",name[i]));
      if(i==1) hBL[i]->GetYaxis()->SetRangeUser(-50,50);

      c = draw2D(hBL[i],Form("Au+Au %s: %s distribution of projected tracks at MTD per module;backleg*5+module",trigName[kTrigType],title[i]));
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Projection/%s.track_%s_vs_Module_%s.png",run_config,name[i],trigName[kTrigType]));

      hAll[i] = (TH1F*)hBL[i]->ProjectionY(Form("Track_%s",name[i]));    

      c = draw1D(hAll[i],Form("Au+Au %s: %s distribution of projected tracks at MTD",trigName[kTrigType],title[i]),kFALSE,kFALSE);
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Projection/%s.track_%s_%s.png",run_config,name[i],trigName[kTrigType]));
    }

  TH2F *hTrkYZ = (TH2F*)hTrkProjYZ->Projection(1,2);
  hTrkYZ->SetName(Form("Track_y_vs_z"));
  hTrkYZ->GetYaxis()->SetRangeUser(-50,50);
  c = draw2D(hTrkYZ,Form("Au+Au %s: local y vs global z of projected tracks at MTD;global z (cm);local y (cm)",trigName[kTrigType]),0.04,kFALSE);
  gPad->SetRightMargin(0.13);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Projection/%s.track_y_vs_z_%s.png",run_config,trigName[kTrigType]));
  hTrkProjYZ->GetAxis(5)->SetRange(0,-1);

  // non-overlapping
  hTrkProjYZ->GetAxis(5)->SetRange(1,1);
  TH2F *hTrkYvsZ = (TH2F*)hTrkProjYZ->Projection(1,2);
  hTrkYvsZ->SetName(Form("Track_y_vs_z_nonOverlapping"));
  hTrkYvsZ->GetYaxis()->SetRangeUser(-50,50);
  c = draw2D(hTrkYvsZ,Form("Au+Au %s:local y vs global z of projected tracks (non-overlapping);global z (cm);local y (cm)",trigName[kTrigType]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Projection/%s.track_y_vs_z_nonOverlapping_%s.png",run_config,trigName[kTrigType]));

  TH1F *hTrkY = (TH1F*)hTrkYvsZ->ProjectionY(Form("Track_y_nonOverlapping"));    
  c = draw1D(hTrkY,Form("Au+Au %s: %s of projected tracks at MTD (non-overlapping)",trigName[kTrigType],title[1]),kFALSE,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Projection/%s.track_%s_nonOverlapping_%s.png",run_config,name[1],trigName[kTrigType]));
  hTrkProjYZ->GetAxis(5)->SetRange(0,-1);

  // pt dependence
  hTrkProjYZ->GetAxis(5)->SetRange(2,2);
  TH2F *hTrkYvsPt = (TH2F*)hTrkProjYZ->Projection(1,0);
  hTrkYvsPt->SetName(Form("Track_y_vs_pt"));
  hTrkYvsPt->GetYaxis()->SetRangeUser(-50,50);
  c = draw2D(hTrkYvsPt,Form("Au+Au %s:local y vs p_{T} of projected tracks at MTD radius;p_{T} (GeV/c);local y (cm)",trigName[kTrigType]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Projection/%s.track_y_vs_pt_%s.png",run_config,trigName[kTrigType]));
  Double_t pt_cut[5] = {1,1.5,2,5,20};
  TH1F *hTrkYinPt[4];
  list->Clear();
  TString legName_pt[4];
  for(Int_t i=0; i<4; i++)
    {
      Int_t low_bin = hTrkYvsPt->GetXaxis()->FindFixBin(pt_cut[i]+0.1);
      Int_t up_bin  = hTrkYvsPt->GetXaxis()->FindFixBin(pt_cut[i+1]-0.1);
      hTrkYinPt[i] = (TH1F*)hTrkYvsPt->ProjectionY(Form("hTrkY_pt_%1.1f_to_%1.1f",pt_cut[i],pt_cut[i+1]),low_bin,up_bin);
      //hTrkYinPt[i]->Scale(1./hTrkYinPt[i]->Integral());
      list->Add(hTrkYinPt[i]);
      legName_pt[i] = Form("%1.1f < p_{T} < %1.1f GeV/c",pt_cut[i],pt_cut[i+1]);
    }
  c = drawHistos(list,"hTrkY_in_pt",Form("Au+Au %s: local y of projected tracks at MTD radius;y (cm);counts",trigName[kTrigType]),kFALSE,0,100,kTRUE,1e4,1e8,kTRUE,kTRUE,legName_pt,kTRUE,"",0.35,0.55,0.65,0.85,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Projection/%s.track_y_InPtBin_%s.png",run_config,trigName[kTrigType]));
  hTrkProjYZ->GetAxis(5)->SetRange(0,-1);
  list->Clear();

  // charge dependence
  hTrkProjYZ->GetAxis(5)->SetRange(2,2);
  TH2F *hTrkYvsCharge = (TH2F*)hTrkProjYZ->Projection(1,4);
  hTrkYvsCharge->SetName(Form("Track_y_vs_charge"));
  hTrkYvsCharge->GetYaxis()->SetRangeUser(-50,50);
  c = draw2D(hTrkYvsCharge,Form("Au+Au %s:local y vs charge of projected tracks at MTD radius;charge;local y (cm)",trigName[kTrigType]),0.04,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Projection/%s.track_y_vs_charge_%s.png",run_config,trigName[kTrigType]));

  TH1F *hTrkYinCharge[2];
  TString legName_charge[2] = {"negative","positive"};
  for(Int_t i=0; i<2; i++)
    {
      hTrkYinCharge[i] = (TH1F*)hTrkYvsCharge->ProjectionY(Form("hTrkY_%s",legName_charge[i].Data()),1+i*2,1+i*2);
      list->Add(hTrkYinCharge[i]);
    }
  c = drawHistos(list,"hTrkY_in_charge",Form("Au+Au %s: local y of projected tracks with different charges;y (cm);counts",trigName[kTrigType]),kFALSE,0,100,kTRUE,1e-1,6e6,kFALSE,kTRUE,legName_charge,kTRUE,"p_{T} > 1 GeV/c",0.4,0.6,0.65,0.85,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Projection/%s.track_y_InCharge_%s.png",run_config,trigName[kTrigType]));
  list->Clear();
  hTrkProjYZ->GetAxis(5)->SetRange(0,-1);
}


//================================================
void EnergyLoss(const Int_t save = 1)
{
  TFile *f1 = TFile::Open(Form("~/Work/STAR/analysis/output/jpsi.AuAu200.Run14.Match.Global.HLT.root"),"read");

  const char *name[2] = {"wEloss","woEloss"};
  THnSparseF *hTrkProjYZ[2]; 
  TH1F *hTrkY[2];

  // With energy loss
  hTrkProjYZ[0] = (THnSparseF*)f1->Get(Form("hTrkProjYZ_qa_%s",trigName[kTrigType]));
  hTrkProjYZ[0]->SetName(Form("%s_%s",hTrkProjYZ[0]->GetName(),name[0]));
  hTrkY[0] = (TH1F*)hTrkProjYZ[0]->Projection(0);

  // Without energy loss
  hTrkProjYZ[1] = (THnSparseF*)f->Get(Form("hTrkProjYZ_qa_%s",trigName[kTrigType]));
  hTrkProjYZ[1]->SetName(Form("%s_%s",hTrkProjYZ[1]->GetName(),name[1]));
  hTrkProjYZ[1]->GetAxis(5)->SetRange(2,2);
  hTrkY[1] = (TH1F*)hTrkProjYZ[1]->Projection(1);

  TList *list = new TList;
  TString legName[2] = {"With energy loss during projection","Without energy loss during projection"};
  for(Int_t i=0; i<2; i++)
    {
      hTrkY[i]->SetName(Form("hTrkY_%s",name[i]));
      hTrkY[i]->Scale(1./hTrkY[i]->Integral());
      list->Add(hTrkY[i]);
    }
  c = drawHistos(list,"hTrkY_Eloss",Form("Au+Au %s: local y of projected tracks at MTD radius;y (cm)",trigName[kTrigType]),kFALSE,0,100,kTRUE,1e-5,0.043,kFALSE,kTRUE,legName,kTRUE,"p_{T} > 1 GeV/c",0.15,0.6,0.7,0.88,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Projection/%s.track_y_Eloss_%s.png",run_config,trigName[kTrigType]));
}
