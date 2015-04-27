TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;

const char *run_config = "DeltaZ.";
//const char *run_config = "CutRanking.";
const Bool_t iPico = 1;
TString run_cfg_name;

//================================================
void study_Match()
{						
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  

  TString cut_name = run_config;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  if(iPico)
    f = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.%sroot",run_config),"read");
  else
    f = TFile::Open(Form("./output/Run13.pp500.jpsi.%sroot",run_config),"read");

  if(run_config=="CutRanking.")
    run_cfg_name = Form("%s",run_config);
  else
    run_cfg_name = "";

  if(iPico) run_cfg_name = Form("Pico.%s",run_cfg_name.Data());

  StudyDeltaZ();

}

//================================================
void StudyDeltaZ(const Int_t save = 0)
{
  THnSparseF *hn = (THnSparseF*)f->Get(Form("mhTrkDzDyStudy_%s",trigName[kTrigType]));
  TList *list = new TList;

  // dz vs mtd z
  TH2F *hDzVsMtdZ = (TH2F*)hn->Projection(1,6);
  c = draw2D(hDzVsMtdZ,Form("%s: #Deltaz of matched track-hit pairs vs MTD hit z",trigName[kTrigType]));
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/study_Match/%sDeltaZ_vs_MtdZ_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/study_Match/%sDeltaZ_vs_MtdZ_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  list->Clear();
  TString legName[5];
  TH1F *hMtdZInMod[5];
  for(Int_t i=0; i<5; i++)
    {
      hn->GetAxis(4)->SetRange(i+1,i+1);
      hMtdZInMod[i] = (TH1F*)hn->Projection(6);
      hMtdZInMod[i]->SetName(Form("hMtdZ_Mod%d",i+1));
      legName[i] = Form("Module %d",i+1);
      //hMtdZInMod[i]->SetLineColor(color[i]);
      list->Add(hMtdZInMod[i]);
    }
  hn->GetAxis(4)->SetRange(0,-1);
  c = drawHistos(list,"MtdZInMod",Form("%s: #Deltaz of matched track-hit pairs in backleg;#Deltaz (cm)",trigName[kTrigType]),kFALSE,-100,100,kTRUE,1e-1,0.6*hMtdZInMod[3]->GetMaximum(),kFALSE,kTRUE,legName,kTRUE,"",0.15,0.25,0.6,0.88,kFALSE);

  TH1F *hBadMtdZInMod[5];
  hn->GetAxis(1)->SetRangeUser(50,200);
  for(Int_t i=0; i<5; i++)
    {
      hn->GetAxis(4)->SetRange(i+1,i+1);
      hBadMtdZInMod[i] = (TH1F*)hn->Projection(6);
      hBadMtdZInMod[i]->SetName(Form("hMtdZ_Mod%d",i+1));
      hBadMtdZInMod[i]->SetLineColor(color[i]);
      hBadMtdZInMod[i]->SetLineStyle(2);
      hBadMtdZInMod[i]->SetLineWidth(2);
      hBadMtdZInMod[i]->Draw("same");
    }
  hn->GetAxis(4)->SetRange(0,-1);

  //gPad->SetLogy();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sDeltaZ_in_BL_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sDeltaZ_in_BL_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }
  return;
}
