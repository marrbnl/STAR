const char *triggerName[3] = {"di-muon","single-muon","e-mu"};

//================================================
void vertex(const char *day = "077", const Int_t save = 1)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.15);                
  gStyle->SetStatH(0.13); 
  TFile *f = TFile::Open(Form("Output/output.%s.root",day),"read");

  TH2F *hVpdVsTpc[3];
  TH2F *hDiffVsTpc[3];
  TH1F *hDiff[3];
  TList *list = new TList;
  TString legName[3];
  for(Int_t i=0; i<3; i++)
    {
      hVpdVsTpc[i] = (TH2F*)f->Get(Form("hVtxVpdTpc_%s",triggerName[i]));
      c = draw2D(hVpdVsTpc[i],Form("Vertex z: VPD vs TPC (%s trigger);vz_{TPC} (cm);vz_{VPD} (cm)",triggerName[i]));
      if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/Vertex/%s_VertexZ_VPD_vs_TPC_%s.png",day,triggerName[i]));

      hDiffVsTpc[i] = (TH2F*)f->Get(Form("hVtxDiff_%s",triggerName[i]));
      c = draw2D(hDiffVsTpc[i],Form("Vertex z: TPC-VPD vs TPC (%s trigger);vz_{TPC} (cm);vz_{TPC}-vz_{VPD} (cm)",triggerName[i]));
      if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/Vertex/%s_VertexZ_TPC-VPD_vs_TPC_%s.png",day,triggerName[i]));

      hDiff[i] = (TH1F*)hDiffVsTpc[i]->ProjectionY(Form("hDiff_%s",triggerName[i]));
      list->Add(hDiff[i]);
      legName[i] = triggerName[i];
    }
  c = drawHistos(list,"hVtxDiff_TpcVpd","Vertex z difference between TPC and VPD",kTRUE,-10,10,kFALSE,0.01,1.2e5,kTRUE,kTRUE,legName,kTRUE,"Trigger",0.2,0.4,0.6,0.85,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/Vertex/%s_VertexZ_TPC-VPD.png",day));

  // for(Int_t i=0; i<3; i++)
  //   {
  //     hDiff[i]->GetXaxis()->SetRangeUser(-10,10);
  //     hDiff[i]->SetLineColor(2);
  //     draw1D(hDiff[i],Form("%s: Vertex z difference between TPC and VPD",triggerName[i]),kFALSE,kFALSE);
  //     TF1 *func = new TF1(Form("func_%d",i),"gaus",-2,2);
  //     hDiff[i]->Fit(func,"IR");
  //     func->Draw("sames");
  //   }
}
