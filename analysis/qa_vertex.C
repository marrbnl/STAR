//================================================
void qa_vertex(const Int_t save = 1)
{
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.DzCut.root","read");
  TH1F *hCuts = (TH1F*)f->Get("mhAnalysisCuts");
  Double_t scale = hCuts->GetBinContent(3)/1e4;
  cout << "scale = " << scale << endl;


  TFile *f1 = TFile::Open("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.NoDzCut.root","read");
  TH2F *hDvzVsTpc[3];
  Double_t dzCut = hCuts->GetBinContent(11)/scale;
  for(Int_t i=0; i<3; i++)
    {
      hDvzVsTpc[i] = (TH2F*)f1->Get(Form("hDvzVsTpc_%s",trigName[i]));
      c = draw2D(hDvzVsTpc[i]);
      Double_t low_x = hDvzVsTpc[i]->GetXaxis()->GetBinLowEdge(1);
      Double_t up_x  = hDvzVsTpc[i]->GetXaxis()->GetBinUpEdge(hDvzVsTpc[i]->GetNbinsX());
      TLine *line = GetLine(low_x,-1*dzCut,up_x,-1*dzCut);
      line->Draw();
      line = GetLine(low_x,dzCut,up_x,dzCut);
      line->Draw();
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_vertex/DvzVsTpcVz_%s.png",trigName[i]));
    }

  TH1F *hVzTpc[3][2];
  TH1F *hVrTpc[3];
  TH2F *hVxVyTpc[3];
  TList *list = new TList;
  TString legName[3] = {"di-muon","single-muon","e-mu"};
  for(Int_t i=0; i<3; i++)
    {
      hVxVyTpc[i] = (TH2F*)f->Get(Form("hVxVyTpcWithCut_%s",trigName[i]));
      c = draw2D(hVxVyTpc[i]);
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_vertex/VxVyTpcWithCut_%s.png",trigName[i]));

      hVrTpc[i] = (TH1F*)f->Get(Form("hVrTpcWithCut_%s",trigName[i]));
      list->Add(hVrTpc[i]);
    }
  c = drawHistos(list,"hVrTpc","Au+Au: vr distribution of primary vertex",kFALSE,0,pi/2,kFALSE,-0.005,0.02,kTRUE,kTRUE,legName,kTRUE,Form("|vz_{TPC}-vz_{VPD}|<%1.0f cm",dzCut),0.5,0.68,0.6,0.85,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_vertex/VrTpcWithCut.png"));

  TString legName1[2] = {"w/o |vz_{TPC}-vz_{VPD}| cut","w   |vz_{TPC}-vz_{VPD}| cut"};
  for(Int_t i=0; i<3; i++)
    {
      hVzTpc[i][0] = (TH1F*)f->Get(Form("hVzTpc_%s",trigName[i]));
      hVzTpc[i][1] = (TH1F*)f->Get(Form("hVzTpcWithCut_%s",trigName[i]));
      list->Clear();
      list->Add(hVzTpc[i][0]);
      list->Add(hVzTpc[i][1]);
      c = drawHistos(list,hVzTpc[i][0]->GetName(),"Au+Au: vz distribution of primary vertex",kFALSE,0,pi/2,kFALSE,-0.005,0.02,kFALSE,kTRUE,legName1,kTRUE,trigName[i],0.15,0.3,0.6,0.85,kFALSE);
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_vertex/VzTpc_CompareCut_%s.png",trigName[i]));
    }
}
