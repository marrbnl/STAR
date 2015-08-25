TFile *f;
const char *run_config = "QA.";
const Bool_t iPico = 1;
const int year = 2014;
TString run_cfg_name;


//================================================
void ana_PartialTracking()
{
  gStyle->SetOptStat(0);
  TString fileName;

  if(year==2013)
    {
      run_type = "Run13_pp500";
      if(iPico) fileName = Form("Pico.Run13.pp500.jpsi.%sroot",run_config);
      else      fileName = Form("Run13.pp500.jpsi.%sroot",run_config);
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
      if(iPico) fileName = Form("Pico.Run14.AuAu200.jpsi.%sroot",run_config);
      else      fileName = Form("Run14.AuAu200.jpsi.%sroot",run_config);
    }

  f = TFile::Open(Form("./output/%s",fileName.Data()),"read");

  run_cfg_name = run_config;

  efficiency();
}

//================================================
void efficiency(const Int_t save = 0)
{
  TH1F *hEventCount = (TH1F*)f->Get("mPartialTrkEvents");
  for(Int_t ibin=1; ibin<=6; ibin++)
    {
      printf("%s: %d\n",hEventCount->GetXaxis()->GetBinLabel(ibin),hEventCount->GetBinContent(ibin));
    }
  c = draw1D(hEventCount,"",kTRUE,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PartialTracking/Event_Count.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PartialTracking/Event_Count.png",run_type));
    }

  const int nbins = 19;
  const double xbins[nbins+1] = {0,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.5,3.0,4.0,5.0,6.0,8,9,10};
  TH2F *hPartialTrkMuonPt = (TH2F*)f->Get("mPartialTrkMuonPt");
  TH1F *hMuonPt[2];
  TList *list = new TList();
  for(int i=0; i<2; i++)
    {
      TH1F *htmp = (TH1F*)hPartialTrkMuonPt->ProjectionX(Form("hMuonPt_%d_tmp",i),i+1,i+1);
      hMuonPt[i] = (TH1F*)htmp->Rebin(nbins,Form("hMuonPt_%d",i),xbins);
      scaleHisto(hMuonPt[i], 1, 1,kTRUE,kFALSE, kFALSE);
      list->Add(hMuonPt[i]);
    }
  TString legName[2] = {"All events with >=2 muons","Filtered events with >=2 muons"};
  c = drawHistos(list,"MuonPt","p_{T} distribution of muon candidates;p_{T} (GeV/c)",kFALSE,0,5,kFALSE,0.1,10,kTRUE,kTRUE,legName,kTRUE,"Au+Au @ 200 GeV",0.45,0.65,0.55,0.75);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PartialTracking/MuonPt.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PartialTracking/MuonPt.png",run_type));
    }

  TGraphAsymmErrors *ratio = new TGraphAsymmErrors();
  ratio->SetName("ratio");
  ratio->Divide(hMuonPt[1],hMuonPt[0],"cl=0.683 b(1,1) mode");
  ratio->SetMarkerStyle(21);
  ratio->GetYaxis()->SetRangeUser(0.98,1.01);
  c = drawGraph((TGraph*)ratio,"Efficiency of partial tracking;p_{T} (GeV/c);efficiency");
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PartialTracking/Efficiency_vs_pt.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PartialTracking/Efficiency_vs_pt.png",run_type));
    }

}
