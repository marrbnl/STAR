const TString system = "Run14_AuAu200";
TFile *f = 0x0;

//================================================
void make_QA()
{
  gStyle->SetOptStat(0);

  f = TFile::Open(Form("Output/%s.jpsi.root",system.Data()),"read");
  if(!f->IsOpen())
    {
      printf("Can't open file %s\n",f->GetName());
      return;
    }
  cout << "+++ make QA histograms for file: " << f->GetName() << endl;

  const char* outPDFName = Form("PDF/%s.FinalQA.pdf",system.Data());
  const char* trigname = "di_mu";
  const char* trigtitle = "dimuon";

  TCanvas *c1 = new TCanvas("1pad","1pad",800,600);
  SetPadMargin(gPad);

  //----------------------------------------------------------------------------
  // Title page
  TPaveText *t1 = new TPaveText(0.28,0.5,0.7,0.7,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.07);
  t1->AddText("QA of MTD analysis");
  if(system=="Run14_AuAu200") t1->AddText("for Au+Au 200 GeV from Run14");
  t1->Draw();
  t1  = new TPaveText(0.28,0.25,0.7,0.4,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.05);
  TDatime time;
  t1->AddText(Form("%04d-%02d-%02d %02d:%02d:%02d",time.GetYear(),time.GetMonth(),time.GetDay(),
		   time.GetHour(),time.GetMinute(),time.GetSecond()));
  t1->Draw();
  c1->Print(Form("%s(",outPDFName));

  //----------------------------------------------------------------------------
  // Get TPDF
  TPDF *pdf = 0;
  TSeqCollection *col = gROOT->GetListOfSpecials();
  for(Int_t i=0; i<col->GetEntries(); i++)
    {
      TObject *obj = (TObject*)col->At(i);
      if( obj->IsA()->InheritsFrom("TPDF") && strcmp(obj->GetName(),outPDFName)==0 )
        pdf = dynamic_cast<TPDF*>obj;
    }
  if(!pdf) 
    {
      printf("No pointer to PDF file is available.\n");
      return;
    }
  pdf->Off();

  //----------------------------------------------------------------------------
  // Page 2
  title  = new TPaveText(0.1,0.85,0.9,0.95,"brNDC");
  title->SetFillStyle(0);
  title->SetBorderSize(0);
  title->SetTextFont(62);
  title->SetTextSize(0.07);
  title->AddText("Event selection");
  c1->Clear();
  title->Draw();
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  
  TH1F *hRun = (TH1F*)f->Get(Form("mhEvtRun_%s",trigname));
  int firstRun= 0;
  int lastRun;
  for(int bin=1; bin<=hRun->GetNbinsX(); bin++)
    {
      if(hRun->GetBinContent(bin)>0)
	{
	  int run = int(hRun->GetBinCenter(bin));
	  lastRun = run;
	  if(firstRun==0)
	    firstRun = run;
	  if(hRun->GetBinContent(bin)<1000) cout << run << endl;
	}
    }
  t1  = new TPaveText(0.05,0.6,0.5,0.85,"brNDC");
  t1->SetTextAlign(11);
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.045);
  t1->SetTextColor(2);
  t1->AddText(Form("Trigger: %s",trigtitle));
  t1->AddText(Form("Statistics: %4.3eM",(int)(hStat->GetBinContent(3))/1e6));
  t1->AddText(Form("Run range: [%d, %d]",firstRun,lastRun));
  t1->Draw();

  TH1F *hEventCuts = (TH1F*)f->Get("hAnalysisCuts");
  Int_t counter = (Int_t)(hEventCuts->GetBinContent(3)/1e4);
  t1  = new TPaveText(0.05,0.1,0.5,0.6,"brNDC");
  t1->SetTextAlign(11);
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.045);
  Double_t vtxz = hEventCuts->GetBinContent(1)/counter;
  Double_t dz = hEventCuts->GetBinContent(11)/counter;
  t1->AddText("Vertex selection:");
  t1->AddText(Form("     Select highest-ranked vertex that is within %1.0fcm",dz));
  t1->AddText("     of vpdVz when vpdVz is available; otherwise, the");
  t1->AddText("     default vertex is selected.");
  t1->AddText(Form("Vertex cut:"));
  t1->AddText(Form("     |tpcVz| < %1.0f cm",vtxz));
  t1->AddText(Form("     |tpcVz-vpdVz| < %1.0f cm",dz));
  t1->AddText(Form("Statistics after vertex cut: %4.3eM",(int)(hStat->GetBinContent(10))/1e6));
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  vector<TString> histoName;
  TH1F *h1 = 0x0;
  TH2F *h2 = 0x0;
  ///========================
  /// vertex
  TPaveText *t1 = new TPaveText(0.28,0.5,0.7,0.7,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.07);
  t1->AddText("Vertex");
  c1->Clear();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  histoName.clear();
  histoName.push_back(Form("mhTpcVx_%s",trigname));
  histoName.push_back(Form("mhTpcVy_%s",trigname));
  histoName.push_back(Form("mhTpcVz_%s",trigname));
  histoName.push_back(Form("mhTpcVr_%s",trigname));
  PaintHistosToPDF(histoName,pdf);

  histoName.clear();
  histoName.push_back(Form("mhTpcVyVx_%s",trigname));
  histoName.push_back(Form("mhTpcVxVz_%s",trigname));
  histoName.push_back(Form("mhTpcVyVz_%s",trigname));
  PaintHistosToPDF(histoName,pdf);

  histoName.clear();
  histoName.push_back(Form("mhVpdVz_%s",trigname));
  histoName.push_back(Form("mhVzDiffVsTpcVz_%s",trigname));
  PaintHistosToPDF(histoName,pdf);

  histoName.clear();
  histoName.push_back(Form("mhTpcVzWithCut_%s",trigname));
  histoName.push_back(Form("mhTpcVyVxWithCut_%s",trigname));
  histoName.push_back(Form("mhDiffVzWithCut_%s",trigname));
  PaintHistosToPDF(histoName,pdf);

  ///========================
  /// vertex
  TPaveText *t1 = new TPaveText(0.28,0.5,0.7,0.7,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.07);
  t1->AddText("Global property");
  c1->Clear();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  histoName.clear();
  histoName.push_back(Form("mhZdcRate_%s",trigname));
  histoName.push_back(Form("mhBbcRate_%s",trigname));
  histoName.push_back(Form("mhRefMult_%s",trigname));
  histoName.push_back(Form("mhTofMult_%s",trigname));
  PaintHistosToPDF(histoName,pdf);

  histoName.clear();
  histoName.push_back(Form("mhgRefMult_%s",trigname));
  histoName.push_back(Form("mhgRefMultCorr_%s",trigname));
  histoName.push_back(Form("mhCentrality_%s",trigname));
  histoName.push_back(Form("mhCentWeight_%s",trigname));
  PaintHistosToPDF(histoName,pdf);

  ///========================
  /// Primary tracks
  title  = new TPaveText(0.1,0.85,0.9,0.95,"brNDC");
  title->SetFillStyle(0);
  title->SetBorderSize(0);
  title->SetTextFont(62);
  title->SetTextSize(0.07);
  title->AddText("Track selection");
  t1  = new TPaveText(0.05,0.3,0.5,0.8,"brNDC");
  t1->SetTextAlign(11);
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.045);
  t1->AddText("Primary tracks");
  t1->AddText(Form("p_{T} > %1.1f GeV/c",hEventCuts->GetBinContent(2)/counter));
  t1->AddText(Form("|#eta| < %1.1f",hEventCuts->GetBinContent(4)/counter));  
  t1->AddText(Form("NHitsFit >= %1.0f",hEventCuts->GetBinContent(5)/counter));     
  t1->AddText(Form("NHitsDedx >= %1.0f",hEventCuts->GetBinContent(6)/counter));
  t1->AddText(Form("DCA < %1.1f cm",hEventCuts->GetBinContent(7)/counter));
  c1->Clear();
  title->Draw();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  histoName.clear();
  histoName.push_back(Form("mhTrkN_%s",trigname));
  histoName.push_back(Form("mhTrkPt_%s",trigname));
  histoName.push_back(Form("mhTrkEtaPhi_%s",trigname));
  PaintHistosToPDF(histoName,pdf);

  histoName.clear();
  histoName.push_back(Form("mhTrkDedx_%s",trigname));
  histoName.push_back(Form("mhTrkNsigmaPi_%s",trigname));
  PaintHistosToPDF(histoName,pdf,0.04);
  
  ///========================
  /// MTD hit
  TPaveText *t1 = new TPaveText(0.28,0.5,0.7,0.7,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.07);
  t1->AddText("MTD hits");
  c1->Clear();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  histoName.clear();
  histoName.push_back(Form("mhMtdHitN_%s",trigname));
  histoName.push_back(Form("mhMtdHitMap_%s",trigname));
  histoName.push_back(Form("mhMthMtdHitN_%s",trigname));
  histoName.push_back(Form("mhMthMtdHitMap_%s",trigname));
  c = DrawHistosOnCanvas(histoName,pdf);
  for(int i=0; i<histoName.size(); i++)
    {
      if(i==1 || i==3) continue;
      c->cd(i+1);
      h1 = (TH1F*)f->Get(histoName[i]);
      TPaveText *t1 = GetPaveText(0.6,0.8,0.6,0.8);
      t1->AddText(Form("<N> = %3.2f",h1->GetMean()));
      t1->Draw();
    }
  PrintCanvasToPDF(c,pdf);

  histoName.clear();
  histoName.push_back(Form("mhMtdTrigHitN_%s",trigname));
  histoName.push_back(Form("mhMtdTrigHitMap_%s",trigname));
  c = DrawHistosOnCanvas(histoName,pdf);
  for(int i=0; i<histoName.size(); i++)
    {
      if(i==1 || i==3) continue;
      c->cd(i+1);
      h1 = (TH1F*)f->Get(histoName[i]);
      TPaveText *t1 = GetPaveText(0.6,0.8,0.6,0.8);
      t1->AddText(Form("<N> = %3.2f",h1->GetMean()));
      t1->Draw();
    }
  PrintCanvasToPDF(c,pdf);

  ///========================
  /// muon PID
  title  = new TPaveText(0.1,0.85,0.9,0.95,"brNDC");
  title->SetFillStyle(0);
  title->SetBorderSize(0);
  title->SetTextFont(62);
  title->SetTextSize(0.07);
  title->AddText("Muon identification");
  
  t1  = new TPaveText(0.05,0.3,0.5,0.8,"brNDC");
  t1->SetTextAlign(11);
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.045);
  t1->AddText("Pimary track");
  t1->AddText("Fire MTD trigger");
  t1->AddText(Form("%1.1f < n#sigma_{#pi} < %1.1f",hEventCuts->GetBinContent(8)/counter,hEventCuts->GetBinContent(9)/counter));
  t1->AddText(Form("|#Deltay| < 2(2.5)#sigma for p_{T} <(>) 3 GeV/c"));
  t1->AddText(Form("|#Deltaz| < 2(2.5)#sigma for p_{T} <(>) 3 GeV/c"));
  t1->AddText(Form("#Deltatof < %1.2f ns",hEventCuts->GetBinContent(18)/counter));
  c1->Clear();
  title->Draw();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  histoName.clear();
  histoName.push_back(Form("mhMthTrkN_%s",trigname));
  histoName.push_back(Form("mhNMthTrkVsRun_%s",trigname));
  histoName.push_back(Form("mhMthTrkPt_%s",trigname));
  histoName.push_back(Form("mhMthTrkPtVsRun_%s",trigname));
  PaintHistosToPDF(histoName,pdf);

  histoName.clear();
  histoName.push_back(Form("mhMthTrkEtaPhi_%s",trigname));
  histoName.push_back(Form("mhMthTrkEtaVsRun_%s",trigname));
  histoName.push_back(Form("mhMthTrkPhiVsRun_%s",trigname));
  PaintHistosToPDF(histoName,pdf);

  histoName.clear();
  histoName.push_back(Form("mhMthTrkNsigmaPi_%s",trigname));
  histoName.push_back(Form("mhMthTrkNsigmaPiVsRun_%s",trigname));
  histoName.push_back(Form("mhDtofVsPt_%s",trigname));
  histoName.push_back(Form("mhMthTrkDtofVsRun_%s",trigname));
  PaintHistosToPDF(histoName,pdf);

  histoName.clear();
  histoName.push_back(Form("mhDyVsPt_%s",trigname));
  histoName.push_back(Form("mhMthTrkDyVsRun_%s",trigname));
  histoName.push_back(Form("mhDzVsPt_%s",trigname));
  histoName.push_back(Form("mhMthTrkDzVsRun_%s",trigname));
  PaintHistosToPDF(histoName,pdf);

  histoName.clear();
  histoName.push_back(Form("mhNMuonCandidate_%s",trigname));
  histoName.push_back(Form("mhMuonPhiVsPt_%s",trigname));
  histoName.push_back(Form("mhNMuonPosVsNeg_%s",trigname));
  PaintHistosToPDF(histoName,pdf);

  histoName.clear();
  histoName.push_back(Form("mhNMuonVsRun_%s",trigname));
  histoName.push_back(Form("mhMuonPtVsRun_%s",trigname));
  histoName.push_back(Form("mhMuonEtaVsRun_%s",trigname));
  histoName.push_back(Form("mhMuonPhiVsRun_%s",trigname));
  PaintHistosToPDF(histoName,pdf);

  pdf->On();
  pdf->Close();
}



//-----------------------------------------
void PaintHistosToPDF(vector<TString> histoName, TPDF *pdf, const float size = 0.045)
{
  TCanvas *c = DrawHistosOnCanvas(histoName,pdf,size);
  PrintCanvasToPDF(c,pdf);
}

//-----------------------------------------
TCanvas *DrawHistosOnCanvas(vector<TString> histoName, TPDF *pdf, const float size = 0.045)
{
  TCanvas *c;
  int nHisto = histoName.size();
  c = new TCanvas(histoName[0],histoName[0],800,600);
  if(nHisto==1)
    {
      //c = new TCanvas(histoName[0],histoName[0],800,600);
    }
  else if(nHisto<=4)
    {
      //c = new TCanvas(histoName[0],histoName[0],1200,750);

      c->Divide(2,2);
    }
  else if(nHisto<=8)
    {
      //c = new TCanvas(histoName[0],histoName[0],1400,750);
      c->Divide(4,2);
    }

  for(int i=0; i<nHisto; i++)
    {
      TObject *obj = f->Get(histoName[i]);
      c->cd(i+1);
      SetPadMargin(gPad,0.12,0.12,0.1,0.1);
      gPad->SetGrid(0,0);
      
      TH1 *h1 = dynamic_cast<TH1*>obj;
      TPaveText *t = GetTitleText(h1->GetTitle(),size);
      h1->SetTitle("");
      //ScaleHistoTitle(h1);
      if(h1->InheritsFrom("TH2"))
	{
	  TH2 *h2 = (TH2*)h1;
	  gPad->SetLogz();
	  h2->Draw("colz");
	}
      else if(h1->InheritsFrom("TProfile"))
	{
	  h1->SetMarkerStyle(24);
	  h1->Draw();
	}
      else
	{
	  gPad->SetLogy();
	  h1->Draw();
	}
      t->Draw();
    }
  return c;
}
