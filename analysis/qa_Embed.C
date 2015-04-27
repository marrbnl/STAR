const char *run_config = "EmbedQA.";
const Bool_t iPico = 1;
const int year = 2013;
TString run_cfg_name;

TFile *fdata, *fmc;

//================================================
void qa_Embed()
{
  gStyle->SetOptStat(0);

  if(year==2013)
    {
      run_type = "Run13_pp500";
      fmc   = TFile::Open(Form("./output/Run13.pp500.jpsi.%sMC.root",run_config),"read");
      fdata = TFile::Open(Form("./output/Run13.pp500.jpsi.%sPicoData.root",run_config),"read");
    }
  run_cfg_name = Form("%s",run_config);
  if(iPico) run_cfg_name = Form("Pico.%s",run_cfg_name.Data());

  makePDF();
}


//================================================
void makePDF(char *outPDFName="")
{
  if(year==2013) outPDFName = "EmbedQA_Run13_pp500_Jpsi.pdf";
  TDatime time;
  Int_t year  = time.GetYear();
  Int_t month = time.GetMonth();
  Int_t day   = time.GetDay();

  TCanvas *c1 = new TCanvas("1pad","1pad",800,600);
  SetPadMargin(gPad);

  //----------------------------------------------------------------------------
  // Title page
  TPaveText *t1 = new TPaveText(0.28,0.5,0.7,0.7,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.07);
  t1->AddText("Embedding QA for J/#psi #rightarrow #mu^{+}#mu^{-}");
  t1->AddText("in pp 500 GeV from Run13");
  t1->Draw();
  t1  = new TPaveText(0.28,0.25,0.7,0.4,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.05);
  t1->AddText(Form("%02d/%02d/%04d",month,day,year));
  t1->AddText("Rongrong Ma");
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
  // analysis cuts
  title  = new TPaveText(0.1,0.8,0.9,0.9,"brNDC");
  title->SetFillStyle(0);
  title->SetBorderSize(0);
  title->SetTextFont(62);
  title->SetTextSize(0.07);
  title->AddText("Analysis cuts");
  
  t1  = new TPaveText(0.05,0.1,0.5,0.8,"brNDC");
  t1->SetTextAlign(11);
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.045);

  TH1F *h1 = (TH1F*)fmc->Get("hAnalysisCuts");
  Int_t counter = (Int_t)(h1->GetBinContent(3)/1e4);
  Double_t vtxz = h1->GetBinContent(1)/counter;
  if(TMath::Abs(vtxz)<1e3)
    {
      t1->AddText(Form("Vertex cut: |vtx_z| < %1.0f cm",vtxz));
    }
  else
    {
      t1->AddText(Form("Vertex cut: none"));
    }

  t1->AddText(Form("Track selection:"));
  t1->AddText("    Primary tracks");
  t1->AddText(Form("    p_{T} > %1.1f GeV/c",h1->GetBinContent(2)/counter));
  t1->AddText(Form("    |#eta| < %1.1f",h1->GetBinContent(4)/counter));  
  t1->AddText(Form("    NHitsFit > %1.0f",h1->GetBinContent(5)/counter));     
  t1->AddText(Form("    NHitsDedx > %1.0f",h1->GetBinContent(6)/counter));
  t1->AddText(Form("    NHitsFit/NHitsPoss > %1.2f",h1->GetBinContent(10)/counter));
  if(TMath::Abs(h1->GetBinContent(7)/counter)<10)
    t1->AddText(Form("    global dca < %1.1f cm",h1->GetBinContent(7)/counter));  
  c1->Clear();
  title->Draw();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  title  = new TPaveText(0.1,0.8,0.9,0.9,"brNDC");
  title->SetFillStyle(0);
  title->SetBorderSize(0);
  title->SetTextFont(62);
  title->SetTextSize(0.07);
  title->AddText("PID cuts in data");
  
  t1  = new TPaveText(0.05,0.1,0.5,0.8,"brNDC");
  t1->SetTextAlign(11);
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.045);
  t1->AddText("Pion PID");
  t1->AddText("    0.016 < m^{2} < 0.021 (GeV/c^{2})^{2}");
  t1->AddText("    |n#sigma_{#pi}| < 4");
  t1->AddText("Muon PID");
  t1->AddText(Form("    %1.1f < n#sigma_{#pi} < %1.1f",h1->GetBinContent(8)/counter,h1->GetBinContent(9)/counter));
  t1->AddText(Form("    |#Deltaz| < %1.0f cm",h1->GetBinContent(14)/counter));
  t1->AddText(Form("    Under J/#psi mass peak [3.0,3.2] GeV/c^{2}));
  c1->Clear();
  title->Draw();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);


  pdf->On();
  pdf->Close();
}
