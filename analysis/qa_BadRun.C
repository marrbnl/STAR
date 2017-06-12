const int year = YEAR;
TFile *f;

//================================================
void qa_All()
{
  gStyle->SetOptStat(0);
  if(year==2013)
    {
      f = TFile::Open("./output/Run13.pp500.jpsi.RunQA.root","read");
    }
  else if(year==2014)
    {
      f = TFile::Open("output/Pico.Run14.AuAu200.jpsi.root ","read");
    }

  findBadRuns();
}

//================================================
void findBadRuns(const Bool_t save = 0)
{
  char outPDFName[256];
  sprintf(outPDFName,"PDF/%s.RunDepQA.pdf",run_type);

  //----------------------------------------------------------------------------
  // Title page
  TDatime time;
  Int_t year  = time.GetYear();
  Int_t month = time.GetMonth();
  Int_t day   = time.GetDay();

  TCanvas *c1 = new TCanvas("1pad","1pad",800,600);
  SetPadMargin(gPad);


  TPaveText *t1 = new TPaveText(0.28,0.6,0.7,0.7,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.06);
  t1->AddText("MTD run-wise QA");
  t1->Draw();

  t1 = new TPaveText(0.28,0.5,0.7,0.6,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.06);
  if(year==2013) t1->AddText("for pp 500 GeV in Run13");
  if(year==2014) t1->AddText("for AuAu 200 GeV in Run14");
  t1->Draw();

  t1  = new TPaveText(0.28,0.3,0.7,0.4,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.05);
  t1->AddText(Form("%02d/%02d/%04d",month,day,year));
  t1->Draw();
  c1->Print(Form("%s(",outPDFName));

  // ----------------------------------------------------------------------------
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
  title  = new TPaveText(0.1,0.85,0.9,0.95,"brNDC");
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

  t1->AddText("Trigger: di-muon");

  TH1F *h1 = (TH1F*)fmc->Get("hAnalysisCuts");
  Int_t counter = (Int_t)(h1->GetBinContent(3)/1e4);
  Double_t vtxz = h1->GetBinContent(1)/counter;
  if(1)
    {
      t1->AddText(Form("Vertex cut: |vz_{TPC}| < %d cm, |vz_{TPC}-vz_{VPD}| < 3 cm, vr < 2 cm",100));
    }
  else
    {
      t1->AddText(Form("Vertex cut: none"));
    }

  t1->AddText(Form("Track selection:"));
  t1->AddText("    Primary tracks");
  //t1->AddText(Form("    p_{T} > %1.1f GeV/c",h1->GetBinContent(2)/counter));
  t1->AddText(Form("    p_{T} > 1 GeV/c"));
  t1->AddText(Form("    |#eta| < %1.1f",h1->GetBinContent(4)/counter));  
  t1->AddText(Form("    NHitsFit >= %1.0f",h1->GetBinContent(5)/counter));     
  t1->AddText(Form("    NHitsDedx >= %1.0f",h1->GetBinContent(6)/counter));
  t1->AddText(Form("    NHitsFit/NHitsPoss >= %1.2f",h1->GetBinContent(10)/counter));
  if(TMath::Abs(h1->GetBinContent(7)/counter)<10)
    t1->AddText(Form("    global dca <= %1.1f cm",h1->GetBinContent(7)/counter));  
  c1->Clear();
  title->Draw();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  title  = new TPaveText(0.1,0.85,0.9,0.95,"brNDC");
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
  if(year==2013) 
    t1->AddText("    0.016 < m^{2} < 0.021 (GeV/c^{2})^{2}");
  t1->AddText("    |n#sigma_{#pi}| < 2");
  t1->AddText("Muon PID");
  t1->AddText(Form("    %1.1f < n#sigma_{#pi} < %1.1f",h1->GetBinContent(8)/counter,h1->GetBinContent(9)/counter));
  t1->AddText(Form("    |#Deltaz| < %1.0f cm",h1->GetBinContent(14)/counter));
  t1->AddText(Form("    Under J/#psi mass peak of [3.0,3.2] GeV/c^{2}"));
  c1->Clear();
  title->Draw();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);


  // ----------------------------------------------------------------------------
  // Histograms
  TH1F *hStat = (TH1F*)f->Get(Form("mhNeventVsRun_%s",trigName[kTrigType]));
  c = draw1D(hStat,"",0.04,kFALSE,kFALSE);
  t1 = GetPaveText(0.25,0.4,0.83,0.88);
  t1->AddText(Form("Total # of events = %e",hStat->GetEntries()));
  t1->SetTextFont(62);
  t1->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Run/NeventVsRun_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Run/NeventVsRun_%s.png",run_type,trigName[kTrigType]));
    }
  PaintCanvasToPDF(c, pdf);

  TList *list = (TList*)f->GetListOfKeys();
  Int_t nEntry = list->GetEntries();
  Int_t runCandidate[200], index[200];
  Int_t nCandidate = 0, nVarialbes = 0;
  for(Int_t i=0; i<nEntry; i++)
    {
      TKey *key = (TKey*)list->At(i);
      TObject *obj = key->ReadObj();
      if(!obj->IsA()->InheritsFrom("TH2")) continue;
      TString objname = obj->GetName();
      TH2F *h2 = (TH2F*)obj;

      if(objname.Contains("mhRefMultVsRun"))
	{
	  h2->GetYaxis()->SetRangeUser(0,50);
	}
      else if(objname.Contains("mhBBCrateVsRun"))
	{
	  h2->GetYaxis()->SetRangeUser(0,7);
	}
      else if(objname.Contains("mhZDCrateVsRun"))
	{
	  h2->GetYaxis()->SetRangeUser(0,0.6);
	}
      else if(objname.Contains("mhNsigmaPiVsRun"))
	{
	  h2->GetYaxis()->SetRangeUser(-6,10);
	}
      else if(objname.Contains("mhTrkPtVsRun"))
	{
	  h2->GetYaxis()->SetRangeUser(0,10);
	}

      TString title = h2->GetTitle();
      c = draw2D(h2,"",0.04,kTRUE);
      TString outname = objname;
      outname.Remove(outname.Index("_"));
      outname.ReplaceAll("mh","");
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Run/%s_%s.pdf",run_type,outname.Data(),trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Run/%s_%s.png",run_type,outname.Data(),trigName[kTrigType]));
	}

      if(!objname.Contains("mhBBCrateVsRun") && !objname.Contains("mhZDCrateVsRun"))
	{
	  nVarialbes ++;
	  TProfile *hpro = (TProfile*)h2->ProfileX(Form("%s_pfx",objname.Data()));
	  hpro->SetMarkerStyle(21);

	  Int_t nRun = 0;
	  Double_t variable[1000];
	  for(Int_t ibin=1; ibin<=hpro->GetNbinsX(); ibin++)
	    {
	      if(hpro->GetBinContent(ibin)==0 && hpro->GetBinError(ibin)==0) continue;
	      variable[nRun] = hpro->GetBinContent(ibin);
	      nRun ++;
	    }
	  Double_t rms  = TMath::RMS(nRun,variable);
	  Double_t mean = TMath::Mean(nRun,variable);

	  TH1F *hOutlier = new TH1F(Form("%s_outlier",objname.Data()),hpro->GetTitle(),hpro->GetNbinsX(),0,hpro->GetNbinsX());
	  for(Int_t ibin=1; ibin<=hpro->GetNbinsX(); ibin++)
	    {
	      if(hpro->GetBinContent(ibin)!=0 &&
		 hpro->GetBinError(ibin)!=0
		 && TMath::Abs(hpro->GetBinContent(ibin)-mean)>5*rms
		 && TMath::Abs(hpro->GetBinContent(ibin)-mean)>5*TMath::Abs(hpro->GetBinError(ibin)))
		{
		  hOutlier->SetBinContent(ibin,hpro->GetBinContent(ibin));
		  hOutlier->SetBinError(ibin,hpro->GetBinError(ibin));
		  cout << objname.Data() << " -> " << ibin << endl;
		  runCandidate[nCandidate] = ibin - 1;
		  nCandidate ++;
		}
	      else
		{
		  hOutlier->SetBinContent(ibin,999);
		  hOutlier->SetBinError(ibin,1e-10);
		}
 
	    }
	  hOutlier->SetMarkerStyle(21);
	  hOutlier->SetMarkerColor(2);
	  hOutlier->SetLineColor(2);
	  c->cd();
	  hpro->Draw("sames");
	  hOutlier->Draw("sameP");

	  hpro->GetYaxis()->SetRangeUser(mean-10*rms, mean + 10*rms);
	  c1 = draw1D(hpro,Form("%s;%s;%s",title.Data(),h2->GetXaxis()->GetTitle(),h2->GetYaxis()->GetTitle()));
	  c1->cd();
	  hOutlier->Draw("sameP");
	  if(save) 
	    {
	      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Run/%s_pro_%s.pdf",run_type,outname.Data(),trigName[kTrigType]));
	      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Run/%s_pro_%s.png",run_type,outname.Data(),trigName[kTrigType]));
	    }
	}

      PaintCanvasToPDF(c,pdf);
    }
  pdf->On();
  pdf->Close();

  TMath::Sort(nCandidate,runCandidate,index,kFALSE);
  vector<Int_t> badRunList, counter;
  badRunList.clear();
  for(Int_t i=0; i<nCandidate; i++)
    {
      Int_t run = runCandidate[index[i]];
      Bool_t exist = kFALSE;
      for(Int_t j=0; j<badRunList.size(); j++)
	{
	  if(run == badRunList[j])
	    {
	      exist = kTRUE;
	      break;
	    }
	}
      if(!exist)
	{
	  badRunList.push_back(runCandidate[index[i]]);
	  Int_t count = 0;
	  for(Int_t ii=0; ii<nCandidate; ii++)
	    {
	      if(run==runCandidate[index[ii]]) count++;
	    }
	  counter.push_back(count);
	}
    }

  Int_t nBadRun = badRunList.size();
  printf("================================\n");
  printf("%d bad runs found!\n",nBadRun);
  if(nBadRun>-1)
    {
      for(Int_t i=0; i<nBadRun; i++)
	{
	  cout << runList[badRunList[i]] << ", ";
	}
      cout << endl;

      for(Int_t i=0; i<nBadRun; i++)
	{
	  cout << runList[badRunList[i]] << "  " 
	       << counter[i] << "/" << nVarialbes << "  " 
	       << hStat->GetBinContent(badRunList[i]+1)
	       << endl;
	}
    }
  printf("================================\n");
  
}
