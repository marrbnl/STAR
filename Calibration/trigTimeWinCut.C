
//================================================
void trigTimeWinCut()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  //gStyle->SetStatY(0.9);                
  //gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.28);                
  gStyle->SetStatH(0.28); 
  gStyle->SetTitleSize(0.06,"T");

  //massProduction();
  //fit();
  //compareDays();
  // example();
}

//================================================
void massProduction()
{
  fit("077",1,1);
  fit("088",1,1);
  fit("089",1,1);
  fit("121",1,1);
  fit("141",1,1);
  fit("136",1,1);
  fit("all",1,1);
}


//================================================
void fit(const char * day = "all", const Int_t saveDB = 0, const Int_t saveHisto = 0)
{
  char inFile[256];
  char outFile[256];
  if(day=="all")
    {
      sprintf(inFile,"Output/histos.AuAu200.root");
      sprintf(outFile,"AuAu200_%s",day);
    }
  else
    {
      sprintf(inFile,"Output/output.%s.root",day);
      sprintf(outFile,"AuAu200_%s",day);
    }
  
  TFile *fin = TFile::Open(inFile,"read");

  TH2F *h2 = (TH2F*)fin->Get("hMtdHitTrigTimeWest");
  TH1F *hTrigTime[30][5];
  TF1 *func[30][5];
  TCanvas *c[6];
  TH1F *hCutLow   = new TH1F("trigTimeWin_lowCut","trigTimeWin_lowCut;backleg*5+tray",150,0,150);
  TH1F *hCutHigh  = new TH1F("trigTimeWin_lowHigh","trigTimeWin_lowHigh;backleg*5+tray",150,0,150);
  TH1F *hCutLow2  = new TH1F("trigTimeWin_lowCut_perChan","trigTimeWin_lowCut_perChannel;global channel",1800,0,1800);
  TH1F *hCutHigh2 = new TH1F("trigTimeWin_lowHigh_perChan","trigTimeWin_lowHigh_perChannel;global channel",1800,0,1800);
  TH1F *hMean     = new TH1F("trigTimeWin_mean","Mean of the time distribution of MTH hits;backleg*5+tray",150,0,150);
  TH1F *hSigma    = new TH1F("trigTimeWin_sigma","Sigma of the time distribution of MTH hits;backleg*5+tray",150,0,150);
  TSpectrum *peakfinder = new TSpectrum();
  for(Int_t i=0; i<30; i++)
    {
      //if(i>4) continue;
      if(i%5==0)
	{
	  c[i/5] = new TCanvas(Form("Backleg%d-%d",i+1,i+5),Form("Backleg%d-%d",i+1,i+5),1350,800);
	  c[i/5]->Divide(5,5);
	}
      for(Int_t j=0; j<5; j++)
	{
	  //printf("Back = %d, tray = %d, channel = %d - %d\n",i+1,j+1,i*60+j*12+1,i*60+(j+1)*12);
	  hTrigTime[i][j] = (TH1F*)h2->ProjectionY(Form("MTDhits_triggerTime_BL%d_Tray%d",i+1,j+1),i*60+j*12+1,i*60+(j+1)*12);
	  c[i/5]->cd(i%5*5+j+1);
	  if(i<15) hTrigTime[i][j]->GetXaxis()->SetRangeUser(2800,3000);
	  else     hTrigTime[i][j]->GetXaxis()->SetRangeUser(2850,3100);
	  hTrigTime[i][j]->SetFillColor(3);
	  hTrigTime[i][j]->SetTitle(Form("BL %d - tray %d",i+1,j+1));
	  hTrigTime[i][j]->Draw("");
	  if(hTrigTime[i][j]->GetEntries()<1) continue;

	  TH1F *htmp = (TH1F*)hTrigTime[i][j]->Clone(Form("%s_peak",hTrigTime[i][j]->GetName()));
	  Int_t nfound = peakfinder->Search(htmp,1,"new");
	  TList *functions = htmp->GetListOfFunctions();
	  TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
	  Double_t *xpeaks = pm->GetX();
	  printf("Found %d candidate peaks, mu = %2.4f\n",nfound,xpeaks[0]);

	  Double_t isigma = 10;
	  Double_t min = xpeaks[0] - 3*isigma;
	  Double_t max = xpeaks[0] + 1*isigma;
	  Double_t imean = xpeaks[0];

	  func[i][j] = new TF1(Form("Fit_triggerTime_BL%d_Tray%d",i+1,j+1),"gaus");
	  func[i][j]->SetLineColor(2);
	  func[i][j]->SetLineWidth(1);
	  func[i][j]->SetRange(min,max);
	  func[i][j]->SetParameter(1,imean);
	  func[i][j]->SetParameter(2,isigma);
	  func[i][j]->SetParName(0,"C");
	  func[i][j]->SetParName(1,"#mu");
	  func[i][j]->SetParName(2,"#sigma");
	  hTrigTime[i][j]->Fit(func[i][j],"IRQ");

	  Double_t mean = func[i][j]->GetParameter(1);
	  Double_t sigma = func[i][j]->GetParameter(2);
	  Double_t yield = func[i][j]->GetParameter(0)*0.8;
	  Int_t nsigma = 3;

	  TLine *line = GetLine(mean-nsigma*sigma,0,mean-nsigma*sigma,yield,4);
	  line->Draw("same");
	  line = GetLine(mean+nsigma*sigma,0,mean+nsigma*sigma,yield,4);
	  line->Draw();
	  hCutLow->SetBinContent(i*5+j+1,mean-nsigma*sigma);
	  hCutHigh->SetBinContent(i*5+j+1,mean+nsigma*sigma);
	  hMean->SetBinContent(i*5+j+1,func[i][j]->GetParameter(1));
	  //hMean->SetBinError(i*5+j+1,func[i][j]->GetParError(1));
	  hSigma->SetBinContent(i*5+j+1,func[i][j]->GetParameter(2));
	  //hSigma->SetBinError(i*5+j+1,func[i][j]->GetParError(2));

	  for(Int_t k = i*60+j*12+1; k<=i*60+(j+1)*12; k++)
	    {
	      hCutLow2->SetBinContent(k,mean-nsigma*sigma);
	      hCutHigh2->SetBinContent(k,mean+nsigma*sigma);
	    }
	}
    }

  c2 = draw2D(h2);
  c2->cd();
  hCutLow2->SetLineColor(6);
  hCutLow2->SetLineWidth(2);
  hCutHigh2->SetLineColor(6);
  hCutHigh2->SetLineWidth(2);
  hCutLow2->Draw("sames PL");
  hCutHigh2->Draw("sames PL");
  TPaveText *t1 = GetPaveText(0.2,0.6,0.78,0.85);
  t1->SetTextColor(6);
  t1->SetFillStyle(1);
  t1->AddText(Form("Magenta lines indicate %d#sigma cut",nsigma));
  t1->Draw();

  hMean->SetMarkerStyle(20);
  hMean->GetYaxis()->SetRangeUser(2800,3000);
  c11 = draw1D(hMean,"",kFALSE,kTRUE,0.04,"PL");
  hSigma->SetMarkerStyle(20);
  hSigma->GetYaxis()->SetRangeUser(7,12);
  c12 = draw1D(hSigma,"",kFALSE,kTRUE,0.04,"PL");

  if(saveHisto)
    {
      c2->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/AuAu200_%s_TrigTimeWinCut.png",day));
      c11->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/AuAu200_%s_Fit_TrigTime_Mean.png",day));
      c12->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/AuAu200_%s_Fit_TrigTime_Sigma.png",day));
      for(Int_t i=0; i<6; i++)
	c[i]->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/AuAu200_%s_Fit_TrigTime_BL%d-%d.png",day,i*5+1, i*5+5));
    }

  // print out
  cout << "Low limit for trigger time window cut: " << endl;
  for(Int_t i=0; i<30; i++)
    {
      for(Int_t j=0; j<5; j++)
	{
	  cout << hCutLow->GetBinContent(i*5+j+1) << ", "; 
	}
      cout << endl;
    } 

  cout << "Up limit for trigger time window cut: " << endl;
  for(Int_t i=0; i<30; i++)
    {
      for(Int_t j=0; j<5; j++)
	{
	  cout << hCutHigh->GetBinContent(i*5+j+1) << ", "; 
	}
      cout << endl;
    }

  if(saveDB)
    {
      //==============================================
      // root file
      TFile *fout = TFile::Open(Form("DB/%s_trigTimeWinCut.root",outFile),"recreate");
      hCutLow->Write();
      hCutLow2->Write();
      hCutHigh->Write();
      hCutHigh2->Write();
      hMean->Write();
      hSigma->Write();
      for(Int_t i=0; i<30; i++)
	{
	  for(Int_t j=0; j<5; j++)
	    {
	      if(hTrigTime[i][j]) hTrigTime[i][j]->Write();
	      if(func[i][j])      func[i][j]->Write();
	    }
	}
      fout->Close();
      
      // PDF file
      char *outPDFName = Form("DB/%s_trigTimeWinCut.pdf",outFile);
      c2->Print(Form("%s(",outPDFName));
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
      for(Int_t i=0; i<6; i++)
	PaintCanvasToPDF(c[i],pdf);
      PaintCanvasToPDF(c11,pdf);
      PaintCanvasToPDF(c12,pdf);
      pdf->On();
      pdf->Close();
      
      // database
      ofstream out;
      out.open(Form("DB/%s_trigTimeWinCut.dat",outFile),std::ofstream::out);
      for(Int_t i=0; i<30; i++)
	{
	  for(Int_t j=0; j<5; j++)
	    {
	      out << i+1 << "  " << j+1 << "  " << "  " << hCutLow->GetBinContent(i*5+j+1) << "  " << hCutHigh->GetBinContent(i*5+j+1) << endl;
	    }
	}  
      out.close();
    }
}


//================================================
void example(const Int_t saveHisto = 1)
{ 
  const Int_t bl = 1, tray = 1;
  TFile *f = TFile::Open(Form("DB/AuAu200_all_trigTimeWinCut.root"),"read");
  TH1F *h = (TH1F*)f->Get(Form("MTDhits_triggerTime_BL%d_Tray%d",bl,tray));
  TCanvas *c = new TCanvas(h->GetName(),h->GetName(),800,600);
  h->Draw();
  TF1 *func = (TF1*)f->Get(Form("Fit_triggerTime_BL%d_Tray%d",bl,tray));
  func->SetLineColor(2);
  func->SetLineWidth(2);
  func->Draw("sames");
  Double_t mean = func->GetParameter(1);
  Double_t sigma = func->GetParameter(2);
  Double_t yield = func->GetParameter(0)*0.8;
  Int_t nsigma = 3;

  TLine *line = GetLine(mean-nsigma*sigma,0,mean-nsigma*sigma,yield,4);
  line->Draw("same");
  line = GetLine(mean+nsigma*sigma,0,mean+nsigma*sigma,yield,4);
  line->Draw();
  if(saveHisto) c->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/AuAu200_TrigTimeWinCut_BL%d_Tray%d.png",bl,tray));
  
}

//================================================
void compareDays(const Int_t saveHisto = 1)
{
  const Int_t nDays = 7;
  const char *days[nDays] = {"all","077","088","089","121","141","136"};
  TH1F *hTrigCuts[2][7];
  TLegend *leg = new TLegend(0.15,0.55,0.25,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);

  TLegend *leg2 = new TLegend(0.35,0.55,0.45,0.85);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextSize(0.04);

  for(Int_t i=0; i<nDays; i++)
    {
      TFile *f = TFile::Open(Form("DB/AuAu200_%s_trigTimeWinCut.root",days[i]),"read");
      hTrigCuts[0][i] = (TH1F*)f->Get("trigTimeWin_lowCut");
      hTrigCuts[1][i] = (TH1F*)f->Get("trigTimeWin_lowHigh");
      for(Int_t j=0; j<2; j++)
	{
	  hTrigCuts[j][i]->SetName(Form("%s_day_%s",hTrigCuts[j][i]->GetName(),days[i]));
	  hTrigCuts[j][i]->SetMarkerStyle(21+j*4);
	  hTrigCuts[j][i]->SetMarkerColor(color[nDays-i-1]);
	  hTrigCuts[j][i]->SetLineColor(color[nDays-i-1]);
	  hTrigCuts[j][i]->GetYaxis()->SetRangeUser(2800,3000);
	  hTrigCuts[j][i]->GetYaxis()->SetNdivisions(505);
	  hTrigCuts[j][i]->GetYaxis()->SetTitleOffset(1.5);
	  if(i==0 && j==0) 
	    {
	      c = draw1D(hTrigCuts[j][i],"Trigger window cuts for MTD hits;backleg*5+tray;TAC_{MTD}-TAC_{THUB}",kFALSE,kTRUE,0.04,"P");
	      gPad->SetLeftMargin(0.12);
	    }
	  else             hTrigCuts[j][i]->DrawCopy("Psames");
	}
      if(i<nDays/2+1)
	{
	  if(i==0) leg->AddEntry(hTrigCuts[0][i],"Combined","PL");
	  else     leg->AddEntry(hTrigCuts[0][i],Form("Day %s",days[i]),"PL");
	}
      else
	leg2->AddEntry(hTrigCuts[0][i],Form("Day %s",days[i]),"PL");
    }
  c->cd();
  leg->Draw();
  leg2->Draw();
  if(saveHisto) c->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/AuAu200_CompareDays_TrigTimeWinCut.png"));

  char *type[2] = {"lower","higher"};
  // make ratios
  for(Int_t j=0; j<2; j++)
    {
      for(Int_t i=1; i<nDays; i++)
	{
	  TH1F *hRatio = (TH1F*)hTrigCuts[j][i]->Clone(Form("%s_ratio",hTrigCuts[j][i]->GetName()));
	  hRatio->Divide(hTrigCuts[j][0]);
	  hRatio->SetMarkerStyle(21+j*4);
	  hRatio->SetMarkerColor(color[nDays-i-1]);
	  hRatio->SetLineColor(color[nDays-i-1]);
	  hRatio->GetYaxis()->SetRangeUser(0.998,1.005);
	  if(i==1) 
	    {
	      c = draw1D(hRatio,Form("Ratio of %s-bound trigger window cuts for MTD hits;backleg*5+tray;Ratio",type[j]),kFALSE,kTRUE,0.04,"P");
	      gPad->SetLeftMargin(0.12);
	      leg->Draw();
	      leg2->Draw();
	    }
	  else            
	    hRatio->DrawCopy("Psames");
	  
	}
      if(saveHisto) c->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/AuAu200_CompareDays_%s_TrigTimeWinCut.png",type[j]));
    }
 
}
