#include <fstream> 
#include <iostream> 

//const TString year = "Run13";
//const TString year = "Run14";
//const TString year = "Run13_cosmic";
const TString year = "Run14_cosmic";


Int_t nDay;
TString days[100];

//================================================
void trigTimeWinCut()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatW(0.28);                
  gStyle->SetStatH(0.28); 
  gStyle->SetTitleSize(0.06,"T");

  if(year=="Run14")
    {
      nDay = 21;
      char *days_tmp[20] = {"077","088","089","100","104-114","115-120","121","122-129","130-135","136","137-140","141","142-147","148-153","154-162","163","164","165","166","173-186"};
      for(Int_t i=0; i<nDay-1; i++)
	days[i] = days_tmp[i];
    }
  else if (year=="Run13")
    {
      nDay = 35;
      for(Int_t i=0; i<nDay-1; i++)
	{
	  days[i] = Form("%d",i + 130);
	}
      days[nDay-3] = "130-149";
      days[nDay-2] = "151-161";
    }
  else if(year=="Run13_cosmic")
    {
      nDay = 24;
      char *days_tmp[23] = {"93","95","96","97","98","101","102","103","104","106","107","109","112","119","120","125","131","139","141","148","150","156","158"};
      for(Int_t i=0; i<nDay-1; i++)
	days[i] = days_tmp[i];
    }
  else if(year=="Run14_cosmic")
    {
      nDay = 34;
      char *days_tmp[33] = {"50","51","75","79","86","88","89","91","92","94","96","99","100","102","103","106","108","109","115","117","125","128","134","135"
			    ,"137","138","139","154","155","156","157","166","174"};
      for(Int_t i=0; i<nDay-1; i++)
	days[i] = days_tmp[i];
    }

  days[nDay-1] = "All";

  //makeHisto();
  //massProduction();
  processDay();
  //compareDays();
  //example();
}

//================================================
void massProduction()
{
  TFile *fin = TFile::Open(Form("Rootfiles/%s.trigTimeWinCut.root",year.Data()),"read");
  for(Int_t i=0; i<nDay; i++)
    {
      TH2F *h2 = (TH2F*)fin->Get(Form("hMtdHitTrigTime_Day%s",days[i].Data()));
      fit(days[i],h2,1,1,0);
    }
}

//================================================
void processDay(const char * day = "All")
{
  TFile *fin = TFile::Open(Form("Rootfiles/%s.trigTimeWinCut.root",year.Data()),"read");
  TH2F *h2 = (TH2F*)fin->Get(Form("hMtdHitTrigTime_Day%s",day));
  fit(day,h2,0,0,1);
}

//================================================
void fit(const TString day = "all", const TH2F *h2 = 0, const Int_t saveDB = 0, const Int_t saveHisto = 0, const Bool_t debug = 1)
{
  if(!h2) return;

  TH1F *hTrigTime[30][5];
  TF1 *func[30][5];
  TCanvas *c[6];
  TH1F *hCutLow   = new TH1F("trigTimeWin_lowCut","trigTimeWin_lowCut;backleg*5+tray",150,0,150);
  TH1F *hCutHigh  = new TH1F("trigTimeWin_lowHigh","trigTimeWin_lowHigh;backleg*5+tray",150,0,150);
  TH1F *hCutLow2  = new TH1F("trigTimeWin_lowCut_perChan","trigTimeWin_lowCut_perChannel;global channel",1800,0,1800);
  TH1F *hCutHigh2 = new TH1F("trigTimeWin_lowHigh_perChan","trigTimeWin_lowHigh_perChannel;global channel",1800,0,1800);
  TH1F *hMean     = new TH1F("trigTimeWin_mean_draw","Mean of the time distribution of MTH hits;backleg*5+tray",150,0,150);
  TH1F *hSigma    = new TH1F("trigTimeWin_sigma_draw","Sigma of the time distribution of MTH hits;backleg*5+tray",150,0,150);
  TH1F *hMeanE    = new TH1F("trigTimeWin_mean","Mean of the time distribution of MTH hits;backleg*5+tray",150,0,150);
  TH1F *hSigmaE   = new TH1F("trigTimeWin_sigma","Sigma of the time distribution of MTH hits;backleg*5+tray",150,0,150);
  TH1F *hAllHit   = new TH1F("hAllHit","MTD hit multiplicity before trigger window cut;backleg*5+tray",150,0,150);
  TH1F *hGoodHit  = new TH1F("hGoodHit","MTD hit multiplicity after trigger window cut;backleg*5+tray",150,0,150);
  TSpectrum *peakfinder = new TSpectrum();
  Double_t nsigma = 3;
  if(!year.Contains("cosmic") && day=="All") 
    nsigma = 3.5;
  if(year=="Run13" && (day=="130-149"||day=="151-161"))
    nsigma = 3.5;
  
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
	  if(year.Contains("cosmic")) hTrigTime[i][j] = (TH1F*)h2->ProjectionY(Form("MTDhits_triggerTime_BL%d_Tray%d",i+1,j+1),i*5+j+1,i*5+j+1);
	  else hTrigTime[i][j] = (TH1F*)h2->ProjectionY(Form("MTDhits_triggerTime_BL%d_Tray%d",i+1,j+1),i*60+j*12+1,i*60+(j+1)*12);

	  func[i][j] = new TF1(Form("Fit_triggerTime_BL%d_Tray%d",i+1,j+1),"[0]*exp(-0.5*((x-[1])/[2])^2) + [3]");

	  if(hTrigTime[i][j]->GetEntries()<1) continue;
	  TH1F *htmp = (TH1F*)hTrigTime[i][j]->Clone(Form("%s_peak",hTrigTime[i][j]->GetName()));
	  //Int_t nfound = peakfinder->Search(htmp,1,"nodrawnew");
	  Int_t nfound = peakfinder->Search(htmp,2,"nodrawnew",0.5);
	  TList *functions = htmp->GetListOfFunctions();
	  TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
	  Double_t *xpeaks = pm->GetX();
	  if(debug)
	    printf("Found %d candidate peaks, mu = %2.4f\n",nfound,xpeaks[0]);

	  Double_t isigma = 10;
	  if(year.Contains("cosmic")) isigma = 15;

	  Double_t min = xpeaks[0] - 3*isigma;
	  Double_t max = xpeaks[0] + 1*isigma;
	  Double_t imean = xpeaks[0];
	  if(year=="Run14_cosmic" && i==7 && j==0) imean = 2826;
	  if(year=="Run13")
	    {
	      min = imean - 3*isigma;
	      max = imean + 2*isigma;
	    }
	  if(year.Contains("cosmic"))
	    {
	      min = imean - 3*isigma;
	      max = imean + 3*isigma;
	    }
	  
	  hTrigTime[i][j]->GetXaxis()->SetRangeUser(imean - 5*isigma,imean + 10*isigma);
	  hTrigTime[i][j]->SetMaximum(1.5*hTrigTime[i][j]->GetMaximum());
	  hTrigTime[i][j]->SetMinimum(0);
	  hTrigTime[i][j]->SetFillColor(3);
	  hTrigTime[i][j]->SetTitle(Form("BL %d - tray %d",i+1,j+1));
	  c[i/5]->cd(i%5*5+j+1);
	  hTrigTime[i][j]->Draw("");

	  func[i][j]->SetLineColor(2);
	  func[i][j]->SetLineWidth(1);
	  func[i][j]->SetRange(min,max);
	  func[i][j]->SetParameter(1,imean);
	  func[i][j]->SetParameter(2,isigma);
	  func[i][j]->SetParName(0,"C");
	  func[i][j]->SetParName(1,"#mu");
	  func[i][j]->SetParName(2,"#sigma");

	  if(year=="Run14_cosmic" && i==7 && j==0)
	    func[i][j]->FixParameter(2,20);
	  hTrigTime[i][j]->Fit(func[i][j],"IRQ");

	  Double_t mean = func[i][j]->GetParameter(1);
	  Double_t sigma = func[i][j]->GetParameter(2);
	  Double_t yield = func[i][j]->GetParameter(0) + func[i][j]->GetParameter(3);

	  TLine *line = GetLine(mean-nsigma*sigma,0,mean-nsigma*sigma,yield,4);
	  line->Draw("same");
	  line = GetLine(mean+nsigma*sigma,0,mean+nsigma*sigma,yield,4);
	  line->Draw();
	  Int_t bin = i*5+j+1;
	  hCutLow->SetBinContent(bin,mean-nsigma*sigma);
	  hCutHigh->SetBinContent(bin,mean+nsigma*sigma);
	  hMean->SetBinContent(bin,func[i][j]->GetParameter(1));
	  hSigma->SetBinContent(bin,func[i][j]->GetParameter(2));
	  hMeanE->SetBinContent(bin,func[i][j]->GetParameter(1));
	  hMeanE->SetBinError(bin,func[i][j]->GetParError(1));
	  hSigmaE->SetBinContent(bin,func[i][j]->GetParameter(2));
	  hSigmaE->SetBinError(bin,func[i][j]->GetParError(2));
	  hAllHit->SetBinContent(bin,hTrigTime[i][j]->GetEntries());
	  hAllHit->SetBinError(bin,TMath::Sqrt(hTrigTime[i][j]->GetEntries()));
	  Int_t bin1 = hTrigTime[i][j]->FindFixBin(mean-nsigma*sigma);
	  Int_t bin2 = hTrigTime[i][j]->FindFixBin(mean+nsigma*sigma);
	  hGoodHit->SetBinContent(bin,hTrigTime[i][j]->Integral(bin1,bin2));
	  hGoodHit->SetBinError(bin,TMath::Sqrt(hTrigTime[i][j]->Integral(bin1,bin2)));
	  //cout << hTrigTime[i][j]->Integral(bin1,bin2) << " / " << hTrigTime[i][j]->GetEntries() << endl;
	  for(Int_t k = i*60+j*12+1; k<=i*60+(j+1)*12; k++)
	    {
	      hCutLow2->SetBinContent(k,mean-nsigma*sigma);
	      hCutHigh2->SetBinContent(k,mean+nsigma*sigma);
	    }
	}
    }
  c2 = draw2D(h2);
  c2->cd();
  if(!year.Contains("cosmic"))
    {
      hCutLow2->SetLineColor(6);
      hCutLow2->SetLineWidth(2);
      hCutHigh2->SetLineColor(6);
      hCutHigh2->SetLineWidth(2);
      hCutLow2->Draw("sames PL");
      hCutHigh2->Draw("sames PL");
    }
  else
    {
      hCutLow->SetLineColor(6);
      hCutLow->SetLineWidth(2);
      hCutHigh->SetLineColor(6);
      hCutHigh->SetLineWidth(2);
      hCutLow->Draw("sames PL");
      hCutHigh->Draw("sames PL");
    }
  TPaveText *t1 = GetPaveText(0.2,0.6,0.78,0.85);
  t1->SetTextColor(6);
  t1->SetFillStyle(1);
  t1->AddText(Form("Magenta lines indicate %1.1f#sigma cut",nsigma));
  t1->Draw();

  hMean->SetMarkerStyle(20);
  if(year=="Run14")        hMean->GetYaxis()->SetRangeUser(2800,3000);
  if(year=="Run13")        hMean->GetYaxis()->SetRangeUser(2650,2850);
  if(year=="Run14_cosmic") hMean->GetYaxis()->SetRangeUser(2700,3000);
  c11 = draw1D(hMean,"",kFALSE,kTRUE,0.04,"PL");
  hSigma->SetMarkerStyle(20);
  if(year=="Run14")        hSigma->GetYaxis()->SetRangeUser(7,12);
  if(year=="Run13")        hSigma->GetYaxis()->SetRangeUser(7,15);
  if(year=="Run14_cosmic") hSigma->GetYaxis()->SetRangeUser(10,25);
  c12 = draw1D(hSigma,"",kFALSE,kTRUE,0.04,"PL");

  TH1F *hHitRatio = (TH1F*)hGoodHit->Clone("hHitRatio");
  hHitRatio->Divide(hAllHit);
  hHitRatio->SetMarkerStyle(21);
  for(Int_t ibin=1; ibin<=hHitRatio->GetNbinsX(); ibin++)
    {
      hHitRatio->SetBinError(ibin,0);
    }
  hHitRatio->GetYaxis()->SetRangeUser(0,1);
  hHitRatio->SetTitle("Fraction of MTD hits within trigger time window;backleg*5+tray;fraction");
  c13 = draw1D(hHitRatio,"",kFALSE,kTRUE,0.04,"PL");

  cout << day.Data() << endl;
  if(saveHisto)
    {
      c2->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/%s/%s_TrigTimeWinCut.png",year.Data(),day.Data()));
      c11->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/%s/%s_Fit_TrigTime_Mean.png",year.Data(),day.Data()));
      c12->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/%s/%s_Fit_TrigTime_Sigma.png",year.Data(),day.Data()));
      c13->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/%s/%s_GoodMtdHitFraction.png",year.Data(),day.Data()));
      for(Int_t i=0; i<6; i++)
	c[i]->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/%s/%s_Fit_TrigTime_BL%d-%d.png",year.Data(),day.Data(),i*5+1, i*5+5));
    }

  if(debug)
    {
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
    }
  

  if(saveDB)
    {

      //==============================================
      // root file
      TFile *fout = TFile::Open(Form("DB/%s/%s_%s_trigTimeWinCut.root",year.Data(),year.Data(),day.Data()),"recreate");
      hCutLow->Write();
      hCutLow2->Write();
      hCutHigh->Write();
      hCutHigh2->Write();
      hMeanE->Write();
      hSigmaE->Write();
      for(Int_t i=0; i<30; i++)
	{
	  for(Int_t j=0; j<5; j++)
	    {
	      if(hTrigTime[i][j]) hTrigTime[i][j]->Write();
	      if(func[i][j])      func[i][j]->Write();
	    }
	}
      hAllHit->Write();
      hGoodHit->Write();
      hHitRatio->Write();
      fout->Close();

      // PDF file
      char *outPDFName = Form("DB/%s/%s_%s_trigTimeWinFit.pdf",year.Data(),year.Data(),day.Data());
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
      PaintCanvasToPDF(c13,pdf);
      for(Int_t i=0; i<6; i++)
	PaintCanvasToPDF(c[i],pdf);
      PaintCanvasToPDF(c11,pdf);
      PaintCanvasToPDF(c12,pdf);
      pdf->On();
      pdf->Close();
      

      // database
      std::ofstream outfile;
      outfile.open(Form("DB/%s/%s_%s_trigTimeWinCut.dat",year.Data(),year.Data(),day.Data()),std::fstream::out);

      for(Int_t i=0; i<30; i++)
      	{
      	  for(Int_t j=0; j<5; j++)
      	    {
      	      Double_t low = hCutLow->GetBinContent(i*5+j+1);
      	      Double_t hi = hCutHigh->GetBinContent(i*5+j+1);
      	      if(low<1)
      		{
      		  low = 0;
      		  hi = 99999;
      		}

	      if(year=="Run13")
		{
		  if(i+1==7 && j+1==5)
		    {
		      low = 0;
		      hi = 0;
		    }
		}

      	      outfile << i+1 << "  " << j+1 << "  " << "  " << low << "  " << hi << "\n";
      	    }
      	}
      outfile.close();
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
  const Int_t N = nDay;
  TH1F *hTrigCuts[2][N];
  TH1F *hTrigMean[N];
  TH1F *hTrigSigma[N];

  TLegend *leg = new TLegend(0.15,0.7,0.25,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.035);

  TLegend *leg2 = new TLegend(0.35,0.55,0.45,0.88);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextSize(0.035);

  TFile *f = 0x0;
  for(Int_t i=0; i<nDay; i++)
    {
      f = TFile::Open(Form("DB/%s/%s_Day%s_trigTimeWinCut.root",year.Data(),year.Data(),days[i].Data()),"read");
      hTrigCuts[0][i] = (TH1F*)f->Get("trigTimeWin_lowCut");
      hTrigCuts[1][i] = (TH1F*)f->Get("trigTimeWin_lowHigh");
      hTrigMean[i]    = (TH1F*)f->Get("trigTimeWin_mean");
      hTrigSigma[i]   = (TH1F*)f->Get("trigTimeWin_sigma");
      hTrigMean[i]    ->SetName(Form("%s_day_%s",hTrigMean[i]->GetName(),days[i].Data()));
      hTrigSigma[i]   ->SetName(Form("%s_day_%s",hTrigSigma[i]->GetName(),days[i].Data()));

      for(Int_t j=0; j<2; j++)
	{
	  hTrigCuts[j][i]->SetName(Form("%s_day_%s",hTrigCuts[j][i]->GetName(),days[i].Data()));
	  hTrigCuts[j][i]->SetMarkerStyle(25);
	  hTrigCuts[j][i]->SetMarkerColor(color[nDay-i-1]);
	  hTrigCuts[j][i]->SetLineColor(color[nDay-i-1]);
	  if(year=="Run14") hTrigCuts[j][i]->GetYaxis()->SetRangeUser(2800,3000);
	  if(year=="Run13") hTrigCuts[j][i]->GetYaxis()->SetRangeUser(2650,2900);
	  hTrigCuts[j][i]->GetYaxis()->SetNdivisions(505);
	  hTrigCuts[j][i]->GetYaxis()->SetTitleOffset(1.5);
	  if(i==0 && j==0) 
	    {
	      c = draw1D(hTrigCuts[j][i],"Offline trigger time window cuts for MTD hits;backleg*5+tray;TAC_{MTD}-TAC_{THUB}",kFALSE,kTRUE,0.045,"P");
	      gPad->SetLeftMargin(0.12);
	    }
	  else if(i<nDay-1) hTrigCuts[j][i]->DrawCopy("Psames");
	  else 
	    {
	      //hTrigCuts[j][i]->SetLineWidth(2);
	      hTrigCuts[j][i]->SetMarkerStyle(21);
	      hTrigCuts[j][i]->SetMarkerColor(1);
	      hTrigCuts[j][i]->DrawCopy("Psames");
	    }
	}
      // if(i<nDay/2+1)
      // 	leg->AddEntry(hTrigCuts[0][i],Form("Day %s",days[i].Data()),"PL");
      // else
      // 	leg2->AddEntry(hTrigCuts[0][i],Form("Day %s",days[i].Data()),"PL");
      if(i==nDay-2) leg->AddEntry(hTrigCuts[0][i],"Individual days (3#sigma cut)","P");
      if(i==nDay-1) leg->AddEntry(hTrigCuts[0][i],"Combined (3.5#sigma cut)","P");
    }
  c->cd();
  leg->Draw();
  //leg2->Draw();
  if(saveHisto) c->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/%s/CompareDays_TrigTimeWinCut.png",year.Data()));

  TH1F *hMean[150];
  TH1F *hSigma[150];
  TH1F *hLowCut[150];
  TH1F *hHighCut[150];
  TCanvas *cMean[5];
  TCanvas *cSigma[5];
  TCanvas *cCut[30];
  for(Int_t i=0; i<5; i++)
    {
      cMean[i] = new TCanvas(Form("cMean_%d",i),Form("Mean_%d",i),1450,800);
      cMean[i]->Divide(3,2);
      cSigma[i] = new TCanvas(Form("cSigma_%d",i),Form("cSigma_%d",i),1450,800);
      cSigma[i]->Divide(3,2);
    }

  for(Int_t i=0; i<30; i++)
    {
      cCut[i] = new TCanvas(Form("cCut_%d",i+1),Form("cCut_%d",i+1),1450,800);
      cCut[i]->Divide(3,2);
    }

  for(Int_t i=0; i<150; i++)
    {
      Int_t bl = i/5 + 1;
      Int_t tray = i%5 + 1;
      Int_t nbins = nDay-1;
      if(year=="Run13") nbins = nDay - 3;
      hMean[i]  = new TH1F(Form("hMean_BL%d_Mod%d",bl,tray),";;Mean",nbins,0,nbins);
      hSigma[i] = new TH1F(Form("hSigma_BL%d_Mod%d",bl,tray),";;Sigma",nbins,0,nbins);
      hLowCut[i] = new TH1F(Form("hLowCut_BL%d_Mod%d",bl,tray),";;Time cut",nbins,0,nbins);
      hHighCut[i] = new TH1F(Form("hHighCut_BL%d_Mod%d",bl,tray),";;Time cut",nbins,0,nbins);
      for(Int_t j=0; j<nbins; j++)
	{
	  hMean[i]->GetXaxis()->SetBinLabel(j+1,days[j]);
	  hMean[i]->GetXaxis()->LabelsOption("v");
	  hMean[i]->SetBinContent(j+1,hTrigMean[j]->GetBinContent(i+1));
	  hMean[i]->SetBinError(j+1,hTrigMean[j]->GetBinError(i+1));
	  hMean[i]->SetMarkerStyle(21);
	  hMean[i]->SetMarkerColor(color[tray-1]);
	  hMean[i]->SetLineColor(color[tray-1]);
	  hMean[i]->GetYaxis()->SetTitleOffset(1.5);
	  hMean[i]->GetXaxis()->SetLabelSize(0.06);

	  hSigma[i]->GetXaxis()->SetBinLabel(j+1,days[j]);
	  hSigma[i]->GetXaxis()->LabelsOption("v");
	  hSigma[i]->SetBinContent(j+1,hTrigSigma[j]->GetBinContent(i+1));
	  hSigma[i]->SetBinError(j+1,hTrigSigma[j]->GetBinError(i+1));
	  hSigma[i]->SetMarkerStyle(21);
	  hSigma[i]->SetMarkerColor(color[tray-1]);
	  hSigma[i]->SetLineColor(color[tray-1]);
	  hSigma[i]->GetYaxis()->SetTitleOffset(1.5);
	  hSigma[i]->GetXaxis()->SetLabelSize(0.06);

	  hLowCut[i]->GetXaxis()->SetBinLabel(j+1,days[j]);
	  hLowCut[i]->GetXaxis()->LabelsOption("v");
	  hLowCut[i]->SetBinContent(j+1,hTrigCuts[0][j]->GetBinContent(i+1));
	  //hLowCut[i]->SetBinError(j+1,hTrigCuts[0][j]->GetBinError(i+1));
	  hLowCut[i]->SetMarkerStyle(21);
	  hLowCut[i]->GetYaxis()->SetTitleOffset(1.5);
	  hLowCut[i]->GetXaxis()->SetLabelSize(0.06);

	  hHighCut[i]->GetXaxis()->SetBinLabel(j+1,days[j]);
	  hHighCut[i]->GetXaxis()->LabelsOption("v");
	  hHighCut[i]->SetBinContent(j+1,hTrigCuts[1][j]->GetBinContent(i+1));
	  //hHighCut[i]->SetBinError(j+1,hTrigCuts[1][j]->GetBinError(i+1));
	  hHighCut[i]->SetMarkerStyle(21);
	  hHighCut[i]->GetYaxis()->SetTitleOffset(1.5);
	  hHighCut[i]->GetXaxis()->SetLabelSize(0.06);

	  if(year=="Run14")
	    {
	      if(bl<=15) hMean[i]->GetYaxis()->SetRangeUser(2840,2875);
	      else if (bl==24) hMean[i]->GetYaxis()->SetRangeUser(2870,2910);
	      else             hMean[i]->GetYaxis()->SetRangeUser(2890,2930);
	      hSigma[i]->GetYaxis()->SetRangeUser(7,12);
	    }
	  else if(year=="Run13")
	    {
	      hSigma[i]->GetYaxis()->SetRangeUser(6,15);
	      if(bl<=15) hMean[i]->GetYaxis()->SetRangeUser(2720,2770);
	      else             hMean[i]->GetYaxis()->SetRangeUser(2770,2820);
	    }
	}
      cMean[(bl-1)/6]->cd((bl-1)%6+1);
      SetPadMargin(gPad,0.15,0.12,0.05);
      gPad->SetGridy();
     
      if(tray==1)
	{
	  hMean[i]->Draw("P");
	  TPaveText *t = GetTitleText(Form("Backleg %d",bl),0.06);
	  t->Draw();
	}
      else
	hMean[i]->Draw("samesP");

      cSigma[(bl-1)/6]->cd((bl-1)%6+1);
      SetPadMargin(gPad,0.15,0.12,0.05);
      gPad->SetGridy();
     
      if(tray==1)
	{
	  hSigma[i]->Draw("P");
	  TPaveText *t = GetTitleText(Form("Backleg %d",bl),0.06);
	  t->Draw();
	}
      else
	hSigma[i]->Draw("samesP");

      cCut[bl-1]->cd(tray);
      SetPadMargin(gPad,0.15,0.12,0.05);
      gPad->SetGridy();
      if(year.Data()=="Run14")
	{
	  if(bl<=15) hHighCut[i]->GetYaxis()->SetRangeUser(2780,2940);
	  else if (bl==24) hHighCut[i]->GetYaxis()->SetRangeUser(2810,2970);
	  else             hHighCut[i]->GetYaxis()->SetRangeUser(2840,3000);
	}
      else if(year.Data()=="Run13")
	{
	  if(bl<=15) hHighCut[i]->GetYaxis()->SetRangeUser(2660,2830);
	  else       hHighCut[i]->GetYaxis()->SetRangeUser(2710,2880);
	}
      hHighCut[i]->Draw("P");
      hLowCut[i]->Draw("samesP");
      TPaveText *t = GetTitleText(Form("Backleg %d, module %d",bl,tray),0.06);
      t->Draw();
      if(year.Data()=="Run14")
	{
	  TLine *line1 = GetLine(0,hTrigCuts[1][nDay-1]->GetBinContent(i+1),nDay,hTrigCuts[1][nDay-1]->GetBinContent(i+1),2,1,1); 
	  line1->Draw();
	  TLine *line2 = GetLine(0,hTrigCuts[0][nDay-1]->GetBinContent(i+1),nDay,hTrigCuts[0][nDay-1]->GetBinContent(i+1),2,1,1);
	  line2->Draw();
	}
      if(year.Data()=="Run13")
	{
	  TLine *line1 = GetLine(0,hTrigCuts[1][nDay-3]->GetBinContent(i+1),20,hTrigCuts[1][nDay-3]->GetBinContent(i+1),2,1,1);
	  line1->Draw();
	  line1 = GetLine(0,hTrigCuts[0][nDay-3]->GetBinContent(i+1),20,hTrigCuts[0][nDay-3]->GetBinContent(i+1),2,1,1);
	  line1->Draw();

	  TLine *line2 = GetLine(20,hTrigCuts[1][nDay-2]->GetBinContent(i+1),nDay-3,hTrigCuts[1][nDay-2]->GetBinContent(i+1),2,1,1);
	  line2->Draw();
	  TLine *line2 = GetLine(20,hTrigCuts[0][nDay-2]->GetBinContent(i+1),nDay-3,hTrigCuts[0][nDay-2]->GetBinContent(i+1),2,1,1);
	  line2->Draw();
	}
    }

  if(saveHisto)
    {
      for(Int_t i=0; i<5; i++)
	{
	  cMean[i]->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/%s/CompareDays_Mean_BL%d-%d.png",year.Data(),i*6+1,i*6+6));
	  cSigma[i]->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/%s/CompareDays_Sigma_BL%d-%d.png",year.Data(),i*6+1,i*6+6));
	}
      for(Int_t i=0; i<30; i++)
	{
	  cCut[i]->SaveAs(Form("~/Work/STAR/Calibration/Plots/trigTimeWin/%s/CompareDays_TrigTimeCut_BL%d.png",year.Data(),i+1));
	}
    }

  if(saveHisto)
    {
      char *outPDFName = Form("%s_Calib_trigTimeWinCut.pdf",year.Data());
      c->Print(Form("%s(",outPDFName));
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
      for(Int_t i=0; i<5; i++)
	{
	  PrintCanvasToPDF(cMean[i],pdf);
	  PrintCanvasToPDF(cSigma[i],pdf);
	}
      for(Int_t i=0; i<30; i++)
	{
	  PrintCanvasToPDF(cCut[i],pdf);
	}
      pdf->On();
      pdf->Close();
    }
}

//================================================
void makeHisto()
{
  TFile *fout = TFile::Open(Form("Rootfiles/%s.trigTimeWinCut.root",year.Data()),"recreate");

  const char *listName = Form("file.%s.list",year.Data());
  ifstream input;
  input.open(listName);
  
  TString histoName = "hMtdHitTrigTimeWest";
  if(year.Contains("cosmic"))
    histoName = "hTriggerTimeBL";

  char fileName[256];
  while(!input.eof())
    {
      input.getline(fileName,256);
      TString file = fileName;
      TString run = file;
      if(file.Length()>0)
	{
	  run.Remove(run.Index("root")-1,5);
	  run.Remove(0,run.Index("histos")+7);
	  TFile *fin = TFile::Open(file.Data(),"read");
	  TH2F *h2 = (TH2F*)fin->Get(histoName.Data());
	  fout->cd();
	  h2->Write(Form("hMtdHitTrigTime_Day%s",run.Data()));
	  fin->Close();
	}
    }
  fout->cd();
  fout->Close();
  input.close();
}

