

void trigTimeWinCut(const char * day = "121", const Int_t saveDB = 0, const Int_t saveHisto = 0)
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  //gStyle->SetStatY(0.9);                
  //gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.28);                
  gStyle->SetStatH(0.28); 
  gStyle->SetTitleSize(0.08,"T");

  char inFile[256];
  char outFile[256];
  if(day=="all")
    {
      sprintf(inFile,"Output/histos.AuAu200.root");
      sprintf(outFile,"AuAu200");
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
      if(i>4) continue;
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
	  if(i<15) hTrigTime[i][j]->GetXaxis()->SetRangeUser(2800,3100);
	  else     hTrigTime[i][j]->GetXaxis()->SetRangeUser(2850,3100);
	  hTrigTime[i][j]->SetFillColor(3);
	  hTrigTime[i][j]->SetTitle(Form("BL %d - tray %d",i+1,j+1));
	  hTrigTime[i][j]->Draw("");

	  Int_t nfound = peakfinder->Search(hTrigTime[i][j],1,"goff");
	  printf("Found %d candidate peaks\n",nfound);

	  func[i][j] = new TF1(Form("Fit_triggerTime_BL%d_Tray%d",i+1,j+1),"gaus(0)+gaus(3)");
	  func[i][j]->SetLineColor(2);
	  func[i][j]->SetLineWidth(2);
	  if(i<15)
	    {
	      func[i][j]->SetRange(2800,3000);
	      func[i][j]->SetParameter(1,2850);
	      func[i][j]->SetParameter(2,10);
	      func[i][j]->SetParameter(4,2930);
	      func[i][j]->SetParameter(5,15);
	    }
	  else
	    {
	      func[i][j]->SetRange(2870,3070);
	      func[i][j]->SetParameter(1,2920);
	      func[i][j]->SetParameter(2,10);
	      func[i][j]->SetParameter(4,2970);
	      func[i][j]->SetParameter(5,40);
	    }
	  if(i==14 && j==2)
	    {
	      func[i][j]->SetParameter(1,2970);
	      func[i][j]->SetParameter(4,3170);
	    }
	  if(i==24 && j==0)
	    {
	      func[i][j]->SetParameter(1,2870);
	      func[i][j]->SetParameter(4,3170);
	    }
	  if(hTrigTime[i][j]->GetEntries()<1) continue;
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
  hSigma->GetYaxis()->SetRangeUser(7,10);
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
      TFile *fout = TFile::Open(Form("DB/trigTimeWinCut_%s.root",outFile),"recreate");
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
	      hTrigTime[i][j]->Write();
	      func[i][j]->Write();
	    }
	}
      fout->Close();
      
      // PDF file
      char *outPDFName = Form("DB/trigTimeWinCut_%s.pdf",outFile);
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
      out.open(Form("DB/trigTimeWinCut_%s.dat",outFile),std::ofstream::out);
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
