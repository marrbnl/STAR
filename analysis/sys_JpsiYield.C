TFile *f;
const double pt1_cut = 1.5, pt2_cut = 1.0;
const double low_mass = 2.9;
const double high_mass = 3.3;

const char *run_config = "";
const Bool_t iPico = 1;
const int year = 2014;
const int mix_type = 0; // 0 - Shuai; 1 - Rongrong
TString run_cfg_name;
TString file_name;

//================================================
void sys_JpsiYield()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  TString cut_name = run_config;
  if(year==2013)
    {
      run_type = "Run13_pp500";
      if(iPico) file_name = Form("Pico.Run13.pp500.jpsi.%sroot",run_config);
      else      file_name = Form("Run13.pp500.jpsi.%sroot",run_config)
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";

      if(iPico) file_name = Form("Pico.Run14.AuAu200.jpsi.%sroot",run_config);
      else      file_name = Form("Run14.AuAu200.jpsi.%sroot",run_config);
    }
  run_cfg_name = Form("%s",run_config);
  if(iPico) run_cfg_name = Form("Pico.%s",run_cfg_name.Data());

  signalExtraction();
}

//================================================
void signalExtraction(int save = 1)
{
  const int nSys = 7;
  const char *sys_name[nSys] = {"","_LargeScale","_SmallScale","_LargeFit","_SmallFit","_Rebin","_pol1"};
  TString outName = file_name;
  outName.ReplaceAll(".root",".yield.root");

  TString outNameSys = file_name;
  outNameSys.ReplaceAll(".root",".sys.signal.root");

  f = TFile::Open(Form("Rootfiles/%s",outName.Data()),"read");
  TFile *fSys = TFile::Open(Form("Rootfiles/%s",outNameSys.Data()),"read");
  TH1F *hSignal[nCentBins][nSys];

  const double sys_value = 0.06;
  for(int i=0; i<nCentBins; i++)
    {
      TCanvas *c = new TCanvas(Form("Sys_Signal_cent%s",cent_Title[i]),Form("Sys_Signal_cent%s",cent_Title[i]),800,600);
      for(int j=0; j<6; j++)
	{
	  if(j==0) hSignal[i][j] = (TH1F*)f->Get(Form("Jpsi_BinCountYield_cent%s%s",cent_Title[i],sys_name[j]));
	  else     hSignal[i][j] = (TH1F*)fSys->Get(Form("Jpsi_BinCountYield_cent%s%s",cent_Title[i],sys_name[j]));
	  hSignal[i][j]->SetMarkerSize(1.5);
	  hSignal[i][j]->SetMarkerStyle(21);
	  hSignal[i][j]->SetMarkerColor(color[j]);
	  
	  TH1F *htmp = (TH1F*)hSignal[i][j]->Clone(Form("%s_tmp",hSignal[i][j]->GetName()));
	  htmp->Divide(hSignal[i][0]);
	  for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
	    {
	      if(j==0) htmp->SetBinError(bin,hSignal[i][j]->GetBinError(bin)/hSignal[i][j]->GetBinContent(bin));
	      else     htmp->SetBinError(bin,0);
	    }
	  htmp->GetYaxis()->SetRangeUser(0.4,1.7);
	  htmp->SetTitle(";p_{T} (GeV/c);Relative difference");
	  if(j==0) htmp->Draw("P");
	  else htmp->Draw("PL sames");
	  TPaveText *t1 = GetTitleText(Form("Systematic uncertainty of signal extraction (%s%%)",cent_Name[i]));
	  t1->Draw();
	}
      double xmin = 0.15, xmax = 0.4, ymin = 0.7, ymax = 0.88;
      leg1 = new TLegend(xmin,ymin,xmax,ymax);
      leg1->SetBorderSize(0);
      leg1->SetFillColor(0);
      leg1->SetTextSize(0.035);
      leg1->AddEntry(hSignal[i][0],"Default with stat. err.","P");
      leg1->AddEntry(hSignal[i][1],"Larger scale range","P");
      leg1->AddEntry(hSignal[i][2],"Smaller scale range","P");
      leg1->Draw();

      leg2 = new TLegend(xmin+0.3,ymin,xmax+0.3,ymax);
      leg2->SetBorderSize(0);
      leg2->SetFillColor(0);
      leg2->SetTextSize(0.035);
      leg2->AddEntry(hSignal[i][3],"Larger fit range","P");
      leg2->AddEntry(hSignal[i][4],"Smaller fit range","P");
      leg2->AddEntry(hSignal[i][5],"Different binning","P");
      leg2->Draw();

      TLine *line = GetLine(0,1+sys_value,10,1+sys_value);
      line->Draw();
      line = GetLine(0,1-sys_value,10,1-sys_value);
      line->Draw();

      if(save)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_JpsiYield/%sSys_sigExt_cent%s.pdf",run_type,run_cfg_name.Data(),cent_Title[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_JpsiYield/%sSys_sigExt_cent%s.png",run_type,run_cfg_name.Data(),cent_Title[i]));
	}
    }
}
