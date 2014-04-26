#include "Output/mtdQAData.h"
#include "map.h"

const char *triggerName[3] = {"di-muon","single-muon","e-mu"};
const char *mtdHitName[3] = {"NoCut","LooseCut","TightCut"};
const char *tacDiffTitle[9] = {"raw","p_{T,mth} > 1 GeV/c","p_{T,mth} > 1.5 GeV/c","p_{T,mth} > 2 GeV/c","p_{T,mth} > 2 GeV/c && 0<n#sigma_{#pi}<3","p_{T,mth} > 2 GeV/c && 0<n#sigma_{#pi}<3 && |#Deltaz| < 20 cm","p_{T,mth} > 1 GeV/c && 0<n#sigma_{#pi}<3 && |#Deltaz| < 20 cm","p_{T,mth} > 1.5 GeV/c && 0<n#sigma_{#pi}<3 && |#Deltaz| < 20 cm","p_{T} > 3 GeV/c && 0<n#sigma_{#pi}<3 && |#Deltaz| < 20 cm"};
const char *tacDiffName[9] = {"raw","1GeV","1.5GeV","2GeV","2GeV_nSigmaPi0-3","2GeV_nSigmaPi0-3_dZ20","1GeV_nSigmaPi0-3_dZ20","1.5GeV_nSigmaPi0-3_dZ20","3GeV_nSigmaPi0-3_dZ20"};
const char *adcName[4] = {"All","Signal","Pedestal","Zero"};

//================================================
void mtdVpdTac()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.15);                
  gStyle->SetStatH(0.13); 

  //MTDhits();
  MtdVpdTacDiff();
  //CompareDays();
}
  
//================================================
void CompareDays(const Int_t save = 0)
{
  const Int_t nDay = 3;
  const char *days[nDay] = {"077","088","089"};

  TFile *f[nDay];
  TH1F *hTacDiff[9][nDay];
  for(Int_t i=0; i<nDay; i++)
    {
      f[i] = TFile::Open(Form("Rootfiles/output.%s.root",days[i]),"read");
      for(Int_t j=0; j<9; j++)
	{
	  hTacDiff[j][i] = (TH1F*)f[i]->Get(Form("hTacDiff_%s_all",tacDiffName[j]));
	  hTacDiff[j][i]->SetName(Form("%s_%s",hTacDiff[j][i]->GetName(),days[i]));
	}
    }

  TList *list = new TList;
  TString legName[nDay];
  for(Int_t j=0; j<9; j++)
    {
      list->Clear();
      for(Int_t i=0; i<nDay; i++)
	{
	  hTacDiff[j][i]->Rebin(4);
	  //hTacDiff[j][i]->Scale(1./hTacDiff[j][i]->GetEntries());
	  if(j<4)
	    hTacDiff[j][i]->SetMaximum(10*hTacDiff[j][i]->GetMaximum());
	  list->Add(hTacDiff[j][i]);
	  legName[i] = Form("Day %s",days[i]);
	}
      c = drawHistos(list,Form("hTacDiff_%s",tacDiffName[j]),Form("MTD-VPD tac difference (%s);tac_{MTD}-tac_{VPD}+pos.corr.;",tacDiffTitle[j]),kTRUE,-2400,-1400,kFALSE,0.01,1.2e5,kTRUE,kTRUE,legName,kTRUE,"Days",0.2,0.4,0.6,0.85,kFALSE);
      if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/Compare_MtdVpdTacDiff_%s.png",tacDiffName[j]));
    }
}


//================================================
void MTDhits(const char *day = "077", const Int_t save = 0)
{
  TFile *f = TFile::Open(Form("Rootfiles/output.%s.root",day),"read");
  TList *list = new TList;
  TString legName3[3];
  TH1F *hVertexZ[3];
  for(Int_t i=0; i<3; i++)
    {
      hVertexZ[i] = (TH1F*)f->Get(Form("hVertexZ_%s",triggerName[i]));
      list->Add(hVertexZ[i]);
      legName3[i] = triggerName[i];
    }
  c = drawHistos(list,"VertexZ","",kFALSE,2,4.5,kTRUE,0.01,1.2e5,kFALSE,kTRUE,legName3,kTRUE,"Event trigger",0.2,0.4,0.6,0.85,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/VertexZ_%s.png",day,day));
  list->Clear();
  
  TH1F *h1 = (TH1F*)f->Get("hVtxZDiff");
  if(h1)
    {
      c = draw1D(h1);
      if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/VtxZDiffTrkVpd_%s.png",day,day));  
    }

  TH1F *hNTracks[3];
  for(Int_t i=0; i<3; i++)
    {
      hNTracks[i] = (TH1F*)f->Get(Form("hNTracks_%s",triggerName[i]));
      hNTracks[i]->SetLineWidth(2);
      list->Add(hNTracks[i]);
    }
  c = drawHistos(list,"NTracks","",kTRUE,0,80,kFALSE,0.01,1.2e5,kTRUE,kTRUE,legName3,kTRUE,"Event trigger",0.6,0.8,0.6,0.85,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/NGoodGlobalTracks_%s.png",day,day));

  c = drawHistos(list,"NTracks_zoom","# of global tracks per event (zoomed-in)",kTRUE,0,20,kFALSE,0.01,1.2e5,kTRUE,kTRUE,legName3,kTRUE,"Event trigger",0.6,0.8,0.6,0.85,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/NGoodGlobalTracks_zoomed_%s.png",day,day));
  list->Clear();


  TH2F *hMtdHitTrigTimeWest = (TH2F*)f->Get("hMtdHitTrigTimeWest");
  c = draw2D(hMtdHitTrigTimeWest);
  TLine *line = GetLine(0,2750,1800,2750);
  line->Draw();
  line = GetLine(0,3000,1800,3000);
  line->Draw();
  line = GetLine(0,2820,900,2820,1);
  line->Draw();
  line = GetLine(0,2900,900,2900,1);
  line->Draw();
  line = GetLine(900,2860,1800,2860,1);
  line->Draw();
  line = GetLine(900,2950,1800,2950,1);
  line->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/MtdHitTrigTimeWest_%s.png",day,day));

  return;

  TH1F *hTrigTimeBL5 = (TH1F*)hMtdHitTrigTimeWest->ProjectionY("hTrigTimeBL5",241,300);
  TH1F *hTrigTimeBL1 = (TH1F*)hMtdHitTrigTimeWest->ProjectionY("hTrigTimeBL1",1,240);
  hTrigTimeBL1->Scale(0.25);
  list->Add(hTrigTimeBL5);
  list->Add(hTrigTimeBL1);
  TString legName2[2] = {"BL5","BL1-4"};
  c = drawHistos(list,"hTrigTimeBL5","Trigger time of MTD hits",kTRUE,0,80,kFALSE,0.01,1.2e5,kTRUE,kTRUE,legName2,kTRUE,"Trigger time",0.6,0.8,0.6,0.85,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/TrigTimeBL5_%s.png",day,day));
  list->Clear();

  TH2F *h2 = (TH2F*)f->Get("hMtdRawHitMap");
  c = draw2D(h2,"",0.04,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/MtdRawHitMap_%s.png",day,day));

  TH1F *hNMtdHit[3][3];
  TH1F *hNMthMtdHit[3][3];
  TH2F *hMtdHitMap[3];
  TH2F *hMtdHitMapOdd[3];
  for(Int_t j=0; j<3; j++)
    {
      for(Int_t i=0; i<3; i++)
	{
	  hNMtdHit[i][j] = (TH1F*)f->Get(Form("hNMtdHit_%s_%s",triggerName[i],mtdHitName[j]));
	  list->Add(hNMtdHit[i][j]);
	  legName3[i] = Form("%s: mean = %2.4f",triggerName[i],hNMtdHit[i][j]->GetMean());
	}
      c = drawHistos(list,list->At(0)->GetName(),Form("Number of MTD hits per event (%s)",mtdHitName[j]),kFALSE,2,4.5,kFALSE,-0.01,0.03,kTRUE,kTRUE,legName3,kTRUE,"Event trigger",0.5,0.6,0.6,0.85);
      if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/NMtdHit_%s_%s.png",day,mtdHitName[j],day));
      list->Clear();

      for(Int_t i=0; i<3; i++)
	{
	  hNMthMtdHit[i][j] = (TH1F*)f->Get(Form("hNMthMtdHit_%s_%s",triggerName[i],mtdHitName[j]));
	  list->Add(hNMthMtdHit[i][j]);
	  legName3[i] = Form("%s: mean = %2.4f",triggerName[i],hNMthMtdHit[i][j]->GetMean());
	}
      c = drawHistos(list,list->At(0)->GetName(),Form("Number of matched MTD hits per event (%s)",mtdHitName[j]),kFALSE,2,4.5,kFALSE,-0.01,0.03,kTRUE,kTRUE,legName3,kTRUE,"Event trigger",0.5,0.6,0.6,0.85);
      if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/NMthMtdHit_%s_%s.png",day,mtdHitName[j],day));
      list->Clear();

      hMtdHitMap[j] = (TH2F*)f->Get(Form("hMtdHitMap_%s_%s",triggerName[0],mtdHitName[j]));
      hMtdHitMap[j]->SetTitle(Form("%s (%s)",hMtdHitMap[j]->GetTitle(),mtdHitName[j]));
      c = draw2D(hMtdHitMap[j],"",0.04,kFALSE);
      if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/MtdHitMap_%s_%s.png",day,mtdHitName[j],day));

      hMtdHitMapOdd[j] = (TH2F*)f->Get(Form("hMtdHitMapOdd_%s",mtdHitName[j]));
      hMtdHitMapOdd[j]->SetTitle(Form("%s (%s,N_{hits}>30)",hMtdHitMapOdd[j]->GetTitle(),mtdHitName[j]));
      c = draw2D(hMtdHitMapOdd[j],"",0.04,kFALSE);
      if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/MtdHitMapOdd_%s_%s.png",day,mtdHitName[j],day));
    }

  TH1F *hNMtdHitRatio[3];
  for(Int_t i=0; i<3; i++)
    {
      hNMtdHitRatio[i] = (TH1F*)hNMtdHit[i][2]->Clone(Form("hNMtdHitRatio_LooseToTight_%s",triggerName[i]));
      hNMtdHitRatio[i]->Divide(hNMtdHit[i][1]);
      list->Add(hNMtdHitRatio[i]);
      legName3[i] = Form("%s",triggerName[i]);
    }
  c = drawHistos(list,list->At(0)->GetName(),Form("Ratio of the number of MTD hits per event (Tight/Loose);N;ratio"),kTRUE,0,40,kTRUE,-0.5,3,kFALSE,kTRUE,legName3,kTRUE,"Event trigger",0.5,0.6,0.6,0.85);
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/NMtdHitRatio_TightToLoose_%s.png",day,day));
  list->Clear();

  TH1F *hNMthMtdHitVtxCut[3];
  for(Int_t i=0; i<3; i++)
    {
      hNMthMtdHitVtxCut[i]  = (TH1F*)f->Get(Form("hNMthMtdHitVtxCut_%s",triggerName[i]));
      list->Add(hNMthMtdHitVtxCut[i]);
      legName3[i] = Form("%s: mean = %2.4f",triggerName[i],hNMthMtdHitVtxCut[i]->GetMean());
    }
  c = drawHistos(list,list->At(0)->GetName(),Form("Number of matched MTD hits per event (|#Deltaz| < 5 cm)"),kTRUE,0,10,kFALSE,-0.01,0.03,kTRUE,kTRUE,legName3,kTRUE,"Event trigger",0.5,0.6,0.6,0.85);
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/NMthMtdHitVtxCut_%s.png",day,day));
  list->Clear();
}

//================================================
void MtdVpdTacDiff(const char *day = "077", const Int_t save = 1)
{
  TFile *f = TFile::Open(Form("Rootfiles/output.%s.root",day),"read");
  TList *list = new TList;
  TString legName3[3];

  TLegend *leg = new TLegend(0.15,0.65,0.5,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);

  // TAC = 0
  TH1F *hTacFraction[3];
  const char *tacFrac[3] = {"raw","GoodMtdHit","MthMtdHit"};
  for(Int_t i=0; i<3; i++)
    {
      hTacFraction[i] = (TH1F*)f->Get(Form("hTacFraction_%s",tacFrac[i]));
      for(Int_t j=0; j<33; j++)
	{
	  if(j==0) hTacFraction[i]->GetXaxis()->SetBinLabel(j+1,"Sum");
	  else     hTacFraction[i]->GetXaxis()->SetBinLabel(j+1,qtlabel2[j-1]);
	}
      hTacFraction[i]->SetMarkerStyle(20);
      hTacFraction[i]->SetMarkerColor(color[i]);
      hTacFraction[i]->SetLineColor(color[i]);
      hTacFraction[i]->GetYaxis()->SetRangeUser(0,0.5);
      hTacFraction[i]->GetXaxis()->LabelsOption("v");
      hTacFraction[i]->GetXaxis()->SetLabelSize(0.04);
      if(i==0) c = draw1D(hTacFraction[i],"Fraction of zero TAC with ADC>10",kFALSE,kTRUE,0.04,"PL");
      else     hTacFraction[i]->Draw("same PL");
      gPad->SetBottomMargin(0.15);
      leg->AddEntry(hTacFraction[i],tacFrac[i],"PL");
    }
  leg->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/ZeroTacFraction_%s.png",day,day));

  const Int_t qt = 1, pos = 4;
  TH1F *hTac = (TH1F*)f->Get(Form("hQTtac__QT%d_pos%d",qt,pos));
  hTac->GetXaxis()->SetRangeUser(0,2100);
  hTac->SetXTitle("TAC");
  c = draw1D(hTac,"",kTRUE,kFALSE);
  TLine *line = GetLine(20,0,20,hTac->GetMaximum());
  line->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/QTtdc_QT%d_pos%d_%s.png",day,qt,pos,day));

  TH1F *hAdc[4];
  for(Int_t i=0; i<4; i++)
    {
      hAdc[i] = (TH1F*)f->Get(Form("hQTadc_%s_QT%d_pos%d",adcName[i],qt,pos));
      list->Add(hAdc[i]);
    }
  TString legName4[4] = { "No cut on TAC", "100 < TAC < 1300", "20 < TAC < 100", "TAC < 20"};  
  c = drawHistos(list,list->At(0)->GetName(),Form("QT %d - pos %d - J%d: ADC distribution;ADC",qt,(pos-1)/2+1,(pos-1)%2+2),kFALSE,0,10,kFALSE,-0.01,0.03,kTRUE,kTRUE,legName4,kTRUE,"TAC cut",0.5,0.65,0.6,0.85,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/QTadc_QT%d_pos%d_%s.png",day,qt,pos,day));
  list->Clear();

  TH1F *hMtdHitTrigTime[4];
  for(Int_t i=0; i<4; i++)
    {
      hMtdHitTrigTime[i] = (TH1F*)f->Get(Form("hMtdHitTrigTime_%s",adcName[i]));
      list->Add(hMtdHitTrigTime[i]);
    }
  legName4[0] = "No cut on TAC";
  legName4[1] = "TAC > 100 && ADC > 100";  
  legName4[2] = "TAC = 0 && ADC > 100";  
  legName4[3] = "TAC = 0 && ADC < 100";  
  c = drawHistos(list,list->At(0)->GetName(),"Trigger time of MTD hits;hptdc-thub (ns)",kFALSE,0,10,kTRUE,10,1e6,kTRUE,kTRUE,legName4,kTRUE,"TAC cut",0.6,0.75,0.6,0.85,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/MtdHitTrigTime_TACcut_%s.png",day,day));
  list->Clear(); 

  // VPD - MTD tac difference
  TH2F *hMtdVpdTacDiff[9];
  char tmpname[256];
  for(Int_t i=0; i<9; i++)
    {
      sprintf(tmpname,"hMtdVpdTacDiff_%s",tacDiffName[i]);
      hMtdVpdTacDiff[i] = (TH2F*)f->Get(tmpname);
      for(Int_t ibin=1; ibin<=32; ibin++)
      	{
      	  hMtdVpdTacDiff[i]->GetXaxis()->SetBinLabel(ibin,qtlabel2[ibin-1]);
      	}
      hMtdVpdTacDiff[i]->GetXaxis()->LabelsOption("v");
      hMtdVpdTacDiff[i]->GetXaxis()->SetLabelSize(0.04);
      if(i<5) c = draw2D(hMtdVpdTacDiff[i],"");
      else    c = draw2D(hMtdVpdTacDiff[i],"",0.03);
      gPad->SetBottomMargin(0.15);
      if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/MtdVpdTacDiff2D_%s_%s.png",tacDiffName[i],day));
    }

  TH1F *hTacDiff[9][33];
  TH1F *hFraction[9];
  const Int_t lowDiff = -1912;
  const Int_t upDiff  = -1608;
  for(Int_t i=0; i<9; i++)
    {
      hFraction[i] = new TH1F(Form("hFraction_%d",i),Form("Fraction of the signals in the trigger window [%d,%d];;fraction",lowDiff,upDiff),33,0,33);
      for(Int_t j=0; j<33; j++)
	{
	  Int_t lowbin = j, hibin = j;
	  char *name = Form("_QT%d_pos%d", (j-1)/8+1, (j-1)%8+1);
	  if(j==0)
	    {
	      lowbin = 1; hibin = 32;
	      name="_all";
	    }

	  hTacDiff[i][j] = (TH1F*)f->Get(Form("hTacDiff_%s%s",tacDiffName[i],name));

	  if(j==0) hFraction[i]->GetXaxis()->SetBinLabel(j+1,"Sum");
	  else     hFraction[i]->GetXaxis()->SetBinLabel(j+1,qtlabel2[j-1]);
	  if(hTacDiff[i][j]->GetEntries()==0)
	    hFraction[i]->SetBinContent(j+1,0);
	  else
	    hFraction[i]->SetBinContent(j+1,hTacDiff[i][j]->Integral(hTacDiff[i][j]->FindFixBin(lowDiff+0.5),hTacDiff[i][j]->FindFixBin(upDiff-0.5))/ hTacDiff[i][j]->GetEntries());
	}
    }

  const char *hTitle[7] = {"Raw","Matched p_{T} > 1 GeV/c    && -1<n#sigma_{#pi}<3","Matched p_{T} > 1.5 GeV/c && -1<n#sigma_{#pi}<3","Matched p_{T} > 2 GeV/c    && -1<n#sigma_{#pi}<3","Matched p_{T} > 2 GeV/c    && 0<n#sigma_{#pi}<3","Matched p_{T} > 2 GeV/c    && 0<n#sigma_{#pi}<3 && |#Deltaz|<20 cm","Matched p_{T} > 3 GeV/c    && 0<n#sigma_{#pi}<3 && |#Deltaz|<20 cm"};
  TString legName6[6];
  TH1F *hTacDiffClone[6];

  for(Int_t i=0; i<6; i++)
    {
      Double_t fraction = hTacDiff[i][0]->Integral(hTacDiff[i][0]->FindFixBin(lowDiff+0.5),hTacDiff[i][0]->FindFixBin(upDiff-0.5))/ hTacDiff[i][0]->GetEntries();
      TH1F *htmp = (TH1F*)hTacDiff[i][0]->Clone(Form("%s_tmp",hTacDiff[i][0]->GetName()));
      htmp->SetTitle(";tac_{MTD}-tac_{VPD}+pos.corr. (chan);counts");
      ScaleHistoTitle(htmp,26,1,20,26,1.1,22);
      htmp->Rebin(4);
      list->Add(htmp);
      //legName6[i] = Form("%s, f=%2.2f%%",hTitle[i],fraction*100);
      legName6[i] = Form("%s",hTitle[i]);

      hTacDiffClone[i] = (TH1F*)hTacDiff[i][0]->Clone(Form("%s_clone",hTacDiff[i][0]->GetName()));
      hTacDiffClone[i]->Rebin(4);
      hTacDiffClone[i]->Scale(1./hTacDiffClone[i]->GetMaximum());
    }
  c = drawHistos(list,"MtdVpdTacDiff","",kTRUE,-2200,-1600,kTRUE,1,1e9,kTRUE,kTRUE,legName6,kTRUE,Form("Match cuts, trigger window [%d,%d]",lowDiff,upDiff),0.15,0.58,0.62,0.95,kFALSE);
  gPad->SetTopMargin(0.03);
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/Compare_MtdVpdTacDiff_%s.png",day,day));
  

  list->Clear();
  for(Int_t i=0; i<4; i++)
    {
      list->Add(hTacDiffClone[i]);
      legName6[i] = Form("%s",hTitle[i]);
    }
  hTacDiffClone[0]->SetTitle("");
  c = drawHistos(list,"MtdVpdTacDiff_scaled","",kTRUE,-2200,-1600,kTRUE,1e-2,10,kTRUE,kTRUE,legName6,kTRUE,"Match cuts",0.15,0.3,0.7,0.99,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/Compare_MtdVpdTacDiff_Scaled_%s.png",day,day));
  list->Clear();

  for(Int_t i=0; i<6; i++)
    {
      for(Int_t j=0; j<33; j++)
	{
	  hTacDiff[i][j]->GetXaxis()->SetRangeUser(-2400,-1400);
	  if(j>0) hTacDiff[i][j]->Rebin(4);
	  if(i==0 && j==0)    hTacDiff[i][j]->SetMaximum(1e6);
	  hTacDiff[i][j]->SetLineColor(j+1);
	  if(j==0) c = draw1D(hTacDiff[i][j],Form("MTD-VPD tac differece (%s)",hTitle[i]),kTRUE,kFALSE);
	  else     hTacDiff[i][j]->Draw("sames HIST");
	}
      TLine *line = GetLine(lowDiff,0,lowDiff,hTacDiff[i][0]->GetMaximum());
      line->Draw();
      line = GetLine(upDiff,0,upDiff,hTacDiff[i][0]->GetMaximum());
      line->Draw();
      if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/MtdVpdTacDiffPerChannel_%s_%s.png",day,tacDiffName[i],day));
    }  

  TLegend *leg = new TLegend(0.15,0.6,0.5,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  for(Int_t i=0; i<6; i++)
    {
      hFraction[i]->SetMarkerStyle(20);
      hFraction[i]->SetMarkerColor(color[i]);
      hFraction[i]->SetLineColor(color[i]);
      hFraction[i]->GetYaxis()->SetRangeUser(0,1.7);
      hFraction[i]->GetXaxis()->LabelsOption("v");
      hFraction[i]->GetXaxis()->SetLabelSize(0.04);
      if(i==0) c = draw1D(hFraction[i],"",kFALSE,kTRUE,0.04,"PL");
      else     hFraction[i]->Draw("same PL");
      gPad->SetBottomMargin(0.15);
      leg->AddEntry(hFraction[i],hTitle[i],"PL");
    }
  leg->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/SignalFractionPerChannel_%s.png",day,day));

  // Determine fraction of signals in the trigger window
  TH1F *hdiff[3][2];
  TString legName2[2];
  Int_t k = -1;
  for(Int_t i=0; i<3; i++)
    {
      list->Clear();
      Double_t ptCut = 1+i*0.5;
      for(Int_t j=0; j<2; j++)
	{
	  if(j==0) k = i + 1;
	  else 
	    {
	      if(i==0) k = 6;
	      if(i==1) k = 7;
	      if(i==2) k = 5;
	    }
	  hdiff[i][j] = (TH1F*)hTacDiff[k][0]->Clone(Form("tacdiff_%s",tacDiffName[k]));
	  hdiff[i][j]->SetMaximum(1.3*hdiff[i][j]->GetMaximum());
	  list->Add(hdiff[i][j]);
	  if(j==0) legName2[j] = Form("p_{T,mth} > %1.1f GeV/c",ptCut);
	  else     legName2[j] = Form("p_{T,mth} > %1.1f GeV/c && 0<n#sigma_{#pi}<3 && |#Deltaz|<20",ptCut);
	}
      c = drawHistos(list,Form("hTacDiff_All_%1.1fGeV",ptCut),"TAC: MTD-VPD with slewing and position corrections;tacDiff",kTRUE,-2000,-1600,kFALSE,0.01,1.2e5,kFALSE,kTRUE,legName2,kTRUE,"Match cuts",0.15,0.3,0.7,0.88,kFALSE);
      if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/FitTacDiff%1.1fGeV_%s.png",day,ptCut,day));
    }

  Int_t lowBound[4] = {-1912,-1888,-1880,-1864};
  printf("trigger window | [%d,-1600] | [%d,-1600] | [%d,-1600] | [%d,-1600]\n",lowBound[0],lowBound[1],lowBound[2],lowBound[3]);
  printf("Entry counting\n");
  for(Int_t i=0; i<3; i++)
    {
      printf("pt>%1.1f GeV/c   ",1+i*0.5);
      for(Int_t j=0; j<4; j++)
	{
	  printf("| %2.2f%%        ",hdiff[i][1]->Integral(hdiff[i][1]->FindFixBin(lowBound[j]+0.5),hdiff[i][1]->FindFixBin(-1600-0.5))/hdiff[i][1]->Integral(1,3000)*100); 
	}
      printf("\n");
    }

  printf("Gaussian fitting\n");
  for(Int_t i=0; i<3; i++)
    {
      printf("pt>%1.1f GeV/c   ",1+i*0.5);
      TH1F *hFit = (TH1F *)hdiff[i][1]->Clone(Form("%s_Fit",hdiff[i][1]->GetName()));
      TF1 *func1 = new TF1(Form("func_%d",i),"gaus(0)",-1900,-1700);
      func1->SetParameters(100,-1820,60);
      draw1D(hFit);
      hFit->Fit(func1,"RQ");

      for(Int_t j=0; j<4; j++)
	{
	  printf("| %2.2f%%        ",func1->Integral(lowBound[j],-1600)/func1->Integral(func1->GetParameter(1)-3*func1->GetParameter(2),-1600)*100); 
	}
      printf("\n");
    }
      
  hnTacDiff->GetAxis(3)->SetRange(2,4);
  hnTacDiff->GetAxis(2)->SetRangeUser(2.01,5);
  TH1F *hDeltaZ = (TH1F*)hnTacDiff->Projection(4);
  hDeltaZ->SetName("MtdTrkDeltaZ_2GeV_nSigmaPi0-3");
  TF1 *func = new TF1("func","gaus(0)+gaus(3)",-150,150);
  func->SetParameters(100,0,10,100,0,50);
  hDeltaZ->Fit(func,"R0Q");
  c = draw1D(hDeltaZ,"#Deltaz between matched tracks and MTD hits (p_{T}>2 GeV/c, 0<n#sigma_{#pi}<3);#Deltaz (cm)");
  TF1 *func1 = new TF1("func1","gaus",-150,150);
  func1->SetParameters(func->GetParameter(3),func->GetParameter(4),func->GetParameter(5));
  func->SetLineColor(2);
  func->Draw("sames");
  func1->SetLineColor(4);
  func1->Draw("sames");
  Double_t cut = 150;
  TPaveText *t1 = GetPaveText(0.2,0.3,0.65,0.8,0.04);
  t1->AddText(Form("#Deltaz ~ [-%1.0f,%1.0f]",cut,cut));
  t1->AddText(Form("signal ~ %d",(Int_t)(func->Integral(-1*cut,cut)/2)-(Int_t)(func1->Integral(-1*cut,cut)/2)));
  t1->AddText(Form("background ~ %d",(Int_t)(func1->Integral(-1*cut,cut)/2)));
  t1->Draw();
  cut = 20;
  t1 = GetPaveText(0.2,0.3,0.4,0.55,0.04);
  t1->AddText(Form("#Deltaz ~ [-%1.0f,%1.0f]",cut,cut));
  t1->AddText(Form("signal ~ %d",(Int_t)(func->Integral(-1*cut,cut)/2)-(Int_t)(func1->Integral(-1*cut,cut)/2)));
  t1->AddText(Form("background ~ %d",(Int_t)(func1->Integral(-1*cut,cut)/2)));
  t1->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/QA/Plots/MtdVpdTac/%s/FitDeltaZ_%s.png",day,day));
}
