TFile *f;
const int year = 2013;
TString run_cfg_name;

//================================================
void qa_dTof()
{
  gStyle->SetOptStat(0);
  if(year==2013) run_type = "Run13_pp500";
  else if(year==2014) run_type = "Run14_AuAu200";
  
  embed();
}


//================================================
void embed(int save = 0)
{
  char *run_config = "Embed.fix.";
  
  TString fileName;
  if(year==2013) fileName = Form("Run13.pp500.jpsi.%sroot",run_config);

  f = TFile::Open(Form("./output/%s",fileName.Data()),"read");
  THnSparseF *hnDtof = (THnSparseF*)f->Get("mhMcTofQA_di_mu");

  // dtof vs pt
  TH2F *hTofVsPt = (TH2F*)hnDtof->Projection(0,3);
  hTofVsPt->SetName("Embed_dTof_vs_pt");
  hTofVsPt->SetTitle("Embedding: tof_{mc} - tof_{exp} vs p_{T}");
  hTofVsPt->GetXaxis()->SetRangeUser(0,12);
  hTofVsPt->GetYaxis()->SetRangeUser(-1,1);
  c = draw2D(hTofVsPt);
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_dTof/%sdTof_vs_pt.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_dTof/%sdTof_vs_pt.png",run_type,run_config));
    }

  // dtof vs mctof
  TH2F *hTofVsMcTof = (TH2F*)hnDtof->Projection(0,1);
  hTofVsMcTof->SetName("Embed_dTof_vs_mc_tof");
  hTofVsMcTof->SetTitle("Embedding: tof_{mc} - tof_{exp} vs tof_{mc}");
  hTofVsMcTof->GetXaxis()->SetRangeUser(13,17);
  hTofVsMcTof->GetYaxis()->SetRangeUser(-1,1);
  c = draw2D(hTofVsMcTof);
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_dTof/%sdTof_vs_mctof.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_dTof/%sdTof_vs_mctof.png",run_type,run_config));
    }

  // dtof vs exptof
  TH2F *hTofVsExpTof = (TH2F*)hnDtof->Projection(0,2);
  hTofVsExpTof->SetName("Embed_dTof_vs_exp_tof");
  hTofVsExpTof->SetTitle("Embedding: tof_{mc} - tof_{exp} vs tof_{exp}");
  hTofVsExpTof->GetXaxis()->SetRangeUser(13,17);
  hTofVsExpTof->GetYaxis()->SetRangeUser(-1,1);
  c = draw2D(hTofVsExpTof);
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_dTof/%sdTof_vs_exptof.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_dTof/%sdTof_vs_exptof.png",run_type,run_config));
    }

  // dtof vs module
  TH2F *hTofVsMod = (TH2F*)hnDtof->Projection(0,4);
  hTofVsMod->SetName("Embed_dTof_vs_module");
  hTofVsMod->SetTitle("Embedding: tof_{mc} - tof_{exp} vs backleg");
  hTofVsMod->GetYaxis()->SetRangeUser(-1,1);
  c = draw2D(hTofVsMod);
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_dTof/%sdTof_vs_backleg.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_dTof/%sdTof_vs_backleg.png",run_type,run_config));
    }

  // dtof intervals
  TList *list = new TList;
  TString legName[3];
  TH1F *histo[3][3];
  double low_bounds[3] = {-5, 0.15, -0.16};
  double up_bounds[3] = {5, 0.16, -0.15};
  char *name[3] = {"mc_tof","exp_tof","pt"};
  double scale[3] = {2,2.5,1.5};
  double min[3] = {12,12,0};
  double max[3]= {18,18,12};
  for(int i=0; i<3; i++)
    {
      list->Clear();
      for(int j=0; j<3; j++)
	{
	  hnDtof->GetAxis(0)->SetRangeUser(low_bounds[j]+0.001,up_bounds[j]-0.001);
	  histo[i][j] = (TH1F*)hnDtof->Projection(i+1);
	  histo[i][j]->SetName(Form("%s_dTofBin%d",name[i],j));
	  histo[i][j]->Scale(1./histo[i][j]->Integral());
	  histo[i][j]->SetLineWidth(2);
	  histo[i][j]->SetMaximum(scale[i]*histo[i][j]->GetMaximum());
	  legName[j] = Form("%2.2f #leq #Deltatof #leq %2.2f ns",low_bounds[j],up_bounds[j]);
	  list->Add(histo[i][j]);
	}
      c = drawHistos(list,name[i],Form("Embedding: %s distributions",name[i]),kTRUE,min[i],max[i],kFALSE,0.1,10,kFALSE,kTRUE,legName,kTRUE,"",0.5,0.7,0.6,0.85,kFALSE);
      if(save)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_dTof/%s%s_in_dTofBin.pdf",run_type,run_config,name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_dTof/%s%s_in_dTofBin.png",run_type,run_config,name[i]));
	}
    }
  hnDtof->GetAxis(0)->SetRange(0,-1);

  // dtof in module
  TH1F *hTofInMod[5];
  TH2F *hTofVsProjMod[5];
  for(int i=0; i<5; i++)
    {
      hnDtof->GetAxis(5)->SetRange(i+1,i+1);
      hTofInMod[i] = (TH1F*)hnDtof->Projection(0);
      hTofInMod[i]->SetName(Form("hTof_Mod%d",i+1));
      hTofVsProjMod[i] = (TH2F*)hnDtof->Projection(0,6);
      hTofVsProjMod[i]->SetName(Form("hTofVsProjMod_Mod%d",i+1));
    }
  hnDtof->GetAxis(5)->SetRange(0,-1);
  TCanvas *c = new TCanvas("hTof_in_mod","hTof_in_mod",1100,750);
  c->Divide(3,2);
  for(int i=0; i<5; i++)
    {
      c->cd(i+1);
      gPad->SetLogy();
      hTofInMod[i]->Draw();
      hTofInMod[i]->SetTitle("");
      hTofInMod[i]->GetXaxis()->SetRangeUser(-1,2);
      TPaveText *t1 = GetTitleText(Form("Module %d",i+1),0.05);
      t1->Draw();
    }
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_dTof/%sdTof_in_module.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_dTof/%sdTof_in_module.png",run_type,run_config));
    }

  TCanvas *c = new TCanvas("hTofVsProjMod_in_mod","hTofVsProjMod_in_mod",1100,750);
  c->Divide(3,2);
  for(int i=0; i<5; i++)
    {
      c->cd(i+1);
      gPad->SetLogz();
      hTofVsProjMod[i]->SetTitle(";track module");
      hTofVsProjMod[i]->GetYaxis()->SetRangeUser(-0.5,0.5);
      hTofVsProjMod[i]->Draw("colz");
      TPaveText *t1 = GetTitleText(Form("Hit in module %d",i+1),0.05);
      t1->Draw();

    }
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_dTof/%sdTofVsProjMod_in_module.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_dTof/%sdTofVsProjMod_in_module.png",run_type,run_config));
    }

}
