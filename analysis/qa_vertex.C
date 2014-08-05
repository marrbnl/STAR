TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;
Int_t vtx_index = 0;
const char *run_config = "PrimTrk.ClosePrimVtx";

//================================================
void qa_vertex(const Int_t save = 1)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.15);                
  gStyle->SetStatH(0.2); 

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  if(cut_name.Contains("ClosePrimVtx"))
    vtx_index = 1;
  else if (cut_name.Contains("MtdVtx"))
    vtx_index = 2;

  TPaveText *vtx_text = 0x0;
  if(vtx_index==0) vtx_text = GetPaveText(0.13,0.4,0.82,0.87);
  else             vtx_text = GetPaveText(0.2,0.6,0.82,0.87);
  //vtx_text->SetFillStyle(1);
  vtx_text->AddText(vtx_name[vtx_index]);
  vtx_text->SetTextColor(kMagenta+1);

  f = TFile::Open(Form("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.%s.root",run_config),"read");
  TH1F *hCuts = (TH1F*)f->Get("hAnalysisCuts");
  Double_t scale = hCuts->GetBinContent(3)/1e4;
  cout << "scale = " << scale << endl;

  // VPD vz
  TH1F *hVzVpd = (TH1F*)f->Get(Form("hVzVpd_%s",trigName[kTrigType]));
  c = draw1D(hVzVpd,Form("Au+Au %s: vertex z distribution reconstructed using VPD%s;v_{z} (cm)",trigName[kTrigType],hlt_name[hlt_index]),kFALSE,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_vertex/%s.VzVpd_%s.png",run_config,trigName[kTrigType]));

  // TPC vz
  TH1F *hVzTpc = (TH1F*)f->Get(Form("hVzTpc_%s",trigName[kTrigType]));
  hVzTpc->SetMaximum(1.2*hVzTpc->GetMaximum());
  c = draw1D(hVzTpc,Form("Au+Au %s: vertex z distribution reconstructed using TPC%s;v_{z} (cm)",trigName[kTrigType],hlt_name[hlt_index]),kFALSE,kFALSE);
  vtx_text->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_vertex/%s.VzTpc_Type%d_%s.png",run_config,vtx_index,trigName[kTrigType]));

  TH1F *hVzTpc_Fit = (TH1F*)hVzTpc->Clone(Form("%s_fit",hVzTpc->GetName()));
  TF1 *func = new TF1("fit_VzTpc","gaus",-30,30);
  func->SetParameters(1,0,10);
  hVzTpc_Fit->Fit(func,"0IR");
  c = draw1D(hVzTpc_Fit,Form("Au+Au %s: vertex z distribution reconstructed using TPC%s;v_{z} (cm)",trigName[kTrigType],hlt_name[hlt_index]),kFALSE,kFALSE);
  vtx_text->Draw();
  func->SetLineColor(2);
  func->Draw("sames");
  Double_t constant = func->GetParameter(0);
  Double_t sigma = func->GetParameter(2);
  TLine *line = GetLine(-40,0,-40,0.5*constant);
  line->Draw();
  line = GetLine(40,0,40,0.5*constant);
  line->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_vertex/%s.Fit_VzTpc_Type%d_%s.png",run_config,vtx_index,trigName[kTrigType]));

  // VPD vs TPC vz
  TH2F *hVzVpdVsTpc = (TH2F*)f->Get(Form("hVzVpdVsTpc_%s",trigName[kTrigType]));
  c = draw2D(hVzVpdVsTpc,Form("Au+Au %s: vertex z VPD vs TPC%s",trigName[kTrigType],hlt_name[hlt_index]));
  vtx_text->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_vertex/%s.VzVpdVsTpc_Type%d_%s.png",run_config,vtx_index,trigName[kTrigType]));

  // TPC-VPD vs TPC vz
  TH2F *hDvzVsTpc = (TH2F*)f->Get(Form("hDvzVsTpc_%s",trigName[kTrigType]));
  c = draw2D(hDvzVsTpc,Form("Au+Au %s: vertex z TPC-VPD vs TPC%s",trigName[kTrigType],hlt_name[hlt_index]));
  vtx_text->Draw();
  Double_t dzCut = hCuts->GetBinContent(11)/scale;
  Double_t low_x = hDvzVsTpc->GetXaxis()->GetBinLowEdge(1);
  Double_t up_x  = hDvzVsTpc->GetXaxis()->GetBinUpEdge(hDvzVsTpc->GetNbinsX());
  TLine *line = GetLine(low_x,-1*dzCut,up_x,-1*dzCut);
  line->Draw();
  line = GetLine(low_x,dzCut,up_x,dzCut);
  line->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_vertex/%s.DvzVsTpc_Type%d_%s.png",run_config,vtx_index,trigName[kTrigType]));

  // TPC vz with cuts
  TH1F *hVzTpcWithCut = (TH1F*)f->Get(Form("hVzTpcWithCut_%s",trigName[kTrigType]));
  hVzTpcWithCut->SetMaximum(1.2*hVzTpcWithCut->GetMaximum());
  c = draw1D(hVzTpcWithCut,Form("Au+Au %s: vertex z distribution reconstructed using TPC%s;v_{z} (cm)",trigName[kTrigType],hlt_name[hlt_index]),kFALSE,kFALSE);
  vtx_text->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_vertex/%s.VzTpcWithCut_Type%d_%s.png",run_config,vtx_index,trigName[kTrigType]));

  // TPC vy vs vx with cuts
  TH2F *hVxVyTpcWithCut = (TH2F*)f->Get(Form("hVxVyTpcWithCut_%s",trigName[kTrigType]));
  c = draw2D(hVxVyTpcWithCut,Form("Au+Au %s: v_{y} vs v_{x} of vertex reconstructed using TPC%s;v_{x} (cm);v_{y} (cm)",trigName[kTrigType],hlt_name[hlt_index]));
  vtx_text->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_vertex/%s.VxVyTpcWithCut_Type%d_%s.png",run_config,vtx_index,trigName[kTrigType]));

  TH1F *hVxTpcWithCut = (TH1F*)hVxVyTpcWithCut->ProjectionX(Form("hVxTpcWithCut_%s",trigName[kTrigType]));
  c = draw1D(hVxTpcWithCut,Form("Au+Au %s: vertex x distribution reconstructed using TPC%s;v_{z} (cm)",trigName[kTrigType],hlt_name[hlt_index]),kFALSE,kFALSE);
  vtx_text->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_vertex/%s.VxTpcWithCut_Type%d_%s.png",run_config,vtx_index,trigName[kTrigType]));

  // TPC vr with cuts
  TH1F *hVrTpcWithCut = (TH1F*)f->Get(Form("hVrTpcWithCut_%s",trigName[kTrigType]));
  hVrTpcWithCut->SetMaximum(1.2*hVrTpcWithCut->GetMaximum());
  c = draw1D(hVrTpcWithCut,Form("Au+Au %s: vertex r distribution reconstructed using TPC%s;v_{r} (cm)",trigName[kTrigType],hlt_name[hlt_index]),kTRUE,kFALSE);
  vtx_text->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_vertex/%s.VrTpcWithCut_Type%d_%s.png",run_config,vtx_index,trigName[kTrigType]));
}
