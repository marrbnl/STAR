void fit_Upsilon()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  /*
  TFile *fin = TFile::Open("output/Pico.Run14.AuAu200.jpsi.mix.pt1.5.pt1.0.root","read");
  TH3D *hULMmumuvsPtCen = (TH3D*)fin->Get("hULMmumuvsPtCen");
  TH3D *hLPosMmumuvsPtCen = (TH3D*)fin->Get("hLPosMmumuvsPtCen");
  TH3D *hLNegMmumuvsPtCen = (TH3D*)fin->Get("hLNegMmumuvsPtCen");
  
  TH1F *hULMmumu = (TH1F*)hULMmumuvsPtCen->ProjectionZ("hULMmumu");
  TH1F *hLSMmumu = (TH1F*)hLPosMmumuvsPtCen->ProjectionZ("hLPosMmumu");
  TH1F *hLNegMmumu = (TH1F*)hLNegMmumuvsPtCen->ProjectionZ("hLNegMmumu");  
  hLSMmumu->Sumw2();
  hLSMmumu->Add(hLNegMmumu);
  
  hULMmumu->Rebin(40);
  hLSMmumu->Rebin(40);
  hULMmumu->GetXaxis()->SetRangeUser(8,11.5);
  hULMmumu->SetMarkerStyle(21);
  hULMmumu->SetMarkerColor(2);
  hULMmumu->SetLineColor(2);
  //hULMmumu->Add(hLSMmumu, -1);
  */

  TFile *fin = TFile::Open("Rootfiles/Ups.root","read");
  TH1F *hULMmumu = (TH1F*)fin->Get("SE_UL");
  TH1F *hLSMmumu = (TH1F*)fin->Get("SE_LS;2");
  TH1F *htmp = (TH1F*)hLSMmumu->Clone("hLSMmumu_fit");
  hULMmumu->Add(hLSMmumu,-1);

  /*
  TF1 *bkg = new TF1("bkg","pol3",8,12);
  c = draw1D(htmp);
  htmp->Fit(bkg,"IR0");

  //return;
  TF1 *fun = new TF1("fun","[0]/sqrt(2*TMath::Pi())/[2]*exp(-pow((x-[1])/sqrt(2)/[2],2))*[3] + [4]*[0]*0.6888361/sqrt(2*TMath::Pi())/[6]*exp(-pow((x-[5])/sqrt(2)/[6],2))*[3] + [4]*[0]*0.3111639/sqrt(2*TMath::Pi())/[8]*exp(-pow((x-[7])/sqrt(2)/[8],2))*[3]+[9]*expo(10)+[12]*pol3(13)",-20,20);
  for(int i=0; i<4; i++) fun->FixParameter(13+i, bkg->GetParameter(i));
  */

  TF1 *fun = new TF1("fun","[0]/sqrt(2*TMath::Pi())/[2]*exp(-pow((x-[1])/sqrt(2)/[2],2))*[3] + [4]*[0]*0.6888361/sqrt(2*TMath::Pi())/[6]*exp(-pow((x-[5])/sqrt(2)/[6],2))*[3] + [4]*[0]*0.3111639/sqrt(2*TMath::Pi())/[8]*exp(-pow((x-[7])/sqrt(2)/[8],2))*[3]+[9]*expo(10)",-20,20);
  fun->SetLineColor(1);
  fun->SetLineWidth(2);

  fun->FixParameter(1,9.46);
  fun->FixParameter(5,10.023);
  fun->FixParameter(7,10.34);

  fun->FixParameter(2,0.18);
  fun->FixParameter(6,0.2);
  fun->FixParameter(8,0.2);

  fun->FixParameter(10, -1.411);
  fun->FixParameter(11, -0.5938);
  fun->SetParName(0,"#varUpsilon(1S) Yield");
  fun->SetParName(4,"#varUpsilon(2S+3S)/#varUpsilon(1S)");
  //fun->SetParName(12,"Scale LS");

  TF1 *r1s = new TF1("r1s","[0]/sqrt(2*TMath::Pi())/[2]*exp(-pow((x-[1])/sqrt(2)/[2],2))*[3]",8,12);
  r1s->SetLineColor(2);
  r1s->SetLineStyle(2);
  r1s->SetLineWidth(2);

  TF1 *r2s = new TF1("r2s","[0]/sqrt(2*TMath::Pi())/[2]*exp(-pow((x-[1])/sqrt(2)/[2],2))*[3]",8,12);
  r2s->SetLineColor(4);
  r2s->SetLineStyle(2);
  r2s->SetLineWidth(2);

  TF1 *r3s = new TF1("r3s","[0]/sqrt(2*TMath::Pi())/[2]*exp(-pow((x-[1])/sqrt(2)/[2],2))*[3]",8,12);
  r3s->SetLineColor(6);
  r3s->SetLineStyle(2);
  r3s->SetLineWidth(2);

  TF1 *dy = new TF1("dy","[0]*expo(1)",8,12);
  dy->SetLineColor(kGreen+2);
  dy->SetLineStyle(2);
  dy->SetLineWidth(2);

  Double_t binWidth = hULMmumu->GetXaxis()->GetBinWidth(1);
  fun->FixParameter(3,binWidth);
  fun->SetParameter(0,30);
  hULMmumu->Fit(fun,"IR0","",8,12);
  r1s->SetParameters(fun->GetParameter(0),fun->GetParameter(1),fun->GetParameter(2),binWidth);
  r2s->SetParameters(fun->GetParameter(0)*fun->GetParameter(4)*0.6888361,fun->GetParameter(5),fun->GetParameter(6),binWidth);
  r3s->SetParameters(fun->GetParameter(0)*fun->GetParameter(4)*0.3111639,fun->GetParameter(7),fun->GetParameter(8),binWidth);
  dy->SetParameters(fun->GetParameter(9),fun->GetParameter(10),fun->GetParameter(11),binWidth);

  hULMmumu->GetYaxis()->SetRangeUser(-40,120);
  c = draw1D(hULMmumu,"Invariant mass distribution of muon pairs");
  //hLSMmumu->SetMarkerStyle(21);
  //hLSMmumu->DrawCopy("sames");
  //bkg->Draw("sames");

  r1s->Draw("sames");
  r2s->Draw("sames");
  r3s->Draw("sames");
  dy->Draw("sames");
  TLegend *leg = new TLegend(0.3,0.6,0.4,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.035);
  leg->AddEntry(hULMmumu,"Unlike-sign","P");
  //leg->AddEntry(hLSMmumu,"Like-sign","P");
  //leg->AddEntry(bkg,"Fit to LS","l");
  leg->AddEntry(r1s,"#varUpsilon(1S)","l");
  leg->AddEntry(r2s,"#varUpsilon(2S)","l");
  leg->AddEntry(r3s,"#varUpsilon(3S)","l");
  leg->AddEntry(dy,"Drell-Yan+bb","l");
  leg->Draw("same");
  fun->SetLineColor(2);
  fun->Draw("sames");

  cout<<"all:"<<fun->Integral(8.4,11.2)/binWidth<<endl;
  cout<<"1s:"<<r1s->Integral(8.4,11.2)/binWidth<<endl;
  cout<<"2s:"<<r2s->Integral(8.4,11.2)/binWidth<<endl;
  cout<<"3s:"<<r3s->Integral(8.4,11.2)/binWidth<<endl;
  
  double n1s = fun->GetParameter(0);
  double e1s = fun->GetParError(0);
  double n23s = fun->GetParameter(4) * fun->GetParameter(0);
  double e23s = fun->GetParError(4) * fun->GetParameter(0);
  double ns = n1s + n23s;
  double es = TMath::Sqrt(e1s*e1s+e23s*e23s);
  TPaveText *t1 = GetPaveText(0.65,0.8,0.4,0.45,0.05);
  t1->AddText(Form("N(#Upsilon) = %1.0f #pm %1.0f",ns,es));
  t1->Draw();

  c->SaveAs("Upsilon1S2S3S.png");
}
