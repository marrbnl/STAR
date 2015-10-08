//================================================
void ana_v2(const int savePlot = 1)
{
  // default method
  double pt_val[4] = { 1.01, 2.91, 4.98, 7.16};

  double pt_err[4] = {0,0,0,0};

  double v2_val[3][3][4] = { {{0.05, 0.01, 0.05, 0.08}, {0.09, -0.01, 0.04, 0.1}, {0.09, -0.01, 0.05, 0.09}},
			     {{0.11, 0.02, 0.07, 0.14}, {0.08, 0.02, 0.06, 0.2}, {0.11, 0.02, 0.05, 0.21}},
			     {{0.08, 0.015, 0.06, 0.10}, {0.084, 0.005, 0.052, 0.14}, {0.106, 0.0083, 0.05, 0.13}}};

  double v2_statErr[3][3][4] = { {{0.03, 0.04, 0.04,0.05}, {0.04, 0.04, 0.05, 0.05}, {0.04, 0.05, 0.05, 0.05}},
				 {{0.03, 0.04, 0.04, 0.06}, {0.03, 0.04, 0.04, 0.06}, {0.02, 0.04, 0.05, 0.07}},
				 {{0.021, 0.028, 0.028, 0.038}, {0.024, 0.028, 0.031, 0.038}, {0.018, 0.03, 0.035, 0.04}}};

  TGraphErrors *gJpsiV2[3][3];
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
	{
	  gJpsiV2[i][j] = new TGraphErrors(4, pt_val, v2_val[i][j], pt_err, v2_statErr[i][j]);
	  gJpsiV2[i][j]->SetName(Form("JpsiV2_%d_%d",i,j));
	  gJpsiV2[i][j]->SetMarkerStyle(20+i+j/2*4);
	  gJpsiV2[i][j]->SetMarkerColor(color[2-i]);
	  gJpsiV2[i][j]->SetLineColor(color[2-i]);
	  gJpsiV2[i][j]->SetMarkerSize(1.5);
	}
      offset_x(gJpsiV2[i][0],0.1*i-0.1);
    }
  offset_x(gJpsiV2[1][2],0.1);

  c = drawGraph(gJpsiV2[0][0],";p_{T} (GeV/c);v_{2}",kFALSE,0.04,"APEZ");
  gJpsiV2[0][0]->GetYaxis()->SetRangeUser(-0.1,0.24);
  for(int i=1; i<3; i++)
    {
      gJpsiV2[i][0]->Draw("samesPEZ");
    }
  TLegend* leg1 = new TLegend(0.4, 0.65, 0.84, 0.88);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetHeader("AuAu 0-80%");
  leg1->AddEntry(gJpsiV2[2][0], "Run10+11", "PL");
  leg1->AddEntry(gJpsiV2[0][0], "Run10", "PL");
  leg1->AddEntry(gJpsiV2[1][0], "Run11", "PL");
  leg1->SetTextSize(0.04);
  leg1->Draw("same");
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/v2/JpsiV2_combined.pdf"));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/v2/JpsiV2_combined.png"));
    }

  const char *name[3] = {"Run10","Run11","Run10+11"};
  for(int i=0; i<3; i++)
    {
      gJpsiV2[i][0]->GetYaxis()->SetRangeUser(-0.1,0.24);
      gJpsiV2[i][0]->SetName(Form("JpsiV2_%d_%d_2",i,j));
      c = drawGraph(gJpsiV2[i][0],";p_{T} (GeV/c);v_{2}",kFALSE,0.04,"APEZ");
      gJpsiV2[i][2]->Draw("samesPEZ");
      TLegend* leg1 = new TLegend(0.2, 0.65, 0.64, 0.88);
      leg1->SetFillStyle(0);
      leg1->SetBorderSize(0);
      leg1->SetHeader(Form("AuAu 0-80%% %s",name[i]));
      leg1->AddEntry(gJpsiV2[i][0], "CyrstalBall+pol2 to unlike", "PL");
      leg1->AddEntry(gJpsiV2[i][2], "CyrstalBall+pol2 to unlike-like", "PL");
      leg1->SetTextSize(0.04);
      leg1->Draw("same");

      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/v2/JpsiV2_Compare_%s.pdf",name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/v2/JpsiV2_Compare_%s.png",name[i]));
	}
    }
}
