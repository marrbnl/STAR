{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  const char *dir = "~/Work/ALICE/Utility";
  gROOT->LoadMacro(Form("%s/definitions.C",dir));
  gROOT->LoadMacro(Form("%s/drawHistos.C",dir));
  gROOT->LoadMacro(Form("%s/sysCompare.C",dir));
  gROOT->LoadMacro(Form("%s/Utility.C",dir));
  gROOT->LoadMacro(Form("%s/HistUtility.C",dir));
  gROOT->LoadMacro(Form("%s/projectToPtSlice.C",dir));
} 
