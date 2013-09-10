{
//=========Macro generated from canvas: CrossSectionCanvas/CrossSection
//=========  (Wed Jan 19 12:25:07 2011) by ROOT version5.25/04
   TCanvas *CrossSectionCanvas = new TCanvas("CrossSectionCanvas", "CrossSection",10,32,700,500);
   gStyle->SetOptStat(0);
   CrossSectionCanvas->Range(-0.125,-5.480796,1.125,6.360405);
   CrossSectionCanvas->SetBorderSize(2);
   CrossSectionCanvas->SetLogy();
   CrossSectionCanvas->SetFrameFillColor(0);
   
   TF1 *3jet = new TF1("3jet","(2*x*x)/((1-x)*(1-x))",0,1);
   3jet->SetFillColor(19);
   3jet->SetFillStyle(0);
   3jet->SetLineWidth(3);
   3jet->Draw("");
   
   TPaveText *pt = new TPaveText(0.01,0.9390678,0.3029885,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(2);
   pt->SetFillColor(19);
   TText *text = pt->AddText("(2*x*x)/((1-x)*(1-x))");
   pt->Draw();
   CrossSectionCanvas->Modified();
   CrossSectionCanvas->cd();
   CrossSectionCanvas->SetSelected(CrossSectionCanvas);
}
