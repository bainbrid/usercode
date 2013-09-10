{
//=========Macro generated from canvas: CanvasMeff450/
//=========  (Mon Dec 13 16:30:29 2010) by ROOT version5.25/04
   TCanvas *CanvasMeff450 = new TCanvas("CanvasMeff450", "",10,32,700,500);
   CanvasMeff450->Range(0,0,1,1);
   CanvasMeff450->SetBorderSize(2);
   CanvasMeff450->SetFrameFillColor(0);
  
// ------------>Primitives in pad: PadMeff450
   TPad *PadMeff450 = new TPad("PadMeff450", "",0,0,1,1);
   PadMeff450->Draw();
   PadMeff450->cd();
   PadMeff450->Range(-0.125,-0.125,1.125,1.125);
   PadMeff450->SetBorderSize(2);
   PadMeff450->SetLogz();
   PadMeff450->SetGridx();
   PadMeff450->SetGridy();
   PadMeff450->SetFrameFillColor(0);
   
   TH1F *hframe__3 = new TH1F("hframe__3","M_{eff}=450 GeV, p_{T1}=90 GeV, p_{T2}=90 GeV, p_{T3}=45 GeV",1000,0,1);
   hframe__3->SetMinimum(0);
   hframe__3->SetMaximum(1);
   hframe__3->SetDirectory(0);
   hframe__3->SetStats(0);
   hframe__3->GetXaxis()->SetTitle("x_{2}");
   hframe__3->GetYaxis()->SetTitle("x_{1}");
   hframe__3->Draw(" ");
   
   TH2D *HistoMeff450 = new TH2D("HistoMeff450","",200,0,1,200,0,1);
   HistoMeff450->SetMinimum(1);
   HistoMeff450->SetMaximum(1e+09);
   HistoMeff450->SetEntries(40000);
   HistoMeff450->SetContour(20);
   HistoMeff450->SetContourLevel(0,1);
   HistoMeff450->SetContourLevel(1,2.818383);
   HistoMeff450->SetContourLevel(2,7.943282);
   HistoMeff450->SetContourLevel(3,22.38721);
   HistoMeff450->SetContourLevel(4,63.09573);
   HistoMeff450->SetContourLevel(5,177.8279);
   HistoMeff450->SetContourLevel(6,501.1872);
   HistoMeff450->SetContourLevel(7,1412.538);
   HistoMeff450->SetContourLevel(8,3981.072);
   HistoMeff450->SetContourLevel(9,11220.18);
   HistoMeff450->SetContourLevel(10,31622.78);
   HistoMeff450->SetContourLevel(11,89125.09);
   HistoMeff450->SetContourLevel(12,251188.6);
   HistoMeff450->SetContourLevel(13,707945.8);
   HistoMeff450->SetContourLevel(14,1995262);
   HistoMeff450->SetContourLevel(15,5623413);
   HistoMeff450->SetContourLevel(16,1.584893e+07);
   HistoMeff450->SetContourLevel(17,4.466836e+07);
   HistoMeff450->SetContourLevel(18,1.258925e+08);
   HistoMeff450->SetContourLevel(19,3.548134e+08);
   
   TPaletteAxis *palette = new TPaletteAxis(1.00625,0,1.0625,1,HistoMeff450);
palette->SetLabelColor(1);
palette->SetLabelFont(62);
palette->SetLabelOffset(0.005);
palette->SetLabelSize(0.04);
palette->SetTitleOffset(1);
palette->SetTitleSize(0.04);
   palette->SetFillColor(100);
   palette->SetFillStyle(1001);
   HistoMeff450->GetListOfFunctions()->Add(palette,"br");
   HistoMeff450->Draw("COLZsame");
   
   TPaveText *pt = new TPaveText(0.01,0.9212712,0.71,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(2);
   pt->SetFillColor(19);
   TText *text = pt->AddText("M_{eff}=450 GeV, p_{T1}=90 GeV, p_{T2}=90 GeV, p_{T3}=45 GeV");
   pt->Draw();
   PadMeff450->Modified();
   CanvasMeff450->cd();
  
// ------------>Primitives in pad: OverlayMeff450
   OverlayMeff450 = new TPad("OverlayMeff450", "",0,0,1,1);
   OverlayMeff450->Draw();
   OverlayMeff450->cd();
   OverlayMeff450->Range(-0.125,-0.125,1.125,1.125);
   OverlayMeff450->SetFillColor(0);
   OverlayMeff450->SetFillStyle(4000);
   OverlayMeff450->SetBorderSize(2);
   OverlayMeff450->SetFrameFillStyle(4000);
   OverlayMeff450->SetFrameFillStyle(4000);
   
   TH1F *hframe__4 = new TH1F("hframe__4","",1000,0,1);
   hframe__4->SetMinimum(0);
   hframe__4->SetMaximum(1);
   hframe__4->SetDirectory(0);
   hframe__4->SetStats(0);
   hframe__4->Draw(" ");
   OverlayMeff450->Modified();
   CanvasMeff450->cd();
   CanvasMeff450->Modified();
   CanvasMeff450->cd();
   CanvasMeff450->SetSelected(CanvasMeff450);
}
