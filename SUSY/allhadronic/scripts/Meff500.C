{
//=========Macro generated from canvas: CanvasMeff500/
//=========  (Mon Dec 13 16:30:30 2010) by ROOT version5.25/04
   TCanvas *CanvasMeff500 = new TCanvas("CanvasMeff500", "",10,32,700,500);
   CanvasMeff500->Range(0,0,1,1);
   CanvasMeff500->SetBorderSize(2);
   CanvasMeff500->SetFrameFillColor(0);
  
// ------------>Primitives in pad: PadMeff500
   TPad *PadMeff500 = new TPad("PadMeff500", "",0,0,1,1);
   PadMeff500->Draw();
   PadMeff500->cd();
   PadMeff500->Range(-0.125,-0.125,1.125,1.125);
   PadMeff500->SetBorderSize(2);
   PadMeff500->SetLogz();
   PadMeff500->SetGridx();
   PadMeff500->SetGridy();
   PadMeff500->SetFrameFillColor(0);
   
   TH1F *hframe__5 = new TH1F("hframe__5","M_{eff}=500 GeV, p_{T1}=100 GeV, p_{T2}=100 GeV, p_{T3}=50 GeV",1000,0,1);
   hframe__5->SetMinimum(0);
   hframe__5->SetMaximum(1);
   hframe__5->SetDirectory(0);
   hframe__5->SetStats(0);
   hframe__5->GetXaxis()->SetTitle("x_{2}");
   hframe__5->GetYaxis()->SetTitle("x_{1}");
   hframe__5->Draw(" ");
   
   TH2D *HistoMeff500 = new TH2D("HistoMeff500","",200,0,1,200,0,1);
   HistoMeff500->SetMinimum(1);
   HistoMeff500->SetMaximum(1e+09);
   HistoMeff500->SetEntries(40000);
   HistoMeff500->SetContour(20);
   HistoMeff500->SetContourLevel(0,1);
   HistoMeff500->SetContourLevel(1,2.818383);
   HistoMeff500->SetContourLevel(2,7.943282);
   HistoMeff500->SetContourLevel(3,22.38721);
   HistoMeff500->SetContourLevel(4,63.09573);
   HistoMeff500->SetContourLevel(5,177.8279);
   HistoMeff500->SetContourLevel(6,501.1872);
   HistoMeff500->SetContourLevel(7,1412.538);
   HistoMeff500->SetContourLevel(8,3981.072);
   HistoMeff500->SetContourLevel(9,11220.18);
   HistoMeff500->SetContourLevel(10,31622.78);
   HistoMeff500->SetContourLevel(11,89125.09);
   HistoMeff500->SetContourLevel(12,251188.6);
   HistoMeff500->SetContourLevel(13,707945.8);
   HistoMeff500->SetContourLevel(14,1995262);
   HistoMeff500->SetContourLevel(15,5623413);
   HistoMeff500->SetContourLevel(16,1.584893e+07);
   HistoMeff500->SetContourLevel(17,4.466836e+07);
   HistoMeff500->SetContourLevel(18,1.258925e+08);
   HistoMeff500->SetContourLevel(19,3.548134e+08);
   
   TPaletteAxis *palette = new TPaletteAxis(1.00625,0,1.0625,1,HistoMeff500);
palette->SetLabelColor(1);
palette->SetLabelFont(62);
palette->SetLabelOffset(0.005);
palette->SetLabelSize(0.04);
palette->SetTitleOffset(1);
palette->SetTitleSize(0.04);
   palette->SetFillColor(100);
   palette->SetFillStyle(1001);
   HistoMeff500->GetListOfFunctions()->Add(palette,"br");
   HistoMeff500->Draw("COLZsame");
   
   TPaveText *pt = new TPaveText(0.01,0.9212712,0.71,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(2);
   pt->SetFillColor(19);
   TText *text = pt->AddText("M_{eff}=500 GeV, p_{T1}=100 GeV, p_{T2}=100 GeV, p_{T3}=50 GeV");
   pt->Draw();
   PadMeff500->Modified();
   CanvasMeff500->cd();
  
// ------------>Primitives in pad: OverlayMeff500
   OverlayMeff500 = new TPad("OverlayMeff500", "",0,0,1,1);
   OverlayMeff500->Draw();
   OverlayMeff500->cd();
   OverlayMeff500->Range(-0.125,-0.125,1.125,1.125);
   OverlayMeff500->SetFillColor(0);
   OverlayMeff500->SetFillStyle(4000);
   OverlayMeff500->SetBorderSize(2);
   OverlayMeff500->SetFrameFillStyle(4000);
   OverlayMeff500->SetFrameFillStyle(4000);
   
   TH1F *hframe__6 = new TH1F("hframe__6","",1000,0,1);
   hframe__6->SetMinimum(0);
   hframe__6->SetMaximum(1);
   hframe__6->SetDirectory(0);
   hframe__6->SetStats(0);
   hframe__6->Draw(" ");
   OverlayMeff500->Modified();
   CanvasMeff500->cd();
   CanvasMeff500->Modified();
   CanvasMeff500->cd();
   CanvasMeff500->SetSelected(CanvasMeff500);
}
