{
//=========Macro generated from canvas: hDalitzMTj1MTjj_IC5Calo_Zinv/
//=========  (Thu Apr 15 10:08:54 2010) by ROOT version5.26/00
   TCanvas *hDalitzMTj1MTjj_IC5Calo_Zinv = new TCanvas("hDalitzMTj1MTjj_IC5Calo_Zinv", "",10,54,700,500);
   hDalitzMTj1MTjj_IC5Calo_Zinv->Range(-130.7288,-120.4696,1176.559,1084.227);
   hDalitzMTj1MTjj_IC5Calo_Zinv->SetBorderSize(2);
   hDalitzMTj1MTjj_IC5Calo_Zinv->SetRightMargin(0.1350575);
   hDalitzMTj1MTjj_IC5Calo_Zinv->SetTopMargin(0.06991526);
   hDalitzMTj1MTjj_IC5Calo_Zinv->SetFrameFillColor(0);
   
   
   TPaveStats *ptstats = new TPaveStats(0.6566092,0.6779661,0.8563218,0.9173729,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(2);
   ptstats->SetFillColor(19);
   ptstats->SetTextAlign(12);
   TText *text = ptstats->AddText("hDalitzMTj1MTjj");
   text->SetTextSize(0.03670904);
   text = ptstats->AddText("Entries = 508    ");
   text = ptstats->AddText("Mean x =  208.1");
   text = ptstats->AddText("Mean y =  562.5");
   text = ptstats->AddText("RMS x =  73.25");
   text = ptstats->AddText("RMS y =  154.7");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   hDalitzMTj1MTjj->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(hDalitzMTj1MTjj->GetListOfFunctions());
   
   TPaletteAxis *palette = new TPaletteAxis(1006.536,0,1065.364,1000,hDalitzMTj1MTjj);
palette->SetLabelColor(1);
palette->SetLabelFont(62);
palette->SetLabelOffset(0.005);
palette->SetLabelSize(0.04);
palette->SetTitleOffset(1);
palette->SetTitleSize(0.04);
   palette->SetFillColor(100);
   palette->SetFillStyle(1001);
   hDalitzMTj1MTjj->GetListOfFunctions()->Add(palette,"br");
   hDalitzMTj1MTjj->GetXaxis()->SetTitle("M_{T}^{2}(j_{1},#slash{E}_{T})");
   hDalitzMTj1MTjj->GetYaxis()->SetTitle("M^{2}(j_{1},j_{2})");
   hDalitzMTj1MTjj->Draw("COLZ ");
   
   TPaveText *pt = new TPaveText(0.01,0.9314407,0.1449425,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(2);
   pt->SetFillColor(19);
   text = pt->AddText("100 pb^{-1}");
   pt->Draw();
   hDalitzMTj1MTjj_IC5Calo_Zinv->Modified();
   hDalitzMTj1MTjj_IC5Calo_Zinv->cd();
   hDalitzMTj1MTjj_IC5Calo_Zinv->SetSelected(hDalitzMTj1MTjj_IC5Calo_Zinv);
}
