
void unfold(int test1 = 0, int test2 = -1, int test3 = -1, int itermax = 10000)
{
    int nbinsx = 5010;
    double xmin = 0;
    double xmax = nbinsx;
    
    int nbinsy = 64;
    double ymin = 0;
    double ymax = nbinsy;
    
    double ** response = new double * [nbinsy];
    for(int i=0; i<nbinsy; i++) response[i] = new double[nbinsx];
    TH1F * d = new TH1F("d","Decomposition - unfolding",nbinsy,ymin,ymax);

    TFile * f = TFile::Open("~/data/hists2012.root");
    TH1F ** hists = new TH1F * [nbinsy];
    for(int i=0; i<nbinsy; i++) { 
        hists[i] = (TH1F*)f->Get(Form("ScionixCal%d",i)); 
        //hists[i] = (TH1F*)f->Get(Form("ProtonCal%d",i)); 
        if(hists[i]->GetNbinsX() == 50100) hists[i]->Rebin(10);
    }
    for(int i=0; i<nbinsy; i++) {
        for(int j=0; j<nbinsx; j++) {
            //if( i != 0 ) response[i][j] = hists[i]->GetBinContent(j);
            //else response[i][j] = 0;
            response[i][j] = hists[i]->GetBinContent(j);
        }
    }
    

    double * source = new double[nbinsx];
    TH1F * h = (TH1F*)hists[test1]->Clone();
    if(test2 != -1 ) h->Add(hists[test2]);
    if(test3 != -1 ) h->Add(hists[test3]);
    for(int i=0; i<nbinsx; i++) source[i] = h->GetBinContent(i);
    
    TCanvas * canv = 0;
    canv = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("canv");
    if(!canv) {canv = new TCanvas("canv","canv");
    //TCanvas * canv = new TCanvas();
    canv->Divide(1,2);}
    canv->cd(1);
    
    h->SetLineColor(kBlack);
    h->SetLineWidth(3);
    h->Scale(1.1);
    h->Draw();
    hists[test1]->SetLineColor(kRed);
    hists[test1]->Draw("same");
    if(test2 != -1) { hists[test2]->SetLineColor(kBlue); hists[test2]->Draw("same"); }
    if(test3 != -1) { hists[test3]->SetLineColor(kGreen); hists[test3]->Draw("same"); }

    TSpectrum * spec = new TSpectrum();
    //spec->Unfolding(source,(const double **)response, nbinsx, nbinsy, 1000, 1, 1);
    spec->Unfolding(source,(const double **)response, nbinsx, nbinsy, itermax, 1, 1);
    for(int i=0; i<nbinsy; i++) d->SetBinContent(i,source[i]);
    //d->Draw();
    
    double energy_vector[] =
    {
        0.66201 ,0.80879    ,1.05264    ,1.21767    ,2.96445    ,2.86088    ,2.69813    ,2.48962    ,2.25152    ,
        2.00070 ,1.75268    ,1.52003    ,2.27882    ,2.60832    ,2.95748    ,3.30884    ,3.64108    ,3.93118    ,
        4.15713 ,4.30074    ,1.98127    ,1.72261    ,1.50509    ,1.50509    ,1.72261    ,1.98127    ,0.81657    ,
        0.95882 ,1.11930    ,1.29208    ,1.46825    ,1.63655    ,1.78465    ,1.90062    ,1.97458    ,0.69544    ,
        0.59575 ,0.51611    ,0.45421    ,0.38915    ,0.68765    ,0.65189    ,0.59641    ,0.48874    ,0.37184    ,
        0.26609 ,0.18371    ,0.12684    ,0.10089    ,0.07506    ,0.05996    ,0.05996    ,7.90894    ,7.64353    ,
        7.22603 ,6.69019    ,6.07682    ,5.42853    ,4.78473    ,4.17760    ,3.62976    ,3.15389    ,2.75400    ,
        2.42779
    };
    int run_vector[] =
    {
        50  ,49  ,48  ,47  ,46  ,45  ,44  ,39  ,38  ,43  ,37  ,36  ,42  ,41  ,0   ,40  ,35  ,1   ,26  ,27  ,2   ,
        28  ,3   ,29  ,30  ,22  ,23  ,11  ,31  ,21  ,24  ,10  ,32  ,33  ,34  ,20  ,25  ,9   ,8   ,12  ,63  ,7   ,
        13  ,6   ,62  ,5   ,14  ,4   ,61  ,15  ,60  ,16  ,17  ,18  ,59  ,19  ,58  ,57  ,56  ,55  ,54  ,53  ,52
    };

    TH1F * energy = new TH1F("energy","energy",1000,0,10);
    for(int i=0; i<nbinsy; i++) {
        int bin = energy->FindBin(energy_vector[i]);
        energy->SetBinContent(bin,energy->GetBinContent(bin)+d->GetBinContent(i));
    }   
    canv->cd(2);
    energy->Draw();
    


}
