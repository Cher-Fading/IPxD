#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "../atlasstyle-00-04-02/AtlasUtils.h"
#include "../atlasstyle-00-04-02/AtlasStyle.h"
#include "../atlasstyle-00-04-02/AtlasLabels.h"
#include "../atlasstyle-00-04-02/AtlasStyle.C"

#ifdef __CLING__
// these are not headers - do not treat them as such - needed for ROOT6
#include "../atlasstyle-00-04-02/AtlasLabels.C"
#include "../atlasstyle-00-04-02/AtlasUtils.C"
#endif

#ifdef __CINT__
gROOT->LoadMacro("../atlasstyle-00-04-02/AtlasLabels.C");
gROOT->LoadMacro("../atlasstyle-00-04-02/AtlasUtils.C");
#endif
const Float_t Weight[] = {6.7890E+07 * 6.1692E-06, 6.3996E+05 * 5.8420E-05, 4.7192E+03 * 1.1270E-04, 9.2038E-05 * 2.6602E+01};
Float_t FCal_range[] = {0, 0.063719, 0.14414, 0.289595, 0.525092, 0.87541, 1.36875, 2.04651, 2.98931, 5}; // fcal_cuts options

const int drbin = 20;

const Float_t eta_selection = 2.1;
const bool unique_B = false;
const bool jetTruth = true;

const float min_distx = -0.5;
const float max_distx = -0.25;
const float min_disty = -0.95;
const float max_disty = -0.75;
const float min_distz = -200;
const float max_distz = 200;
const int myColor[] = {kBlue, kViolet, kMagenta, kPink, kOrange, kYellow, kSpring, kGreen, kTeal, kCyan, kAzure, kGray, kGray + 1, kGray + 3};
//const int myColor2[] = {kBlue, kGreen, kRed};

const float min_distxs = -20;
const float max_distxs = 20;
const float min_distys = -20;
const float max_distys = 20;
const float min_distzs = -200;
const float max_distzs = 200;
const int cet[] = {0, 2, 2, 5, 5, 8}; //selected centrality sections
const int cet_N = (sizeof(cet) / sizeof(int)) / 2;

const int dist_bin3 = 100;
const int nQuality = 14;
const int nData = 8;

char Qualitytitle[nQuality][300] = {"No hits in first two layers; expected hit in IBL and b-layer", "No hits in first two layers; expected hit in IBL and no expected hit in b-layer", "No hits in first two layers; no expected hit in IBL and expected hit in b-layer", "No hits in first two layers; no expected hit in IBL and b-layer", "No hit in IBL; expected hit in IBL", "No hit in IBL; no expected hit in IBL", "No hit in b-layer; expected hit in b-layer", "No hit in b-layer; no expected hit in b-layer", "Shared hit in both IBL and b-layer", "At least one shared pixel hits", "Two or more shared SCT hits", "Split hits in both IBL and b-layer", "Split pixel hit", "Good"};
char files[nData][200] = {"WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_50.00_0.5.root", "WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_50.00_1.5.root", "WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_50.00_2.0.root", "WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_50.00_4.0.root", "WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_100.00_0.5.root", "WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_100.00_1.5.root", "WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_100.00_2.0.root", "WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_100.00_4.0.root"};
char data[nData][50] = {"pT_trk > 0.5 GeV; pT_jet > 50 GeV", "pT_trk > 1.5 GeV; pT_jet > 50 GeV", "pT_trk > 2.0 GeV; pT_jet > 50 GeV", "pT_trk > 4.0 GeV; pT_jet > 50 GeV", "pT_trk > 0.5 GeV; pT_jet > 100 GeV", "pT_trk > 1.5 GeV; pT_jet > 100 GeV", "pT_trk > 2.0 GeV; pT_jet > 100 GeV", "pT_trk > 4.0 GeV; pT_jet > 100 GeV"};
//char files[nData][200] = {"BTagCalibRUN2Onl-08-40.root"};
//char data[nData][50] = {"Calibration, pT_jet >  20 GeV"};
const int nFlav = 3;

const char leg[nFlav][10] = {"B", "C", "U"};
const char legends[nFlav][100] = {"b-jet", "c-jet", "Light Jet"};

//const char* dataType = "WorkingDefault";
//const bool PbPb = true;
char Type[2][10] = {"pp", "PbPb"};
//for WorkingDefault it's the reverse PbPb for pp and pp for PbPb
char Quality[nQuality][80] = {"0HitIn0HitNInExp2", "0HitIn0HitNInExpIn", "0HitIn0HitNInExpNIn", "0HitIn0HitNIn", "0HitInExp", "0HitIn", "0HitNInExp", "0HitNIn", "InANDNInShared", "PixShared", "SctShared", "InANDNInSplit", "PixSplit", "Good"};
const float cutMax = 10;
const float cutMin = -10;
const int numCuts = 20;
const float cutInt = (cutMax - cutMin) / numCuts;

const float d0sig = 100;
const float z0sig = 100;

void Draw_templates()
{
    TFile *fi;
    TH1D *graph;
    TH2D *graph2;
    TH1D *graph3;
    TPad *pad1[nFlav + 1];
    TH1F *h1[nFlav + 1];
    TH1D *flav;

    TCanvas *c0 = new TCanvas("c0", "c0");
    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->Divide(2, 2);

    for (int m = 0; m < 2; m++)
    {
        //int m = 0;
        int cent_N = m ? cet_N : 1;
        for (int c = 0; c < cent_N; c++)
        {
            std::string Centrality = m ? Form("Overlay %d %% to %d %%", 10 * cet[2 * c], 10 * cet[2 * c + 1]) : "pp";
            cout << Centrality << endl;
            cout << "m: " << m << "c: " << c << endl;
            //for (int d = 0; d < nData; d++)
            //{
            int d = 0;
            fi = TFile::Open(Form(files[d], Type[m], c), "READ");
            //fi = TFile::Open(files[d], "READ");
            cout << Form(files[d], Type[m], c) << endl;
            //cout << fi << endl;

            //for (int q = 0; q < nQuality; q++)
            //{
            int q = 13;
            TPad *pad0 = (TPad *)c0->cd();
            TH1F *h0 = (TH1F *)pad0->DrawFrame(-20, 1e-5, 40, 1);
            h0->GetXaxis()->SetTitle("d0 significance");
            h0->GetYaxis()->SetTitle("Normalized Fraction");
            h0->SetTitle(Form("D0 Significance Templates %s %s", Type[m], data[d]));
            h0->Draw();
            myText(0.6, 0.8, kBlack, Centrality.c_str(), 0.025);
            myText(0.6, 0.7, kBlack, data[d], 0.025);
            myText(0.6, 0.75, kBlack, Form("Track Grade: %s", Qualitytitle[q]), 0.025);

            pad1[nFlav] = (TPad *)c1->cd(nFlav + 1);
            h1[nFlav] = (TH1F *)pad1[nFlav]->DrawFrame(-40, 1e-6, 60, 10);
            h1[nFlav]->GetXaxis()->SetTitle("z0 significance");
            h1[nFlav]->GetYaxis()->SetTitle("Normalized Fraction");
            h1[nFlav]->GetZaxis()->SetTitle("Normalized Fraction");
            h1[nFlav]->SetTitle(Form("Z0 Significance Templates %s %s", Type[m], data[d]));
            h1[nFlav]->Draw();
            gPad->SetLogy(1);
            myText(0.6, 0.8, kBlack, Centrality.c_str(), 0.025);
            myText(0.6, 0.7, kBlack, data[d], 0.025);
            myText(0.6, 0.75, kBlack, Form("Track Grade: %s", Qualitytitle[q]), 0.025);

            for (int f = 0; f < nFlav; f++)
            {
                c0->cd();
                graph = (TH1D *)(fi->Get(Form("IP2D/AntiKt4HI/%s/%s/SipA0", leg[f], Quality[q]))->Clone());
                //graph = (TH1D *)(fi->Get(Form("IP2D/AntiKt4EMTopo/%s/%s/SipA0", leg[f], Quality[q]))->Clone());

                graph->Scale(1. / graph->GetSumOfWeights());
                graph->SetMarkerColor(myColor[f * 3]);
                graph->SetLineColor(myColor[f * 3]);
                graph->SetMarkerStyle(1);
                graph->Draw("SAME");
                myBoxText(0.7, 0.65 - 0.03 * f, 0.1, myColor[f * 3], 0, legends[f], 0.3, myColor[f * 3], 1, true, 0.03);

                c1->cd(f + 1);
                pad1[f] = (TPad *)c1->cd(f + 1);
                gPad->SetLogy(0);
                gPad->SetLogz(1);
                h1[f] = (TH1F *)pad1[f]->DrawFrame(-40, -40, 60, 60);
                //gStyle->SetTextFont(7);
                h1[f]->GetXaxis()->SetTitle("d0 significance");
                h1[f]->GetYaxis()->SetTitle("z0 significance");
                h1[f]->SetTitle(Form("Z0 and d0 Significance Templates of %s %s", Type[m], data[d]));
                //h1[f]->Draw();

                graph2 = (TH2D *)(fi->Get(Form("IP3D/AntiKt4HI/%s/%s/Sip3D", leg[f], Quality[q]))->Clone());
                //graph2 = (TH2D *)(fi->Get(Form("IP3D/AntiKt4EMTopo/%s/%s/Sip3D", leg[f], Quality[q]))->Clone());

                graph2->Scale(1. / graph2->GetSumOfWeights());
                //graph2->SetMarkerColor(myColor[f]);
                //graph2->SetLineColor(myColor[f]);
                //graph2->SetMarkerStyle(1);
                graph2->Draw("SAME colz");
                //graph2->GetZaxis()->SetTitle("Normalized Fraction");
                graph2->GetZaxis()->SetRangeUser(1e-4, 1);
                //graph2->GetZaxis()->SetTitleSize(.15);
                graph2->GetZaxis()->SetLabelSize(.02);
                //c1->cd(f + 1);
                myText(0.5, 0.65, kBlack, legends[f], 0.03);
                myText(0.5, 0.8, kBlack, Centrality.c_str(), 0.03);
                myText(0.5, 0.75, kBlack, data[d], 0.03);
                myText(0.2, 0.75, kBlack, Form("Track Grade: %s", Qualitytitle[q]), 0.025);

                c1->cd(nFlav + 1);
                gPad->SetLogy(1);
                gPad->SetLogz(0);
                flav = (TH1D *)(graph2->ProjectionY()->Clone());
                flav->Scale(1. / flav->GetSumOfWeights());
                flav->SetMarkerColor(myColor[f * 3]);
                flav->SetLineColor(myColor[f * 3]);
                flav->SetMarkerStyle(1);
                flav->Draw("SAME");
                myBoxText(0.7, 0.65 - 0.03 * f, 0.1, myColor[f * 3], 0, legends[f], 0.3, myColor[f * 3], 1, true, 0.03);
                //template[c][d][q][f]->SetName(Form("%s_%s_%s");
            }
            c0->SetLogy(1);
            c0->SaveAs(Form("new_%s_%s_SipD0_%s.pdf", data[d], Quality[q], Centrality.c_str()));

            c1->SaveAs(Form("new_%s_%s_SipZ0_3D_%s.pdf", data[d], Quality[q], Centrality.c_str()));
            //}

            fi->Close();
            //}
        }
    }
}
