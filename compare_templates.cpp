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
//implemented now is 2 selections combining
const int n1 = 2; //jet pt
const int n2 = 4; //trk pt
//the first 8 are test data and the last is calibration.

char Qualitytitle[nQuality][300] = {"No hits in first two layers; \n expected hit in IBL and b-layer", "No hits in first two layers; \n expected hit in IBL and no expected hit in b-layer", "No hits in first two layers; \n no expected hit in IBL and expected hit in b-layer", "No hits in first two layers; \n no expected hit in IBL and b-layer", "No hit in IBL; \n expected hit in IBL", "No hit in IBL; \n no expected hit in IBL", "No hit in b-layer; \n expected hit in b-layer", "No hit in b-layer; \n no expected hit in b-layer", "Shared hit in both IBL and b-layer", "At least one shared pixel hits", "Two or more shared SCT hits", "Split hits in both IBL and b-layer", "Split pixel hit", "Good"};
//const char *files = "%s%s%srapidity%.1f_tuningIPEval_%.2f_%.1f.root";
const char *files = "%s%s_cent_%d_ip3d_tuning_hi_50k_%.2f_%.1f.root";

//files[n1][n2][200] = {{"WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_50.00_0.5.root", "WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_50.00_1.5.root", "WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_50.00_2.0.root", "WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_50.00_4.0.root"}, {"WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_100.00_0.5.root", "WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_100.00_1.5.root", "WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_100.00_2.0.root", "WorkingDefault%s_cent_%d_ip3d_tuning_hi_50k_100.00_4.0.root"}}

const char *calibration = "BTagCalibRUN2Onl-08-40.root";
const char *datal1 = "pT_jet > %.0f GeV";
const char *datal2 = "pT_trk > %.1f GeV";

const char *dataf1 = "pT_jet_%.0f_GeV";
const char *dataf2 = "pT_jet_%.1f_GeV";

const float ptLim[n1] = {50., 100.};
const float trkLim[n2] = {0.5, 1.5, 2.0, 4.0};

const int nFlav = 3;

const char leg[nFlav][10] = {"C", "U", "B"};
const char legends[nFlav][100] = {"c-jet", "Light Jet", "b-jet"};

const char *dataType = "WorkingDefault";
//const bool PbPb = true;
char Type[2][10] = {"pp", "PbPb"};
char jet[2][20] = {"AntiKt4HI", "AntiKt4EMTopo"};
//for WorkingDefault it's the reverse PbPb for pp and pp for PbPb
char Quality[nQuality][80] = {"0HitIn0HitNInExp2", "0HitIn0HitNInExpIn", "0HitIn0HitNInExpNIn", "0HitIn0HitNIn", "0HitInExp", "0HitIn", "0HitNInExp", "0HitNIn", "InANDNInShared", "PixShared", "SctShared", "InANDNInSplit", "PixSplit", "Good"};
const float cutMax = 10;
const float cutMin = -10;
const int numCuts = 20;
const float cutInt = (cutMax - cutMin) / numCuts;

const float d0sig = 100;
const float z0sig = 100;

void compare_templates(bool comp = true, bool plot = true)
{
    gErrorIgnoreLevel = kWarning;
    TFile *fi;
    TH1D *graph;
    TH1::AddDirectory(kFALSE);
    int msize = 21;

    TCanvas *c0 = new TCanvas("c0", "c0");

    if (plot)
    {
        for (int PbPb = 0; PbPb < 2; PbPb++)
        {
            int cent_N = PbPb ? cet_N : 1;
            for (int c = 0; c < cent_N; c++)
            {
                std::string Centrality;
                std::string titleline;
                std::string textline;
                std::string fileline;
                //cout << Centrality << endl;
                //cout << "m: " << m << "c: " << c << endl;
                int jind;
                for (int d = 0; d < n1 * n2 + 1; d++)
                {
                    int d1 = (int)d / n2;
                    int d2 = (int)d % n2;
                    if (d < n1 * n2)
                    {
                        fi = TFile::Open(Form(files, dataType, Type[PbPb], c, ptLim[d1], trkLim[d2]), "READ");
                        jind = 0;
                        titleline = Form("D0 Significance Templates %s %s %s", Type[PbPb], Form(datal1, ptLim[d1]), Form(datal2, trkLim[d2]));
                        Centrality = PbPb ? Form("Overlay %d %% to %d %%", 10 * cet[2 * c], 10 * cet[2 * c + 1]) : "pp";
                        textline = Form("%s; %s", Form(datal1, ptLim[d1]), Form(datal2, trkLim[d2]));
                        fileline = Form("new_%s_%s_%s_SipA0", Form(dataf1, ptLim[d1]), Form(dataf1, ptLim[d1]), Centrality.c_str());
                    }
                    else
                    {
                        fi = TFile::Open(calibration, "READ");
                        jind = 1;
                        titleline = "D0 Significance Templates for Calibration";
                        Centrality = "Calibration";
                        textline = calibration;
                        fileline = Form("new_calibration_SipA0");
                    }
                    cout << "reading file: " << fi->GetName() << endl;

                    for (int q = 0; q < nQuality; q++)
                    {
                        TPad *pad0 = (TPad *)c0->cd();
                        TH1F *h0 = (TH1F *)pad0->DrawFrame(-40, 1e-6, 60, 0.2);
                        h0->GetXaxis()->SetTitle("d0 significance");
                        h0->GetYaxis()->SetTitle("Normalized Fraction");
                        h0->SetTitle(Form("%s", titleline.c_str()));
                        h0->Draw();
                        myText(0.6, 0.8, kBlack, Centrality.c_str(), 0.025);
                        myText(0.6, 0.7, kBlack, textline.c_str(), 0.025);
                        myText(0.6, 0.75, kBlack, Form("Track Grade: %s", Qualitytitle[q]), 0.025);

                        for (int f = 0; f < nFlav; f++)
                        {
                            c0->cd();
                            //cout << fi->GetName() << endl;
                            //cout << Form("IP2D/%s/%s/%s/SipA0", jet[jind], leg[f], Quality[q]) << endl;
                            graph = (TH1D *)(fi->Get(Form("IP2D/%s/%s/%s/SipA0", jet[jind], leg[f], Quality[q]))->Clone());
                            graph->Scale(1. / graph->GetSumOfWeights());
                            graph->SetMarkerColor(myColor[f * 3]);
                            graph->SetLineColor(myColor[f * 3]);
                            graph->SetMarkerStyle(1);
                            graph->Draw("SAME");
                            myBoxText(0.7, 0.65 - 0.03 * f, 0.1, myColor[f * 3], 0, leg[f], 0.3, myColor[f * 3], 1, true, 0.03);
                        }
                        c0->SetLogy();
                        c0->SaveAs(Form("%s_%s.pdf", fileline.c_str(), Quality[q]));
                    }
                    fi->Close();
                }
            }
        }
        //TFile *cal = TFile::Open(cal, "READ");

        cout << "\n"
             << endl;
        cout << "Now compare" << endl;
    }
    if (comp)
    {
        for (int q = 0; q < nQuality; q++)
        {
            //int q = 13;
            for (int f = 0; f < nFlav; f++)
            {
                //int f = 0;
                for (int d2 = 0; d2 < n2; d2++)
                {
                    for (int PbPb = 0; PbPb < 2; PbPb++)
                    {
                        int cent_N = PbPb ? cet_N : 1;
                        for (int c = 0; c < cent_N; c++)
                        {
                            TPad *pad0 = (TPad *)c0->cd();
                            TH1F *h0 = (TH1F *)pad0->DrawFrame(-40, 1e-6, 60, 0.2);
                            std::string Centrality = PbPb ? Form("Overlay %d %% to %d %%", 10 * cet[2 * c], 10 * cet[2 * c + 1]) : "pp";
                            h0->GetXaxis()->SetTitle("d0 significance");
                            h0->GetYaxis()->SetTitle("Normalized Fraction");
                            h0->SetTitle(Form("Comparing With Calibration for Different jet pT given Trk pT > %.1f GeV", trkLim[d2]));
                            h0->Draw();
                            myText(0.6, 0.75, kBlack, Form("%s %s", Centrality.c_str(),legends[f]), 0.025);
                            myText(0.6, 0.7, kBlack, Form(datal2, trkLim[d2], 0.025));
                            myText(0.6, 0.65, kBlack, Form("Track Grade: %s", Qualitytitle[q]), 0.025);


                            for (int d1 = 0; d1 < n1; d1++)
                            {
                                fi = TFile::Open(Form(files, dataType, Type[PbPb], c, ptLim[d1], trkLim[d2]), "READ");
                                //cout << fi->GetName() << endl;
                                graph = (TH1D *)(fi->Get(Form("IP2D/%s/%s/%s/SipA0", jet[0], leg[f], Quality[q]))->Clone());
                                //cout << Form("IP2D/%s/%s/%s/SipA0", jet[0], leg[f], Quality[q]) << endl;
                                graph->Scale(1. / graph->GetSumOfWeights());
                                //cout << graph->GetEntries() << endl;
                                graph->SetMarkerColor(myColor[d1 * 3]);
                                graph->SetLineColor(myColor[d1 * 3]);
                                graph->SetMarkerStyle(msize);
                                graph->SetMarkerSize(1);
                                c0->cd();
                                graph->Draw("SAME");
                                myBoxText(0.6, 0.8 + 0.03 * d1, 0.1, myColor[d1 * 3], 0, Form(dataf1, ptLim[d1]), 0.3, myColor[d1 * 3], msize, true, 0.03);
                                //c0->SetLogy();
                                fi->Close();
                                //return;
                            }
                            fi = TFile::Open(calibration, "READ");
                            graph = (TH1D *)(fi->Get(Form("IP2D/%s/%s/%s/SipA0", jet[1], leg[f], Quality[q]))->Clone());
                            graph->Scale(1. / graph->GetSumOfWeights());
                            graph->SetMarkerColor(myColor[n1 * 3]);
                            graph->SetLineColor(myColor[n1 * 3]);
                            graph->SetMarkerStyle(msize);
                            graph->SetMarkerSize(1);
                            graph->Draw("SAME");
                            myBoxText(0.6, 0.8 + 0.03 * n1, 0.1, myColor[n1 * 3], 0, "Calibration", 0.3, myColor[n1 * 3], msize, true, 0.03);
                            c0->SetLogy();
                            c0->SaveAs(Form("CalComp_%s_%s_%s_%s_%d.pdf", Centrality.c_str(), Form(dataf2, trkLim[d2]), Quality[q], leg[f],msize));
                            fi->Close();
                        }
                    }
                }
            }
        }
    }
}
