//! ////////////////////////////////////////////////////////////////////////////
//! @file: auxTools.h                                                            //
//! @author: Jacob Smith (smithja)                                            //
//! Last edited: October 15th, 2025                                           //   
//!                                                                           // 
//! @details: As I've done more HEP research, I've compilied these functions  //
//! that I've found myself routinely writing. If you're reading this and have //
//! found the functions in this file useful, please share them with           //
//! colleagues and feel free to add to this file with your own contributions. //
//!                                                                           //
//!               `````             `               ```                 ````` //
//!          `````````````````````````````      ```````````         ````````` //
//! ``````` ```````````````  `````````` ````````````````       ``` ```` ````` //
//! ```````` ```       ````` ++.   .#@@@@  ````````   ````` ``  ```   ``````` //
//! ``` `````   ```````   :@@@@@@@@@@@@@@@@@@@@@@::::::,:,,::::   `````` .@@: //  
//! ` ;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##@@@@@@@ //
//! @@@@@@@@@@@@@@@@@@@@@@@@@#';,,@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ //
//! @@@@@@@@                      @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ //
//! @@@@@@@@ @@@@@@@@@@@@@@@@@@   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ //
//! @@ #@@@@ @@@@@@@@#@@@@@@@@@   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ //
//! @   @@@@ @@@@@`  @@@@@@@ :;;  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ @@ //
//!      @@@ @@@`   +@@@ @@ ;;;:.  ;,#@@@@@@@@@@@@@@@@ _ _ _ _ _ _ _ _ _ _ _  //
//!      @@@ @@`    @@@@  @@:    :;;:#@@@@@@@@@@@@@@@ / Wait, there's a     | //
//!     `@@@ @      @@@    ;;;;;;;;`.::: @@@@@@@@@@@ /    function for that?| //
//!      @@@    '   @@:     ;;;:. ::;;;; ,`+@@@@@@@ /  _ _ _ _ _ _ _ _ _ _ _/ //
//!  '   @,@   '''      '''   :.@++@ : @    .@@@@@@/_/@@@@@@@@  @@@@@@@@@     //
//!  '  @@ ., '''''    ''     ;@,   @ .#@    `      .@@@@@@@:  `@@@@@@@@   ;  //
//! `'  @@  `'''''''; ''      ;@    @`;:.;;;;:       @@@@@@  `  @@@@@@@   ''  //
//!  '  @@   ;'''''''''       :.@@@# ::;:;;;;::     @@@@@#  ',  .@@@@@`  `''. //
//! ''  @'    '''''''';      ;;;;;;;;;:.``;::::;;..@'@@@   ' ......` @   ;''' //
//! '  .@     '''''''.,,   `::;;;;;:;:;;;: @@++#@@@  @@   '' ':   .' ++  '''' //
//! '   @     '''''';.,, , :;;::;::;;;;;;; @@@@@@@   @ ` '''`'++++++ ''@@@@@@ //
//! '   `     ''''''',,, , ::: : :;;;;;;::, @@@@@+      ''@@ '++++++ ,@@@@@@@ //
//! ',     '  ''''';;.,, , :::.: :;;;;;;;::,@@@@@,  ;  :@@@@ ;++++'+ @@@@@@@@ //
//! ''     '  ;;;;;;' ,, ,,,::;.:`;;;:;;;:.;'@@@@;  '' ''@@@@@@@@@@@@@@@@@@@@ //
//! '';   ''  ;;;;;;; ,,,,. ,:::`:,     ,`;` `@@@#  ;''''''''''+''''''''''''' //
//! '''. '''  ::::,,'; ,,,,,, `:::::` :: ,,`;:@'@#  .'''''''''''''''''''''''' //
//!  ''';'';  ;;;;''';;; ,,,,,,,,,,,, :;;`  :;:: +  `''''''''+''+''  '''''''' //
//! :'''''';  ;;;;  '';;'   ````      :;::` `  @@   ''''''''    ``   ```````` //
//! ;'''''''  ,;;   ;;;;;;;`,,.....   ``       `    '''  '''         ```````  //
//! ''''''''   ;     ;';;;;.,,.;;''   ,,;.   ;;;  '': ; ''''         ``````   //
//! '.'''''''     ,  ';' '';,,.;;;;   ,,;,   ;;';;;;;' `.;'`        `````     //
//! ' ''`' '''    '   '  ;;;,,.;;;'   ,,;;   ;;;;;;;;;;;;   ` .'    ````      //
//! ':'';''''''   ''     ;;;,,.;;;;   ,,;;   ;;;;;;;';;`     '''    ``        //
//! '',''`':''''; '''`   ;;;` ;;;;;'; ,,;;  ';;;;;;;;:     :''''           ;' //
//! ''''''''''''''''''' ``;'';;;;;;;;';;;;;;;;;;;;;;'   ` ;''''';        '''' //
//! ''''''''''''''''''', ;;;;;;;;;;;';;;;;;;;;;;;;;'      '''''''',    '''''' //
//! '''''''''''''''''; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;'     ''''''''''':`''''''' //
//!                                                                           //
//!                                                                           //
//! @par Disclaimer                                                           //
//! This code is available for non-commercial use under the Creative Commons  //
//! Attribution-NonCommercial (CC BY-NC) 4.0. You must continue to use this   //
//! Creative Commons license if you copy, repurpose, make additions to, or    //
//! otherwise use any of the contents of this file.                           //             
//! ////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>     
#include <vector>
#include <sstream>    //! to treat strings as streams
#include <fstream>    //! to work with files
#include <iostream>   //! to make use of std::cin and std::cout
#include <iomanip>    //! to format how data/output is displayed
#include <cmath>      //! to do mathematical things, e.g. log10()
#include <sys/stat.h> //! to create directories
#include <errno.h>    //! to diagnose errors

#include "TROOT.h"

//! ROOT objects used in this file:
#include "TCanvas.h"
#include "TColor.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"


/** @brief Create a directory if it doesn't already exist.
 *  @param dirname Pointer to the directory name.
 */
void createDir(const std::string* dirname) {
    const char* char_dirname = dirname->c_str();
    int status = mkdir(char_dirname, 0777);
    
    if (status == 0) {
        std::cout << "Directory " << *dirname << " created using mkdir." << std::endl;
    } else if (errno == EEXIST) {
        std::cerr << "Directory " << *dirname << " already exists." << std::endl;
    } else {
        std::cerr << "Error creating directory " << *dirname << " using mkdir.";
    }
}

/** @brief Draw a legend with efficiency and purity histograms.
 *  @param h1 Pointer to first histogram.
 *  @param name1 Name for first histogram.
 *  @param h2 Pointer to second histogram.
 *  @param name2 Name for second histogram.
 */
void DrawEffPurLegend(TH1* h1, const std::string& name1, TH1* h2, const std::string& name2) {
    TLegend* leg = new TLegend(0.77, 0.15, 0.92, 0.20);
    leg->AddEntry(h1, name1.c_str(), "l");
    leg->AddEntry(h2, name2.c_str(), "l");
    leg->SetTextSize(0.03);
    leg->Draw();
}

/** @brief Draw a legend for integrated efficiency/purity plots.
 *  @param h1 First histogram.
 *  @param name1 Name for h1.
 *  @param h2 Second histogram.
 *  @param name2 Name for h2.
 */
void DrawIntEffPurLegend(TH1* h1, const std::string& name1, TH1* h2, const std::string& name2) {
    TLegend* leg = new TLegend(0.77, 0.15, 0.92, 0.20);
    leg->AddEntry(h1, name1.c_str(), "l");
    leg->AddEntry(h2, name2.c_str(), "l");
    leg->SetTextSize(0.03);
    leg->Draw();
}

/** @brief Draw a legend for a single purity histogram.
 *  @param hist Histogram pointer.
 *  @param name Name for the legend entry.
 */
void DrawPurLegend(TH1* hist, const std::string& name) {
    TLegend* leg = new TLegend(0.77, 0.16, 0.92, 0.20);
    leg->AddEntry(hist, name.c_str(), "l");
    leg->SetTextSize(0.03);
    leg->Draw();
}

/** @brief Show signal/background breakdown and POT tag.
 *  @param hSignal Signal histogram.
 *  @param hBackground Background histogram.
 *  @param POTTag String label for POT.
 *  @param textSize Size of text drawn on canvas.
 */
void DrawSigBkgIntText(TH1* hSignal, TH1* hBackground, std::string POTTag, float textSize) {
    float iSignal = hSignal->Integral();
    float iBackground = hBackground->Integral();
    float pSignal = 100. * iSignal / (iSignal + iBackground);
    float pBackground = 100. * iBackground / (iSignal + iBackground);

    TPaveText* pText1 = new TPaveText(0.45, 0.825, 0.55, 0.875, "brNDC");
    TText* text1 = pText1->AddText(POTTag.c_str());
    text1->SetTextSize(textSize);
    pText1->SetBorderSize(0);
    pText1->SetFillStyle(0);
    pText1->Draw();

    TPaveText* pText2 = new TPaveText(0.15, 0.725, 0.25, 0.80, "brNDC");
    TText* text2 = pText2->AddText(Form("Signal: %2.f = %2.f %%", iSignal, pSignal));
    text2->SetTextAlign(11);
    text2->SetTextSize(textSize);
    TText* text3 = pText2->AddText(Form("Background: %2.f = %2.f %%", iBackground, pBackground));
    text3->SetTextAlign(11);
    text3->SetTextSize(textSize);
    pText2->SetBorderSize(0);
    pText2->SetFillStyle(0);
    pText2->Draw();
}

/** @brief Dim the fill color of a histogram (opacity or saturation).
 *  @param hist Histogram to modify.
 *  @param useAlpha If true, use alpha blending (opacity). Else, reduce saturation.
 *  @param dim Fraction of original color to retain.
 */
void fillWithDimmedColor(TH1* hist, bool useAlpha=false, float dim=0.8) {
    if (useAlpha) {
        hist->SetFillColorAlpha(hist->GetLineColor(), dim);
        return;
    }

    TColor* color = gROOT->GetColor(hist->GetLineColor());
    float R, G, B, hR, hG, hB, hHue, hSat, hVal;
    color->GetRGB(hR, hG, hB);
    color->RGB2HSV(hR, hG, hB, hHue, hSat, hVal);
    color->HSV2RGB(hHue, dim * hSat, hVal, R, G, B);
    hist->SetFillColor(color->GetColor(R, G, B));
}

/** @brief Extracts file name (without extension) from a full path.
 *  @param fileName Input path string.
 *  @return Filename without directory or extension.
 */
std::string extractFileName(const std::string& fileName) {
    std::string res = fileName;
    size_t pathPos = res.rfind('/');
    if (pathPos != std::string::npos) res = fileName.substr(pathPos + 1);
    size_t lastDotPos = res.rfind('.');
    if (lastDotPos != std::string::npos) return res.substr(0, lastDotPos);
    return res;
}

/** @brief Returns histogram showing selection efficiency.
 *  @param hSelSignal Selected signal.
 *  @param hSelBackground Selected background.
 *  @param hTrueSignal Ground truth signal.
 */
TH1D* getEffHistogram(TH1* hSelSignal, TH1* hSelBackground, TH1* hTrueSignal) {
    TH1D* hTotal = (TH1D*)hSelSignal->Clone();
    hTotal->Add(hSelBackground);
    TH1D* hEfficiency = (TH1D*)hSelSignal->Clone();
    hEfficiency->Divide(hTrueSignal);
    return hEfficiency;
}

/** @brief Computes maximum Y value across a vector of histograms.
 *  @param histos Vector of TH1 pointers.
 *  @return Maximum Y value.
 */
float getHistMaxY(std::vector<TH1*> histos) {
    float max = 0.;
    for (unsigned int i = 0; i < histos.size(); i++) {
        float thisMax = histos[i]->GetMaximum();
        if (thisMax > max) max = thisMax;
    }
    return max;
}

/** @brief Get the visually "opposite" color to the input color.
 *  @param colorIndex Input color index (ROOT).
 *  @return New color index with opposite hue.
 */
Int_t getOppositeColor(Int_t colorIndex) {
    TColor* inputColor = gROOT->GetColor(colorIndex);
    if (!inputColor) {
        Warning("getOppositeColor", "Invalid color index %d; returning kBlack", colorIndex);
        return kBlack;
    }

    Float_t r, g, b;
    inputColor->GetRGB(r, g, b);

    const Double_t epsilon = 0.05;
    if ((r < epsilon && g < epsilon && b < epsilon) ||
        (r > 1.0 - epsilon && g > 1.0 - epsilon && b > 1.0 - epsilon)) {
        return kBlack;
    }

    Float_t h, s, v;
    TColor::RGB2HSV(r, g, b, h, s, v);
    h = fmod(h + 0.5, 1.0);
    TColor::HSV2RGB(h, s, v, r, g, b);

    Int_t newIndex = TColor::GetFreeColorIndex();
    new TColor(newIndex, r, g, b);
    return newIndex;
}

/** @brief Return histogram of selection purity.
 *  @param hSelSignal Selected signal.
 *  @param hSelBackground Selected background.
 */
TH1D* getPurHistogram(TH1* hSelSignal, TH1* hSelBackground) {
    TH1D* hTotal = (TH1D*)hSelSignal->Clone();
    hTotal->Add(hSelBackground);
    TH1D* hPurity = (TH1D*)hSelSignal->Clone();
    hPurity->Divide(hTotal);
    return hPurity;
}

/** @brief Fix string to be LaTeX-safe.
 *  @param str String to fix.
 *  @return Fixed string.
 */
TString fixLatexName(TString str) {
    std::vector<TString> in  = {"#", " ", ".", "_"};
    std::vector<TString> out = {"",  "",  "",  "\\_"};

    for (size_t i = 0; i < in.size(); ++i) str.ReplaceAll(in[i], out[i]);
    return str;
}

/** @brief Format 1D histogram appearance.
 *  @param hist Histogram pointer.
 *  @param color Line/marker color.
 *  @param lineStyle Line style.
 *  @param lineWidth Line width.
 *  @param markerStyle Marker style (default 8).
 *  @param markerSize Marker size (default 8).
 */
void formatHist(TH1* hist, Color_t color, Style_t lineStyle, int lineWidth,
                Style_t markerStyle = 8, double markerSize = 8) {
    hist->SetLineColor(color);
    hist->SetLineStyle(lineStyle);
    hist->SetLineWidth(lineWidth);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerSize(markerSize);
}

/** @brief Display "ICARUS Sim." or "ICARUS Data" tag on canvas.
 *  @param isSim If true, label as simulation.
 */
void isSimulation(bool isSim = true) {
    TLatex* prelim = new TLatex(0.95, 0.975, isSim ? "ICARUS Sim." : "ICARUS Data");
    prelim->SetTextColor(kGray + 1);
    prelim->SetNDC();
    prelim->SetTextSize(1 / 20.);
    prelim->SetTextAlign(32);
    prelim->Draw();
}

/** @brief Label top-left of canvas with a tag string.
 *  @param str Label string.
 */
void myCornerLabel(const std::string& str) {
    TLatex* CornLab = new TLatex(0.10, 0.93, str.c_str());
    CornLab->SetTextColor(kGray + 1);
    CornLab->SetNDC();
    CornLab->SetTextSize(1 / 20.);
    CornLab->SetTextAlign(11);
    CornLab->Draw();
}

/** @brief Create a POT string tag with scientific notation.
 *  @param value POT value.
 *  @return Formatted string tag.
 */
std::string makePOTTag(double value) {
    if (value == 0.0) return "0 POT";

    int exp = static_cast<int>(std::floor(std::log10(std::fabs(value))));
    int expGroup = (exp / 1) * 1;
    double mantissa = value / std::pow(10.0, expGroup);

    std::ostringstream ss;
    ss << std::fixed << std::setprecision(2) << mantissa
       << " #times 10^{" << expGroup << "} POT";
    return ss.str();
}

/** @brief Parse a structured histogram title string into components.
 *  @param title Histogram title string.
 *  @return TitleInfo struct.
 */
struct TitleInfo {
    std::string variable;
    std::string unit;
    std::string selection;
    std::string type;
};

TitleInfo parseTitle(const std::string& title) {
    TitleInfo info;
    size_t vpos = title.find("var=");
    size_t upos = title.find("_unit=");
    size_t spos = title.find("_sel=");
    size_t tpos = title.find("_type=");
    size_t epos = title.find_last_of("__");

    if (vpos == std::string::npos) return info;

    if (upos != std::string::npos) //! if title includes "_unit=" tag
        info.variable = title.substr(vpos + 4, upos - (vpos + 4));
    else //! if title does not have a "_unit=" tag; @note this is bad practice!!
        info.variable = title.substr(vpos + 4, spos - (vpos + 4));

    if (upos != std::string::npos && spos != std::string::npos)
        info.unit = title.substr(upos + 6, spos - (upos + 6));

    if (spos != std::string::npos && tpos != std::string::npos)
        info.selection = title.substr(spos + 5, tpos - (spos + 5));

    if (tpos != std::string::npos)
        info.type = title.substr(tpos + 6, epos - (tpos + 6));

        
    return info;
}

/** @brief Calculate slice start/end indices for splitting tasks across ranks.
 *  @param total Total number of items.
 *  @param nRanks Number of workers.
 *  @param rank This worker's rank.
 */
std::pair<int, int> sliceIndices(int total, int nRanks, int rank) {
    int base = total / nRanks;
    int rem  = total % nRanks;

    int start = rank * base + std::min(rank, rem);
    int end   = start + base + (rank < rem ? 1 : 0);
    end = std::min(end, total);
    return {start, end};
}

/** @brief Split string by a given substring delimiter.
 *  @param s String to split.
 *  @param delimiter Delimiter substring.
 *  @return Vector of string tokens.
 */
std::vector<std::string> splitBySubstring(const std::string& s, const std::string& delimiter) {
    std::vector<std::string> tokens;
    size_t prev_pos = 0;
    size_t current_pos = s.find(delimiter, prev_pos);

    while (current_pos != std::string::npos) {
        tokens.push_back(s.substr(prev_pos, current_pos - prev_pos));
        prev_pos = current_pos + delimiter.length();
        current_pos = s.find(delimiter, prev_pos);
    }
    tokens.push_back(s.substr(prev_pos));
    return tokens;
}

/** @brief Create a canvas with two overlaid pads.
 *  @param c1 Pointer to canvas pointer.
 *  @param pad1 Pointer to first pad.
 *  @param pad2 Pointer to second pad.
 */
void splitCanvasInTwo(TCanvas*& c1, TPad*& pad1, TPad*& pad2) {
    c1 = new TCanvas("c1", "c1", 1000, 1400);
    c1->cd();

    pad1 = new TPad("pad1", "pad1", 0, 0, 1, 1);
    pad1->SetTopMargin(0.1);
    pad1->SetBottomMargin(0.4);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.03);
    pad1->SetFillStyle(0);
    pad1->Draw();

    pad2 = new TPad("pad2", "pad2", 0, 0, 1, 1);
    pad2->SetTopMargin(0.6);
    pad2->SetBottomMargin(0.1);
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.03);
    pad2->SetFillStyle(0);
    pad2->Draw();
}
