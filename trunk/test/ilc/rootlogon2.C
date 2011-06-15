{
//  gROOT->ProcessLine(".L ../../phast.7.102/myLib/libMY.so");

//std::cout << endl << "Executing rootlogon.C" << endl << endl;
//std::cout << "For B/W plots use: gROOT->SetStyle(\"bw\");" << endl;
//std::cout << "For col plots use: gROOT->SetStyle(\"col\");" << endl << endl;
//std::cout << "To force objects read from file to update their attributes" << endl;
//std::cout << "add:               gROOT->ForceStyle();" << endl << endl;

// Black/White Style
TStyle *bwStyle= new TStyle("bw","B/W publication plot style");

// use plain black on white colors
bwStyle->SetFrameBorderMode(0);
bwStyle->SetCanvasBorderMode(0);
bwStyle->SetCanvasDefH(600);
bwStyle->SetCanvasDefW(600);
bwStyle->SetPadBorderMode(0);
bwStyle->SetPadColor(0);
bwStyle->SetCanvasColor(0);
bwStyle->SetStatColor(0);

// create b/w palette with 50 gray shades
const int ncol=50;
int palette[ncol];
double r,g,b;
for (int i=0;i<ncol;i++) {
 r = 1-(i/(ncol*1.0));
 g = r;
 b = r;
 if (!gROOT->GetColor(300+i)){
   TColor *color = new TColor(300+i,r,g,b,"");
 } else {
   TColor *color = gROOT->GetColor(300+i);
   color->SetRGB(r,g,b);
 }
 palette[i] = 300+i;
}
bwStyle->SetPalette(ncol,palette);

// set the paper & margin sizes
bwStyle->SetPaperSize(20,26);
bwStyle->SetPadTopMargin(0.07);
bwStyle->SetPadRightMargin(0.08);
bwStyle->SetPadBottomMargin(0.12);
bwStyle->SetPadLeftMargin(0.12);

// set fonts
//Int_t bwFont = 132;  // use variable size times
//Float_t bwSize = 0.07; // size in percent of pad size
Int_t bwFont = 133; // use fixed size times
Float_t bwSize = 28; // size in pixels
bwStyle->SetTextFont(bwFont);
bwStyle->SetTextSize(bwSize);
bwStyle->SetLabelFont(bwFont,"xyz");
bwStyle->SetLabelSize(bwSize-2,"xyz");
bwStyle->SetTitleFont(bwFont,"xyz");
bwStyle->SetTitleSize(bwSize,"xyz");
bwStyle->SetTitleOffset(1.14,"x");
bwStyle->SetTitleOffset(1.23,"y");

// restrict the number of digits in labels
 TGaxis::SetMaxDigits(3); 

// use bold lines and markers
bwStyle->SetMarkerStyle(20);
bwStyle->SetHistLineWidth(1.85);
bwStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
bwStyle->SetFuncWidth(1.85);

// get rid of X error bars and y error bar caps
bwStyle->SetErrorX(0.001);

// do not display any of the standard histogram decorations
bwStyle->SetOptTitle(0);
bwStyle->SetOptStat(0);
bwStyle->SetOptFit(0);

// put tick marks on top and RHS of plots
//bwStyle->SetPadTickX(1);
//bwStyle->SetPadTickY(1);

// Color Style
TStyle *colStyle= new TStyle("col","Color publication plot style");
bwStyle->Copy(*colStyle);
colStyle->SetPalette(0);

// Use plain style as default
gROOT->SetStyle("bw");
gStyle->SetPalette(1);
gStyle->SetOptStat(1111111);
//gStyle->SetPadTickX(1);
//gStyle->SetPadTickY(1);
}
