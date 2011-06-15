#! /usr/bin/env python

import pickle
import sys,math,os,array
from ROOT import gROOT, gStyle, gApplication, gPad, TCanvas, TF1, TH1D, TGraph, TGraphErrors, TFile, TPaveLabel, TLatex, TPaveText, TLegend, TH2D
gROOT.ProcessLine('.x rootlogon2.C')
gROOT.Reset()
gStyle.SetOptFit(1111)
if len(sys.argv) < 2:
    print 'usage',sys.argv[0],'PICKLEFILE'
    quit()

infile = open(sys.argv[1])

data = pickle.load(infile)

aMomRes = {}
aMomResTH = {}
aMomRes2 = {}
aMomRes2TH = {}
aMomResid2 = {}
aMomResid2TH = {}

momVals = []
thVals = []

for key,dat in data.iteritems():
    if momVals.count(dat[0]) == 0: momVals.append(dat[0])
for key,dat in data.iteritems():
    if thVals.count(dat[1]) == 0: thVals.append(dat[1])

momVals.sort()
thVals.sort()
print momVals
print thVals

for mom in momVals:
    aMomRes[mom] = array.array('f')
    aMomResTH[mom] = array.array('f')
for th in thVals:
    aMomRes2[th] = array.array('f')
    aMomRes2TH[th] = array.array('f')
    aMomResid2[th] = array.array('f')
    aMomResid2TH[th] = array.array('f')


for key,dat in data.iteritems():
    mom = dat[0]
    theta = dat[1]
    momRes = math.fabs(dat[-1][1]/dat[-1][0])
    momResid = (dat[-1][0] - mom)/mom
    print mom,theta,momRes
    aMomResTH[mom].append( theta )
    aMomRes[mom].append( momRes )

    aMomRes2TH[theta].append( mom )
    aMomRes2[theta].append( momRes )

    aMomResid2TH[theta].append( mom )
    aMomResid2[theta].append( momResid )

print aMomRes
print aMomResTH
print '############'
print aMomRes2
print aMomRes2TH

momResG = {}
momRes2G = {}
momResid2G = {}

for mom,ARR in aMomRes.iteritems():
    momResG[mom] = TGraph(len(ARR),aMomResTH[mom],ARR)
    momResG[mom].SetName('graph_%f'%mom)
    momResG[mom].SetTitle('graph_%f'%mom)

for th,ARR in aMomRes2.iteritems():
    momRes2G[th] = TGraph(len(ARR),aMomRes2TH[th],ARR)
    momRes2G[th].SetName('graph2_%f'%th)
    momRes2G[th].SetTitle('graph2_%f'%th)
for th,ARR in aMomResid2.iteritems():
    momResid2G[th] = TGraph(len(ARR),aMomResid2TH[th],ARR)
    momResid2G[th].SetName('graphRES2_%f'%th)
    momResid2G[th].SetTitle('graphRES2_%f'%th)

momResC = TCanvas('momResC','')
momResC.SetLogy()
momResDH = TH2D('momResDG','momResDG',2,0,100,2,0.001,0.15)
momResDH.Draw()
momResL = TLegend(0.4,0.6,0.9,0.9)
counter = 0
for mom,GR in momResG.iteritems():
    counter +=1
    GR.SetMarkerStyle(20)
    GR.SetMarkerColor(counter)

    GR.Draw("P")
    momResL.AddEntry(GR.GetName(),'mom=%f'%mom,'p')
momResL.Draw()


momRes2C = TCanvas('momRes2C','')

momRes2C.SetLogx()
momRes2C.SetLogy()
momRes2C.SetGrid(1,1)
momRes2DH = TH2D('momRes2DG','momRes2DG',2,1,120,2,0.00081,0.00399)
momRes2DH.SetStats(0)
momRes2DH.GetYaxis().SetMoreLogLabels()
momRes2DH.GetXaxis().SetTitle('p (GeV/c)')
momRes2DH.GetYaxis().SetTitle('#sigma_{p}/p')
momRes2DH.Draw()
momRes2L = TLegend(0.4,0.6,0.9,0.9)
counter2 = 0
print momRes2G
Gres2Keys = momRes2G.keys()
Gres2Keys.sort()

for th in Gres2Keys:
    counter +=1
    GR = momRes2G[th]
    GR.SetMarkerStyle(20)
    GR.SetMarkerColor(counter)

    GR.Draw("P")
    gPad.Update()
    momRes2L.AddEntry(GR.GetName(),'th=%f'%th,'p')
momRes2L.SetBorderSize(0)
momRes2L.SetFillStyle(0)
momRes2L.SetTextFont(133)
momRes2L.SetTextSize(20)
momRes2L.Draw()


momResid2C = TCanvas('momResid2C','')
momResid2C.SetLogx()
#momResid2C.SetLogy()
momResid2DH = TH2D('momResid2DG','momResid2DG',2,1,55,2,-0.01,0.01)
momResid2DH.SetStats(0)
momResid2DH.Draw()
momResid2L = TLegend(0.4,0.6,0.9,0.9)
counter3 = 0
print momResid2G
Gresid2Keys = momResid2G.keys()
Gresid2Keys.sort()

for th in Gresid2Keys:
    counter3 +=1
    GR = momResid2G[th]
    GR.SetMarkerStyle(20)
    GR.SetMarkerColor(counter3)

    GR.Draw("P")
    gPad.Update()
    momResid2L.AddEntry(GR.GetName(),'th=%f'%th,'p')
momResid2L.SetBorderSize(0)
momResid2L.SetFillStyle(0)
momResid2L.SetTextFont(133)
momResid2L.SetTextSize(20)
momResid2L.Draw()
gPad.Update()


gApplication.Run()
