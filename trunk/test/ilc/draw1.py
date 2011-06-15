#! /usr/bin/env python

import pickle
import sys,math,os,array
import sys,math,os,array
from ROOT import gROOT, gStyle, gApplication, gPad, TCanvas, TF1, TH1D, TGraph, TGraphErrors, TFile, TPaveLabel, TLatex, TPaveText, TLegend
#gROOT.ProcessLine('.x rootlogon2.C')
gROOT.Reset()
gROOT.SetStyle('Plain')
gStyle.SetPadLeftMargin(0.16);
gStyle.SetTitleOffset(1.6,"y");
gStyle.SetOptFit(1111)
if len(sys.argv) < 2:
    print 'usage',sys.argv[0],'INFILES'
    quit()

infiles = {}
canv = {}
canv2 = {}
tmpHists = {}
hists = {}
funcs = {}

data = {}

labels = {}
labels[0] = 'q/p'
labels[1] = 'du/dw'
labels[2] = 'dv/dw'
labels[3] = 'u'
labels[4] = 'v'

for file in sys.argv[1:]:
    mom = float(file.split('/')[-1].split('_')[1])
    theta = float(file.split('/')[-1].split('_')[2].split('.')[0]+'.'+file.split('/')[-1].split('_')[2].split('.')[1])

    id = '%.1f_%.1f'%(mom,theta)
    
    data[id] = []
    data[id].append(mom)
    data[id].append(theta)

    print mom,theta
    print file
    infiles[file] = TFile(file)
    canv[file] = TCanvas('canv_%s'%file,'canv_%s'%file)
    canv[file].Divide(3,2)
    t = gROOT.FindObject('t')
    data[id].append( [] )
    for par in range(0,5):
        canv[file].cd(par+1)
        funcs['f_%.1f_%.1f_pull%d'%(mom,theta,par)] = TF1('f_%.1f_%.1f_pull%d'%(mom,theta,par),'[0]*exp(-0.5*((x-[1])/[2])^2)',-6,6.)
        funcs['f_%.1f_%.1f_pull%d'%(mom,theta,par)].SetParameters(20.,0,1)
        hists['h_%.1f_%.1f_pull%d'%(mom,theta,par)] = TH1D('h_%.1f_%.1f_pull%d'%(mom,theta,par),'',200,-10,10.)    
        t.Draw("(stREC[%d]-stMCT[%d])/covREC[%d]>>h_%.1f_%.1f_pull%d"%(par,par,par,mom,theta,par))
        hists['h_%.1f_%.1f_pull%d'%(mom,theta,par)].Fit('f_%.1f_%.1f_pull%d'%(mom,theta,par),'rl')
        hists['h_%.1f_%.1f_pull%d'%(mom,theta,par)].GetXaxis().SetTitle('pull for %s'%labels[par])
        hists['h_%.1f_%.1f_pull%d'%(mom,theta,par)].GetYaxis().SetTitle('number of tracks')
        data[id][-1].append( [funcs['f_%.1f_%.1f_pull%d'%(mom,theta,par)].GetParameter(1),
                              funcs['f_%.1f_%.1f_pull%d'%(mom,theta,par)].GetParameter(2)] ) 

    gPad.Update()
    canv[file].SaveAs('pulls_%.1f_%.1f.pdf'%(mom,theta))

    canv2[file] = TCanvas('canv2_%s'%file,'canv2_%s'%file)
    canv2[file].Divide(3,2)

    data[id].append( [] )
    for par in range(0,5):
        canv2[file].cd(par+1)
        tmpName = 'temp_%.1f_%.1f_%d'%(mom,theta,par)
        t.Draw("(stREC[%d]-stMCT[%d])>>%s"%(par,par,tmpName))
        tmpHists[tmpName]= gROOT.FindObject(tmpName)
        min = tmpHists[tmpName].GetXaxis().GetXmin()
        max = tmpHists[tmpName].GetXaxis().GetXmax()
        if math.fabs(min)>max: min=-1.0*max
        else: max=-1.0*min
        funcs['f_%.1f_%.1f_res%d'%(mom,theta,par)] = TF1('f_%.1f_%.1f_res%d'%(mom,theta,par),'[0]*exp(-0.5*((x-[1])/[2])^2)',
                                                         min,max)
        funcs['f_%.1f_%.1f_res%d'%(mom,theta,par)].SetParameters(20.,0,(tmpHists[tmpName].GetXaxis().GetXmax()-tmpHists[tmpName].GetXaxis().GetXmin())/10.)
        hists['h_%.1f_%.1f_res%d'%(mom,theta,par)] = TH1D('h_%.1f_%.1f_res%d'%(mom,theta,par),'',100,
                                                          min,max)
        hists['h_%.1f_%.1f_res%d'%(mom,theta,par)].GetXaxis().SetTitle('residual for %s'%labels[par])
        hists['h_%.1f_%.1f_res%d'%(mom,theta,par)].GetYaxis().SetTitle('number of tracks')

        t.Draw("(stREC[%d]-stMCT[%d])>>h_%.1f_%.1f_res%d"%(par,par,mom,theta,par))
        hists['h_%.1f_%.1f_res%d'%(mom,theta,par)].Fit('f_%.1f_%.1f_res%d'%(mom,theta,par),'rl')
        data[id][-1].append( [funcs['f_%.1f_%.1f_res%d'%(mom,theta,par)].GetParameter(1),
                              funcs['f_%.1f_%.1f_res%d'%(mom,theta,par)].GetParameter(2)] ) 


    canv2[file].cd(6)
    tmpName = 'temp_%.1f_%.1f_mom'%(mom,theta)
    t.Draw("1./abs(stREC[0])>>%s"%(tmpName))
    tmpHists[tmpName]= gROOT.FindObject(tmpName)
    funcs['f_%.1f_%.1f_mom'%(mom,theta)] = TF1('f_%.1f_%.1f_mom'%(mom,theta),'[0]*exp(-0.5*((x-[1])/[2])^2)',
                                               tmpHists[tmpName].GetMean()-3.*tmpHists[tmpName].GetRMS(),
                                               tmpHists[tmpName].GetMean()+3.*tmpHists[tmpName].GetRMS())
    funcs['f_%.1f_%.1f_mom'%(mom,theta)].SetParameters(20.,tmpHists[tmpName].GetMean(),tmpHists[tmpName].GetRMS())
    hists['h_%.1f_%.1f_mom'%(mom,theta)] = TH1D('h_%.1f_%.1f_mom'%(mom,theta),'',100,
                                                tmpHists[tmpName].GetMean()-3.*tmpHists[tmpName].GetRMS(),
                                                tmpHists[tmpName].GetMean()+3.*tmpHists[tmpName].GetRMS())
    hists['h_%.1f_%.1f_mom'%(mom,theta)].GetXaxis().SetTitle('reconstructed momentum (GeV/c)')
    hists['h_%.1f_%.1f_mom'%(mom,theta)].GetYaxis().SetTitle('number of tracks')
    t.Draw("1./abs(stREC[0])>>h_%.1f_%.1f_mom"%(mom,theta))
    hists['h_%.1f_%.1f_mom'%(mom,theta)].Fit('f_%.1f_%.1f_mom'%(mom,theta),'rl')
    data[id].append( [funcs['f_%.1f_%.1f_mom'%(mom,theta)].GetParameter(1),
                      funcs['f_%.1f_%.1f_mom'%(mom,theta)].GetParameter(2)] ) 

    

    gPad.Update()
    canv2[file].SaveAs('res_%.1f_%.1f.pdf'%(mom,theta))

print data
outfile = open('pickle.dat','w')
pickle.dump(data,outfile)
outfile.close()

gApplication.Run()

