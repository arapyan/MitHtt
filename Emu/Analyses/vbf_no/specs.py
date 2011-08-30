# Defines some plots that are useful to us in the Htt search
import ROOT

from ROOT import gROOT
from MitHtt.Emu.PlotUtils import PlotStyle,Plot,Box

def main():

    CMSSW_BASE=ROOT.gSystem.Getenv("CMSSW_BASE")
    gROOT.LoadMacro(CMSSW_BASE+"/src/MitHtt/macros/HttStyles.cc+")

    ROOT.SetStyle()
    print "I set the style"
    cms_pave_texts = []
    cms_pave_texts.append(ROOT.cmslab("#tau_{e}#tau_{#mu}", 0.45, 0.75))
    cms_pave_texts.append(ROOT.cmslumi("#tau_{e}#tau_{#mu}",0.45,0.75))
    cms_pave_texts.append(ROOT.cmschan("#tau_{e}#tau_{#mu}",0.45,0.75))

    ru_legend = Box("ru_legend",0.45,0.45,0.9,0.75)

    plot_dict = {}
    vismass_plot = Plot(type="TH1D",name="vismass",xvar="vismass",xmin=0,xmax=500,xbins=25,xunits="GeV",xtitle="m_{vis}",ytitle="Events",ylog=False,legend_box=ru_legend,textboxes=cms_pave_texts)
    plot_dict[vismass_plot.name]=vismass_plot
    svmass_plot = Plot(type="TH1D",name="svmass",xvar="svmass",xmin=0,xmax=500,xbins=25,xunits="GeV",xtitle="m_{vis}",ytitle="Events",ylog=False,legend_box=ru_legend,textboxes=cms_pave_texts)
    plot_dict[svmass_plot.name]=svmass_plot
    leppt_plot = Plot(type="TH1D",name="leppt",xvar="lepton_pt",xmin=0,xmax=100,xbins=25,xunits="GeV",xtitle="m_{vis}",ytitle="Events",ylog=False,legend_box=ru_legend,textboxes=cms_pave_texts)
    plot_dict[leppt_plot.name]=leppt_plot
    return plot_dict

if __name__=='__main__': main()
