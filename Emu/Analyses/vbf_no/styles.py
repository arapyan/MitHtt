# Standard plotting styles.
from MitHtt.Emu.PlotUtils import PlotStyle
import ROOT

def main():
	ROOT.gStyle.SetLineStyleString(11,"20 10")
	signal_stackscale = PlotStyle("signal_stack10")
	signal_stackscale.fillcolor = 1
	signal_stackscale.fillstyle = 0
	signal_stackscale.linecolor = 603
	signal_stackscale.linestyle = 11
	signal_stackscale.linewidth = 3
	signal_stackscale.mode = "linestack"
	signal_stackscale.scale = 10

	signal_overscale = PlotStyle("signal_overlay10")
	signal_overscale.fillcolor = 1
	signal_overscale.fillstyle = 0
	signal_overscale.linecolor = 603
	signal_overscale.linestyle = 11
	signal_overscale.linewidth = 3
	signal_overscale.mode = "mergeline"
	signal_overscale.scale = 10

        merge_fill = PlotStyle("merge_stack")
        merge_fill.fillcolor =1
        merge_fill.fillstyle =1001
        merge_fill.linecolor = 1
        merge_fill.linestyle = 1
        merge_fill.linewidth = 3
        merge_fill.mode = "mergestack"
        merge_fill.scale = 1

	def_fill = PlotStyle("default_stack")
	def_fill.fillcolor =1
	def_fill.fillstyle =1001
	def_fill.linecolor = 1
	def_fill.linestyle = 1
	def_fill.linewidth = 3
	def_fill.mode = "stack"
	def_fill.scale = 1

	data_style = PlotStyle("data_style")
	data_style.mode = "point"
	data_style.scale = 1
	data_style.fillstyle = 0
	data_style.linestyle = 1

        styles_dict={}
        styles_dict[signal_stackscale.name]=signal_stackscale
        styles_dict[signal_overscale.name]=signal_overscale
        styles_dict[def_fill.name]=def_fill
        styles_dict[data_style.name]=data_style
        styles_dict[merge_fill.name] = merge_fill
        return styles_dict

if __name__=='__main__': main()
