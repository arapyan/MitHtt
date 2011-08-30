out_fname = "%s/limit-inputs-mvis.root" % self.directory
out_file = TFile(out_fname,"RECREATE")

prefix = "emu"

for category in self.category_list:
    plot_name = "vismass_%s" % category.name
    dir_name = category.name 
    out_file.mkdir(dir_name)
    out_file.cd(dir_name)
    plot = Plot()
    for a_plot in category.plots:
        if a_plot.name=="vismass":
            plot=a_plot

    for group in plot.groups:
        ps = plot.group_styles[group.name]
        # need to modify this for 2D histograms
        curr_th = ROOT.TH1F()
        for sample in group.samples:
            last_sample = group.samples[-1]
            curr_th = self.container_dict[plot_name][sample.name]
            if "merge" in ps.mode:
                if sample is group.samples[0]:
                    curr_th = self.container_dict[plot_name][sample.name].Clone() # need the cloning to ensure original histograms are preserved
                else:
                    curr_th.Add(self.container_dict[plot_name][last_sample.name])
                if sample is not last_sample:
                    continue
            curr_th.SetName(sample.name)
            curr_th.Write()
    out_file.cd()

out_file.Write()
out_file.Close()

