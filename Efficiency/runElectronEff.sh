#-------------------------------------------------------------------------------------------
# select probes
#-------------------------------------------------------------------------------------------

#
# Electron Trigger Efficiencies
#
#root -l -q selectMuEGElectronLegTrigEffTP.C+\(\"data_smu.conf\",\"data_MuEGEleLeg\",0\)
#root -l -q selectMuEGElectronLegTrigEffTP.C+\(\"s11-tt.conf\",\"s11_tt_MuEGEleLeg\",0\)
#root -l -q selectMuEGElectronLegTrigEffTP.C+\(\"s11-ztt.conf\",\"s11_ztt_MuEGEleLeg\",0,1\)

#
# Electron Selection Efficiencies
#
#root -l -q selectEleMVAWPEffTP.C+\(\"s11-zee.conf\",\"s11_zee_ele_MVAWPEffTP\",1\)
#root -l -q selectEleMVAWPEffTP.C+\(\"data_el.conf\",\"data_ele_MVAWPEffTP\"\)
#root -l -q selectEleMVAIDOnlyEffTP.C+\(\"s11-zee.conf\",\"s11_zee_ele_MVAIDEffTP\",1\)
#root -l -q selectEleMVAIDOnlyEffTP.C+\(\"data_el.conf\",\"data_eleMVAIDEffTP\"\)



#-------------------------------------------------------------------------------------------
# compute efficiencies
#-------------------------------------------------------------------------------------------

#
# Electron Trigger Efficiencies
#
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"s11_ztt_MuEGEleLeg/probes.root\",\"s11_ztt_MuEGEleLeg/basic_76_106\",\"all\",1,0,\"\"\,\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"s11_tt_MuEGEleLeg/probes.root\",\"s11_tt_MuEGEleLeg/basic_76_106\",\"all\",1,0,\"\"\,\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"data_MuEGEleLeg/probes.root\",\"data_MuEGEleLeg/var_76_106\",\"all\",1,0,\"\"\,\"\"\)

#
# Electron Selection Efficiencies
#
#root -l -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"s11_zee_ele_MVAWPEffTP/probes.root\",\"s11_zee_ele_MVAWPEffTP/basic_76_106\",\"all\",1,0,\"\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"data_eleMVAWPEffTP/probes.root\",\"data_eleMVAWPEffTP/basic_76_106\",\"all\",1,0,\"s11_zee_eleMVAWPEffTP/probes.root\"\)
#root -l -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"s11_zee_ele_MVAIDEffTP/probes.root\",\"s11_zee_ele_MVAIDEffTP/basic_76_106\",\"all\",1,0,\"\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"data_eleMVAIDEffTP/probes.root\",\"data_eleMVAIDEffTP/basic_76_106\",\"all\",1,0,\"s11_zee_eleMVAIDEffTP/probes.root\"\)



#-------------------------------------------------------------------------------------------
# compute scale factors
#-------------------------------------------------------------------------------------------

#
# Electron Trigger Efficiencies
#
#root -l -q makeEfficiencyScaleFactors.C+\(\"data_MuEGEleLeg/basic_76_106/eff.root\",\"s11_ztt_MuEGEleLeg/basic_76_106/eff.root\",\"data_MuEGEleLeg/basic_76_106/\",\"h2_results_muegele_trig\",\"hztt\"\)

#
# Electron Selection Efficiencies
#
#root -l -q makeEfficiencyScaleFactors.C+\(\"data_ele_MVAWPEffTP/basic_76_106/eff.root\",\"s11_zee_eleMVAWPEffTP/basic_76_106/eff.root\",\"data_ele_MVAWPEffTP/basic_76_106/\",\"h2_results_ele_selection\",\"htt\"\)
#root -l -q makeEfficiencyScaleFactors.C+\(\"data_eleMVAIDEffTP/basic_76_106/eff.root\",\"s11_zee_eleMVAIDEffTP/basic_76_106/eff.root\",\"data_eleMVAIDEffTP/basic_76_106/\",\"h2_results_ele_selection\",\"htt\"\)


rm *.so *.d
