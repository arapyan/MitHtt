#-------------------------------------------------------------------------------------------
# select probes
#-------------------------------------------------------------------------------------------

#
# Electron Trigger Efficiencies
#
#root -l -q selectMuEGElectronLegTrigEffTP.C+\(\"data_smu.conf\",\"data_MuEGEleLeg\",0\)
#root -l -q selectMuEGElectronLegTrigEffTP.C+\(\"s12-tt.conf\",\"s12_tt_MuEGEleLeg\",0\)
#root -l -q selectMuEGElectronLegTrigEffTP.C+\(\"s12-ztt.conf\",\"s12_ztt_MuEGEleLeg\",0,1\)

#
# Electron Selection Efficiencies
#
#root -l -q selectEleIDIsoEffTP.C+\(\"s12-zll.conf\",\"s12_ele_IDIsoEffTP\",1\)
#root -l -q selectEleIDIsoEffTP.C+\(\"data_el.conf\",\"data_ele_IDIsoEffTP\"\)



#-------------------------------------------------------------------------------------------
# compute efficiencies
#-------------------------------------------------------------------------------------------

#
# Electron Trigger Efficiencies
#
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"s12_ztt_MuEGEleLeg/probes.root\",\"s12_ztt_MuEGEleLeg/pteta_76_106\",\"all\",1,0,\"\"\,\"PUWeights_full2012.root\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"s12_tt_MuEGEleLeg/probes.root\",\"s12_tt_MuEGEleLeg/pteta_76_106\",\"all\",1,0,\"\"\,\"PUWeights_full2012.root\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"data_MuEGEleLeg/probes.root\",\"data_MuEGEleLeg/var_76_106\",\"all\",1,0,\"\"\,\"\"\)

#
# Electron Selection Efficiencies
#
root -l -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"s12_ele_IDIsoEffTP/probes.root\",\"s12_ele_IDIsoEffTP/pteta_76_106\",\"all\",1,0,\"\",\"PUWeights_full2012.root\"\)
root -l -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"data_ele_IDIsoEffTP/probes.root\",\"data_ele_IDIsoEffTP/pteta_76_106\",\"all\",1,0,\"s12_eleIDIsoEffTP/probes.root\"\)



#-------------------------------------------------------------------------------------------
# compute scale factors
#-------------------------------------------------------------------------------------------

#
# Electron Trigger Efficiencies
#
#root -l -q makeEfficiencyScaleFactors.C+\(\"data_MuEGEleLeg/pteta_76_106/eff.root\",\"s12_ztt_MuEGEleLeg/pteta_76_106/eff.root\",\"data_MuEGEleLeg/pteta_76_106/\",\"h2_results_muegele_trig\",\"hztt\"\)

#
# Electron Selection Efficiencies
#
#root -l -q makeEfficiencyScaleFactors.C+\(\"data_ele_IDIsoEffTP/pteta_76_106/eff.root\",\"s12_eleIDIsoEffTP/pteta_76_106/eff.root\",\"data_ele_IDIsoEffTP/pteta_76_106/\",\"h2_results_ele_selection\",\"htt\"\)


rm *.so *.d
