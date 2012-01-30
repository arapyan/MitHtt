#-------------------------------------------------------------------------------------------
# select probes
#-------------------------------------------------------------------------------------------

#
# Muon Trigger Efficiencies
#
#root -l -q selectMuEGMuonLegTrigEffTP.C+\(\"data_sel.conf\",\"data_MuEGMuonLeg\",0\)
#root -l -q selectMuEGMuonLegTrigEffTP.C+\(\"s11-tt.conf\",\"s11_tt_MuEGMuonLeg\",0\)
#root -l -q selectMuEGMuonLegTrigEffTP.C+\(\"s11-ztt.conf\",\"s11_ztt_MuEGMuonLeg\",0,1\)

#
# Muon Selection Efficiencies
#
#root -l -q selectMuonWPEffTP.C+\(\"data_mu.conf\",\"data_muonWPEffTP\"\)
#root -l -q selectMuonWPEffTP.C+\(\"s11-zmm.conf\",\"s11_zmm_muonWPEffTP\",1\)



#-------------------------------------------------------------------------------------------
# compute efficiencies
#-------------------------------------------------------------------------------------------

#
# Muon Trigger Efficiencies
#
#root -l -b -q plotEff.C+\(\"muMuEGMuonTrig.bins\",0,0,0,0,\"s11_ztt_MuEGMuonLeg/probes.root\",\"s11_ztt_MuEGMuonLeg/basic_76_106\",\"all\",1,0,\"\"\,\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"muMuEGMuonTrig.bins\",0,0,0,0,\"s11_tt_MuEGMuonLeg/probes.root\",\"s11_tt_MuEGMuonLeg/basic_76_106\",\"all\",1,0,\"\"\,\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"muMuEGMuonTrig.bins\",0,0,0,0,\"data_MuEGMuonLeg/probes.root\",\"data_MuEGMuonLeg/basic_76_106\",\"all\",1,0,\"\"\,\"\"\)

#
# Muon Selection Efficiencies
#
#root -l -q plotEff.C+\(\"mu0.bins\",0,0,0,0,\"s11_zmm_muonWPEffTP/probes.root\",\"s11_zmm_muonWPEffTP/basic_76_106\",\"all\",1,0,\"\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -q plotEff.C+\(\"mu0.bins\",2,1,2,1,\"data_muonWPEffTP/probes.root\",\"data_muonWPEffTP/basic_76_106\",\"all\",1,0,\"s11_zmm_muonWPEffTP/probes.root\"\)



#-------------------------------------------------------------------------------------------
# compute scale factors
#-------------------------------------------------------------------------------------------

#
# Muon Trigger Efficiencies
#
#root -l -q makeEfficiencyScaleFactors.C+\(\"data_MuEGMuonLeg/basic_76_106/eff.root\",\"s11_ztt_MuEGMuonLeg/basic_76_106/eff.root\",\"data_MuEGMuonLeg/basic_76_106/\",\"h2_results_muegmuon_trig\",\"htt\"\)

#
# Muon Selection Efficiencies
#
#root -l -q makeEfficiencyScaleFactors.C+\(\"data_muonWPEffTP/basic_76_106/eff.root\",\"s11_zmm_muonWPEffTP/basic_76_106/eff.root\",\"data_muonWPEffTP/basic_76_106/\",\"h2_results_mu_selection\",\"htt\"\)


rm *.so *.d
