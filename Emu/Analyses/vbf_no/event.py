
vismass=data.data.mass
svmass=data.data.svfmass

vbfcuts= ((data.data.mjj > mjjMin) and (data.data.jeta1*data.data.jeta2 < 0) and (abs(data.data.jeta1-data.data.jeta2) > dEtaMin))

uni_weight = lumi*data.data.weight

lepton_pt = []
lepton_pt.append(data.data.lpt1)
lepton_pt.append(data.data.lpt2)

from ROOT import TVector3
lv1 = TVector3()
lv2 = TVector3()
met = TVector3()
lv1.SetPtEtaPhi(data.data.lpt1,data.data.leta1,data.data.lphi1)
lv2.SetPtEtaPhi(data.data.lpt2,data.data.leta2,data.data.lphi2)
met.SetPtEtaPhi(data.data.met,0,data.data.metphi)
higgs_vec = lv1 + lv2 + met
higgs_pt = higgs_vec.Pt()

high_jmt = data.data.jpt1>100 and data.data.met>75
high_ht = higgs_pt>100
one_jet = data.data.njets ==1
two_jets = data.data.njets == 2
no_bjets = data.data.nbjets == 0
lt_two_jets = data.data.njets<2
lte_two_jets = data.data.njets<3

vbf_selection = vbfcuts and two_jets
novbf_selection = lte_two_jets and (not vbf_selection)

