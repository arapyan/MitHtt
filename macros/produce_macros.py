from DatacardUtils import parse_dcard
# use parse_dcard to get a dictionary mapping sample name strings to fit weights

class Analysis:
    '''A class designed to insert the proper scale factors into Roger's plotting macros'''
    def __init__(self,process_weight,template_fname,output_fname):
         '''Takes a dictionary (mapping strings representing samples) of fit weights and inserts these into the template macro at template_fname. Output is written to output_fname'''
         self.process_weight = process_weight
         self.template_fname = template_fname
         self.output_fname = output_fname

    def run(self):
         '''Inserts the weights into the macros'''
         input_file = open(self.template_fname,'r')
         output_file = open(self.output_fname,'w')
         
         curr_name = ""
         for line in input_file:
             move_on = False
             word_arr=line.split("\n")
             for process_name in self.process_weight.keys():
                 cand_str = "$ %s" % process_name
                 output_cand = ""
                 if line.strip().startswith(cand_str):
                     print word_arr[0]
                     curr_name = process_name
                     print_me = '''std::cout<< "scaling by %f"<<std::endl;''' % self.process_weight[curr_name]
                     out_line = print_me+"hin->Scale(%f); break; \n" % self.process_weight[curr_name]
                     move_on = True
                     output_file.write(out_line)
                     print out_line
             if not move_on:
                 output_file.write(line)
                 

# specify the fit results text files for sm and mssm cases                    
mssm_results = "fitresults/2108Approval/mike_2208_mssm.txt"
sm_results = "fitresults/2108Approval/mike_2208_sm.txt"    
tmp_dir="templates/"
etau_b_ana = Analysis(parse_dcard("htautau/mssm/eleTau_B_mA120.txt",mssm_results,"ANYBIN"),tmp_dir+"etauAfterFit_b_template.C","etauAfterFit_b.C")
etau_b_ana.run()

etau_nob_ana = Analysis(parse_dcard("htautau/mssm/eleTau_NoB_mA120.txt",mssm_results,"ANYBIN"),tmp_dir+"etauAfterFit_nob_template.C","etauAfterFit_nob.C")
etau_nob_ana.run()

emu_b_ana = Analysis(parse_dcard("htautau/mssm/eleMu_BTag_mA120.txt",mssm_results,"emu_b"),tmp_dir+"emuAfterFit_b_template.C","emuAfterFit_b.C")
emu_b_ana.run()

emu_nob_ana = Analysis(parse_dcard("htautau/mssm/eleMu_BTag_mA120.txt",mssm_results,"emu_nob"),tmp_dir+"emuAfterFit_nob_template.C","emuAfterFit_nob.C")
emu_nob_ana.run()

etau_vbf_ana = Analysis(parse_dcard("htautau/sm/eleTau_SM2_mH120.txt",sm_results,"ANYBIN"),tmp_dir+"etauAfterFit_vbf_template.C","etauAfterFit_vbf.C")
etau_vbf_ana.run()

etau_novbf_ana = Analysis(parse_dcard("htautau/sm/eleTau_SM0_mH120.txt",sm_results,"ANYBIN"),tmp_dir+"etauAfterFit_novbf_template.C","etauAfterFit_novbf.C")
etau_novbf_ana.run()

emu_vbf_ana = Analysis(parse_dcard("htautau/sm/eleMu_VBF_mH120.txt",sm_results,"emu_vbf"),tmp_dir+"emuAfterFit_vbf_template.C","emuAfterFit_vbf.C")
emu_vbf_ana.run()

emu_novbf_ana = Analysis(parse_dcard("htautau/sm/eleMu_VBF_mH120.txt",sm_results,"emu_novbf"),tmp_dir+"emuAfterFit_novbf_template.C","emuAfterFit_novbf.C")
emu_novbf_ana.run()


mutau_b_ana = Analysis(parse_dcard("htautau/mssm/muTau_B_mA120.txt",mssm_results,"ANYBIN"),tmp_dir+"mutauAfterFit_b_template.C","mutauAfterFit_b.C")
mutau_b_ana.run()

mutau_nob_ana = Analysis(parse_dcard("htautau/mssm/muTau_NoB_mA120.txt",mssm_results,"ANYBIN"),tmp_dir+"mutauAfterFit_nob_template.C","mutauAfterFit_nob.C")
mutau_nob_ana.run()

mutau_vbf_ana = Analysis(parse_dcard("htautau/sm/muTau_SM2_mH120.txt",sm_results,"ANYBIN"),tmp_dir+"mutauAfterFit_vbf_template.C","mutauAfterFit_vbf.C")
mutau_vbf_ana.run()

mutau_novbf_ana = Analysis(parse_dcard("htautau/sm/muTau_SM0_mH120.txt",sm_results,"ANYBIN"),tmp_dir+"mutauAfterFit_novbf_template.C","mutauAfterFit_novbf.C")
mutau_novbf_ana.run()


