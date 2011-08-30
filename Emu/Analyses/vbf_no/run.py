#!/usr/bin/env python
import os
from optparse import OptionParser
from MitHtt.Emu.RootUtils import Analysis

parser=OptionParser()
parser.add_option("--p ",dest="do_cplots",default=False,action="store_true",help="make the nice plots?")
parser.add_option("--c ",dest="config_name",default="analysis.config",type="string",help="specify the config file name")
parser.add_option("--n ",dest="ntuple_dir",default="DEFAULT_NO",type="string",help="specify the ntuple_dir name")
(opts,args)=parser.parse_args()

config_name="analysis.config"

mod_path=os.path.dirname(os.path.abspath(__file__))
ana=Analysis(opts.config_name,opts.ntuple_dir)

ana.run(opts.do_cplots)

