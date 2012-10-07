#!/bin/env python

import sys, os, errno, re, optparse, subprocess, datetime
from os.path import expandvars
import ROOT

def main():
    # Parse command-line arguments
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input')
    parser.add_option('-o', '--output')
    parser.add_option('-s', '--script', default='applySVFit.sh')
    parser.add_option('-l', '--logdir', default='batch_logs')
    parser.add_option('-t', '--stepsize', type=int, default=500)
    parser.add_option('-z', '--zeroonly', action='store_true')
    options, args = parser.parse_args()

    if not options.input or not options.output:
        print "--input and --output options are required"
        print "Exiting"
        sys.exit(1)

    inputdir = os.path.abspath(options.input)
    outputdir = os.path.abspath(options.output)

    submitSVFit(options.script, inputdir, outputdir, options.logdir, options.stepsize, options.zeroonly)


def submitSVFit(script, inputdir, outputdir, logdir, stepsize, zeroonly = False):
    # Check CMSSW environment
    if 'CMSSW_BASE' not in os.environ:
        print "CMSSW_BASE is not set. Please set the CMSSW environment."
        print "Aborting..."
        sys.exit(1)

    # Make log directory
    label = script.rsplit('.', 1)[0]
    datestr = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    logsubdir = '%s/%s-%s' % (logdir, label, datestr)

    try:
        os.makedirs(logsubdir)
    except OSError as e:
        if e.errno != errno.EEXIST: raise

    files = os.listdir(inputdir)
    histFiles = filter(lambda x: re.search(r'.*\.root', x), files)

    for filename in histFiles:
        inputfile = os.path.join(inputdir, filename)
        submitFile(script, inputfile, outputdir, logsubdir, stepsize, zeroonly)


def submitFile(script, inputfile, outputdir, logsubdir, stepsize, zeroonly = False):
    # Get number of events
    infile = ROOT.TFile.Open(inputfile)
    ntEvt = infile.Get('ntEvt')
    ntLooseEvt = infile.Get('ntLooseEvt')
    nevents = max(ntEvt.GetEntries(), ntLooseEvt.GetEntries())
    infile.Close()

    # File names
    infilename = os.path.basename(inputfile)
    name, ext = os.path.splitext(infilename)
    os.environ['INPUT'] = inputfile

    if re.search('etau', infilename):
        os.environ['LEPTYPE1'] = 'electron'
    else:
        os.environ['LEPTYPE1'] = 'muon'

    os.environ['LEPTYPE2'] = 'tau'
    os.environ['EVENTS'] = str(stepsize)
    os.environ['DEST'] = outputdir
    os.environ['ZEROONLY'] = 'true' if zeroonly else 'false'

    # Submit jobs for range of events
    for i in xrange(0, nevents, stepsize):
        jobnum = int(i/stepsize)
        outputname = name + ('_%04i' % jobnum) + ext

        os.environ['SKIPEVENTS'] = str(i)
        os.environ['OUTPUT'] = outputname

        label = 'applySVFit_%s' % os.path.splitext(outputname)[0]
        cmd = 'bsub -r -q 8nh -J %(LABEL)s -oo %(LOGDIR)s/%(LABEL)s.log %(SCRIPT)s' % {'SCRIPT': script, 'LOGDIR': logsubdir, 'LABEL': label}
        print cmd
        subprocess.call(cmd, shell=True)


if __name__ == '__main__': main()
