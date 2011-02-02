# .bashrc


# Source global definitions
if [ -f /etc/bashrc ]; then
	        . /etc/bashrc
fi

export PATH="/home/cmsprod/bin:$HOME/bin:${PATH}"
alias cms017='export MIT_VERS=017;export MIT_TAG=Mit_017;c=$HOME/cms;mkdir -p $c/cmssw/$MIT_VERS;cd $c/cmssw/$MIT_VERS;source $c/INIT 3_9_5_patch1; cd $c/root'
alias cms018='export MIT_VERS=018;export MIT_TAG=Mit_018;c=$HOME/cms;mkdir -p $c/cmssw/$MIT_VERS;cd $c/cmssw/$MIT_VERS;source $c/INIT 3_9_7; source $CMSSW_BASE/src/MitAna/bin/mitana.sh; export src=$CMSSW_BASE/src; cd $c/root'
alias mithtt='source $CMSSW_BASE/src/MitAna/bin/mitana.sh;source $CMSSW_BASE/src/MitHtt/bin/mithtt.sh'
alias install='installBambu.sh'
alias compile='cd $CMSSW_BASE/src; scram build MitPlots MitHtt; cd -'
alias edsrc='emacs $src'
