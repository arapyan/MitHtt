if [ -z $CMSSW_BASE ]
then
  echo ""
  echo " Setting up MitHtt failed! (\$CMSSW_BASE = is empty)."
  echo ""
else
  export MIT_HTT_DIR="$CMSSW_BASE/src/MitHtt"
  export PATH="$MIT_HTT_DIR/bin:${PATH}"
  export PYTHONPATH="$MIT_HTT_DIR/python:${PYTHONPATH}"
fi
