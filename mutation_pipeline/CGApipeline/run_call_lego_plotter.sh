#!/bin/bash

args=
while [ $# -gt 0 ]; do
	token=`echo "$1" | sed 's/ /\\\\ /g'`   # Add blackslash before each blank
	args="${args} ${token}" 
	shift
done
echo "args is $args" ; 
source /matlab_source_file_2012a.sh 
	Xvnc :$$ -depth 16&
	XVNC_PID=$!

	export DISPLAY=:$$
set -x
  /usr/local/bin/call_lego_plotter $args
PLOTTER_EXIT_CODE=$? ;
set +x
	kill $XVNC_PID
if [ "0" -ne "$PLOTTER_EXIT_CODE" ] ; 
then
	echo "WARNING : PLOTTER EXIT CODE WAS $PLOTTER_EXIT_CODE" ; 
fi ;

exit $PLOTTER_EXIT_CODE
