 #!/bin/bash 
 COUNTER=22
 while [  $COUNTER -lt 22 ]; do   
     qsub -F "real "$COUNTER -N "real_chr"$COUNTER  MLCB_run.pbs
     let COUNTER=COUNTER+1 
 done