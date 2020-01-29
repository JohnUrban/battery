#!/bin/bash


## NEED: BASE, SAM, REF, MISMATCH, 

echo BASE $BASE
echo SAM $SAM
echo REF $REF
echo MISMATCH $MISMATCH
echo P $P

##Example ONT: MISMATCH=0.3366179

## Note that in its current form, LAP is not entirely great for the BWA long read alignments
## I'd recommend just taking a look at the score if curious, but not using it for final judgements.
## Overall - it seems highly correlated with %aligned and/or MAPQ scores anyway...

##PROBS
 if [ -f ${BASE}.prob ]; then
   c=`cat ${BASE}.prob | wc -l`
   if [ $c -lt 10 ]; then ##it is just a ghost file, so make it for reals
     ## Note: -i $SAM is just a dummy/place holder here - since a SAM is given, it bypasses -i.
     calc_prob.py -p $P -a $REF -i $SAM -s $SAM -c $MISMATCH > ${BASE}.prob
   fi
 else ##DoesNotExist so create
  date; echo prob file DoesNotExist so create
  calc_prob.py -p $P -a $REF -i $SAM -s $SAM -c $MISMATCH > ${BASE}.prob
  date; echo prob file DoesNotExist so created it...
 fi
 
## LAPSCORE
 if [ -f ${BASE}.lapscore ]; then
   c=`cat ${BASE}.lapscore | wc -l`
   if [ $c -lt 1 ]; then ##it is just a ghost file, so make it for reals
      sum_prob.py -i ${BASE}.prob -t 1e-323 -d 3 > ${BASE}.lapscore.detailed
      tail -n 1 ${BASE}.lapscore.detailed | awk 'OFS="\t" {print $3,$4}' > ${BASE}.lapscore
   fi
 else ##DoesNotExist so create
   sum_prob.py -i ${BASE}.prob -t 1e-323 -d 3 > ${BASE}.lapscore.detailed
   tail -n 1 ${BASE}.lapscore.detailed | awk 'OFS="\t" {print $3,$4}' > ${BASE}.lapscore
 fi


if $CLEAN; then bash ${SCRIPTS}/lap.clean.sh ; fi
