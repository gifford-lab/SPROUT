#!/bin/bash

java -Xmx2g -cp /cluster/ccr/sprout_seed/sproutseed.jar edu.mit.csail.cgs.reeder.sproutseed.SproutSeed --genomefile "$genome" --rho $rho --alpha $alpha --beta $beta --a $a --b $b --readfile "$readfile" --dumpfile "$dumpfile" --outfile "$outfile" --directory "$directory" --stage $stage --maxiters $maxiters --mutaufile "$mutaufile" --nodump
