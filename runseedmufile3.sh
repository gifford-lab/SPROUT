#!/bin/bash

java -Xmx4g -cp /cluster/ccr/sprout_seed/sproutseed.jar edu.mit.csail.cgs.reeder.sproutseed.SproutSeed --genomefile "$genome" --rho $rho --alpha $alpha --beta $beta --a $a --b $b --readfile "$readfile" --dumpfile "$dumpfile" --outfile "$outfile" --genericdirectory "$genericdirectory" --stage $stage --maxiters $maxiters --stage2file "$stage2file" --directory "$directory" --eventout "$eventout" --interactionout "$interactionout" --region "$region"
