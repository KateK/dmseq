#!/bin/bash

#PBS -N consolidate_summaries
#PBS -q generic
#PBS -l nodes=1:ppn=12
#PBS -o /home3/struck/scripts/log/
#PBS -e /home3/struck/scripts/log/
#PBS -j oe 

cd /home3/struck/Analysis_Projects/DMseq/bin/

# AS_event_types="ALE AFE nonUTRevents.multi"
AS_event_types="polyA"

groups="tibialis quadricep heart allControls DMseq"
# group=tibialis
# group=quadricep
# group=allControls
# group=heart
# group=DMseq

for group in $groups; do
    dir=/home3/struck/MISO/single_end_mode/summaries/$group
    for event_type in $AS_event_types; do
	summary_dir=$dir/$event_type
	header=TRUE

	file_type=miso_summary
	output=$dir/${group}_${event_type}_consolidatedSummaries.txt
	./createConsolidatedMISOSummaries_foreach.R $summary_dir $file_type $header $output

	# file_type=miso_bf
	# output=$dir/${group}_${event_type}_consolidatedBayesFactors.txt
	# ./createConsolidatedMISOSummaries_foreach.R $summary_dir $file_type $header $output
    done
done
