#!/bin/bash
#
# sam files indexes sam files in current folder
# tested samtools v 1.5
#

FILES=`ls -1 *.sam`

for SamFile in $FILES
do
	i=$(echo $SamFile | sed "s/\..*$//")
	BamFile=$i"_sorted.bam"
	#printf "%s  %s\n" $i,$samfile
	# test wethere soerted.bam exists
	
	if [[ ! -f $BamFile ]]; then
		# sort sam files and convert to bam
		echo "samtools sort -m 4G -@ 4 -o  $BamFile $SamFile"
		samtools sort -m 4G -@ 4 -o  $BamFile $SamFile
		# index bam files
		echo "samtools index $BamFile"
		samtools index $BamFile
	else
		echo "File $BamFile exists!" 	
	fi
	
done
