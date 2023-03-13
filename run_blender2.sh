#!/bin/bash
# Stacia Wyman 22 July 2019
# Bash script to run BLENDER

# sh run_blender2.sh  <path to reference genome> <path to IP bam> <path to control bam> <guide sequence> <output directory>

REF=$1; shift
IP=$1; shift
CTRL=$1; shift
GUIDE=$1; shift
OUTDIR=$1; shift
OPTS="$@"

if [ ! -e $IP ]
then
    echo "$IP does not exist"
    exit
fi
if [ ! -e $CTRL ] 
then
    echo "$CTRL does not exist"
    exit
fi
    
if [ ! -d $OUTDIR ] 
then
   mkdir $OUTDIR
fi
command="blender2.py -f $IP -c $CTRL -g $GUIDE -r $REF $OPTS"
echo "Running $command"
python $command > $OUTDIR/unfiltered_blender_hits.txt

# For pooled samples, comment out the next line and uncomment the filter_pool.pl line
perl filter.pl $OUTDIR/unfiltered_blender_hits.txt > $OUTDIR/filtered_blender_hits.txt
#perl filter_pool.pl $OUTDIR/unfiltered_blender_hits.txt > $OUTDIR/filtered_pooled_blender_hits.txt

# Add PAM to guide for visualization
GUIDE+="NRG"
python draw_blender_fig.py $OUTDIR/filtered_blender_hits.txt $OUTDIR/blender_hits $GUIDE 
echo "Done"