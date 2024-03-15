OUTDIR=/icarus/data/users/gputnam/DMCP2023G/mc-FB
DF=evt

for dir in "$@"
do
  arrPath=(${dir//// })
  sPath=${arrPath[@]:7}
  savename=${sPath// /_}
  find $dir -name "*.flat.caf.root" | grep -v Blind | grep -v Prescaled > ${savename}.list
  echo Making dataframe  $OUTDIR/${savename}_${DF}.df with `wc -l < ${savename}.list` files
  python rundf.py configs/${DF}df.py $OUTDIR/${savename}_${DF}.df ${savename}.list
  rm ${savename}.list
done
