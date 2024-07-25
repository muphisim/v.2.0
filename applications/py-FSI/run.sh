outdir=test-3deg
rm -rf $outdir
mkdir $outdir
rm -rf initialGeo
python3 makeInitialBoundary.py
mv initialGeo $outdir
cp -rf referenceCase  $outdir
python3 runFSI.py -case $outdir  \
                  -aoa 3 -parallel 0 \
                  -young 5e5 -nu 0.3 \
                  -omega 0.5 -maxIter 20 | tee $outdir/log.txt


cp plotGeo.py $outdir
cd $outdir
python3 plotGeo.py 
