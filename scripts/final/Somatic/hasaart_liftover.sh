mamba activate kent-tools

hg19ToHg38chain=/oak/stanford/groups/smontgom/amarder/bin/chain/hg19ToHg38.over.chain.gz
dir=/oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep
inDir=$dir/hg19
outDir=$dir/hg38

rm -r $inDir
rm -r $outDir

mkdir $inDir
mkdir $outDir

for origFile in $dir/*bed;
do

echo $origFile

fileName=${origFile#"$dir"}
inFile=${inDir}${fileName}
outFile=${outDir}${fileName}
unmappFile=${outFile#".bed"}.unmapp.bed

awk '{ $1 = "chr" $1; print }' $origFile | tail -n+2 > $inFile

# echo $outFile
# echo $unmappFile

echo "Running liftOver..."
liftOver $inFile $hg19ToHg38chain $outFile ${outFile%".bed"}.unmapp.bed

done

mamba deactivate 