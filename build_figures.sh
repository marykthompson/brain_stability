#Run notebooks and Rscripts to create the figures
#$ bash build_figures.sh TestFigs/
# run from the conda env pdfcrop
dir="$1"

# Make directory for the cropped pdfs
rm -rf $dir/cropped
rm -rf $dir/temp
rm -rf $dir/uncropped
mkdir $dir/cropped
mkdir $dir/temp
mkdir $dir/uncropped

# Build composites with xelatex
for file in $dir/*.tex; do xelatex -output-directory=$dir -aux-directory=$dir/temp $file; done
for file in $dir/*.aux; do mv $file $dir/temp; done
for file in $dir/*.log; do mv $file $dir/temp; done

for file in $dir/*.pdf
do 
    newname=$dir/cropped/${file##*/}
    newname2=${newname%.*}_cropped.pdf
    pdf-crop-margins -v -s -u -p 3 -o $newname2 $file
done

for file in $dir/*.pdf; do mv $file $dir/uncropped; done