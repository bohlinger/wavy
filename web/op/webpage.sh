#!/bin/bash

# setup
while getopts i:w:m:r:s:Y:M: flag
do
    case "${flag}" in
        i) i=${OPTARG};;
        w) w=${OPTARG};;
        m) m=${OPTARG};;
        r) r=${OPTARG};;
        s) s=${OPTARG};;
        Y) Y=${OPTARG};;
        M) M=${OPTARG};;
    esac
done

img=${i}
web=${w}
mod=${m}
sat=${s}

if [ -z ${Y+x} ] | [ -z ${m+x} ]; then
  echo "Date not supplied correctly, current month is chosen"
  YEAR=`date +%Y`
  MONTH=`date +%m`
else
  echo "Date set manually"
  YEAR=${Y}
  MONTH=${M}
fi

if [ -z ${r+x} ]; then
  echo "Region not supplied correctly, model is chosen as region"
  reg=${m}
else
  echo "Region set manually"
  reg=${r}
fi

# translate input
# path to validation figures
Imgpath=${img}/${mod}/satellites/altimetry/${sat}/ValidationFigures/${YEAR}/${MONTH}/

# Validation figures:
rmseImg=${Imgpath}${mod}_vs_${sat}_for_${reg}_fig_val_ts_rmsd_${YEAR}${MONTH}.png
madImg=${Imgpath}${mod}_vs_${sat}_for_${reg}_fig_val_ts_mad_${YEAR}${MONTH}.png
biasImg=${Imgpath}${mod}_vs_${sat}_for_${reg}_fig_val_ts_bias_${YEAR}${MONTH}.png
corrImg=${Imgpath}${mod}_vs_${sat}_for_${reg}_fig_val_ts_corr_${YEAR}${MONTH}.png
SIImg=${Imgpath}${mod}_vs_${sat}_for_${reg}_fig_val_ts_SI_${YEAR}${MONTH}.png
novImg=${Imgpath}${mod}_vs_${sat}_for_${reg}_fig_val_ts_nov_${YEAR}${MONTH}.png
scatterImg1=${Imgpath}${mod}_vs_${sat}_for_${reg}_fig_val_scatter_lt000h_${YEAR}${MONTH}.png
scatterImg2=${Imgpath}${mod}_vs_${sat}_for_${reg}_fig_val_scatter_lt024h_${YEAR}${MONTH}.png
scatterImg3=${Imgpath}${mod}_vs_${sat}_for_${reg}_fig_val_scatter_lt048h_${YEAR}${MONTH}.png

# Create html document:
htmlFile=${web}/index.html

bName='Bulletin date: '`date +%Y-%m-%d`
echo $bName

cat ${web}/head1.html               >  ${htmlFile}
echo '<title>'${bName}'</title>'       >> ${htmlFile}
cat ${web}/head2.html               >> ${htmlFile}

echo '<h2>Results for '${bName}'</h2>' >> ${htmlFile}
echo '...some additional info...' >> ${htmlFile}
echo '<br>'                            >> ${htmlFile}

echo '<img src="'${rmseImg}'">'        >> ${htmlFile}
echo '<br>'                            >> ${htmlFile}
echo '<img src="'${biasImg}'">'        >> ${htmlFile}
echo '<br>'                            >> ${htmlFile}
echo '<img src="'${corrImg}'">'        >> ${htmlFile}
echo '<br>'                            >> ${htmlFile}
echo '<img src="'${SIImg}'">'          >> ${htmlFile}
echo '<br>'                            >> ${htmlFile}
echo '<img src="'${novImg}'">'         >> ${htmlFile}
echo '<br>'                            >> ${htmlFile}
echo '<img src="'${scatterImg1}'">'    >> ${htmlFile}
echo '<img src="'${scatterImg2}'">'    >> ${htmlFile}
echo '<img src="'${scatterImg3}'">'    >> ${htmlFile}

echo '<p><small>'                      >> ${htmlFile}
echo 'Generated on '`date`             >> ${htmlFile}
echo '</small></p>'                    >> ${htmlFile}

cat ${web}/tail.html                >> ${htmlFile}

exit
