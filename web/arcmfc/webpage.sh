#!/bin/bash

YEAR=`date +%Y`
MONTH=`date +%m`

# Create graph plots based on validation results:

# image files
Imgpath=/lustre/storeB/project/fou/om/ARCMFC/S3a/ValidationFigures/${YEAR}/${MONTH}/
Figtempl=ARCMFC_fig_val_ts
rmseImg=${Imgpath}${Figtempl}_rmsd_lt012h_${YEAR}${MONTH}.png
madImg=${Imgpath}${Figtempl}_mad_lt012h_${YEAR}${MONTH}.png
biasImg=${Imgpath}${Figtempl}_bias_lt012h_${YEAR}${MONTH}.png
corrImg=${Imgpath}${Figtempl}_corr_lt012h_${YEAR}${MONTH}.png
SIImg=${Imgpath}${Figtempl}_SI_lt012h_${YEAR}${MONTH}.png
novImg=${Imgpath}${Figtempl}_nov_lt012h_${YEAR}${MONTH}.png
scatterImg=${Imgpath}ARCMFC_fig_val_scatter_lt012h_${YEAR}${MONTH}.png

Webpath=/lustre/storeB/project/fou/om/ARCMFC/S3a/WebPage/
mkdir -p ${Webpath}${YEAR}/${MONTH}/
cp ${Imgpath}* ${Webpath}${YEAR}/${MONTH}/.

# Create html document:
Imgpath=${YEAR}/${MONTH}/
Figtempl=ARCMFC_fig_val_ts
rmseImg=${Imgpath}${Figtempl}_rmsd_lt012h_${YEAR}${MONTH}.png
madImg=${Imgpath}${Figtempl}_mad_lt012h_${YEAR}${MONTH}.png
biasImg=${Imgpath}${Figtempl}_bias_lt012h_${YEAR}${MONTH}.png
corrImg=${Imgpath}${Figtempl}_corr_lt012h_${YEAR}${MONTH}.png
SIImg=${Imgpath}${Figtempl}_SI_lt012h_${YEAR}${MONTH}.png
novImg=${Imgpath}${Figtempl}_nov_lt012h_${YEAR}${MONTH}.png
scatterImg=${Imgpath}ARCMFC_fig_val_scatter_lt012h_${YEAR}${MONTH}.png

htmlFile=${Webpath}index.html

bName='Bulletin date: '`date +%Y-%m-%d`
echo $bName

cp head1.html ${htmlFile}
echo '<title>'${bName}'</title>'       >> ${htmlFile}
cat head2.html                         >> ${htmlFile}

echo '<h2>Results for '${bName}'</h2>' >> ${htmlFile}

echo '<table><tr><td>'                 >> ${htmlFile}
#echo 'Time series of root mean square error' >> ${htmlFile}
echo '<img src="'${rmseImg}'">'        >> ${htmlFile}
#echo '</td><td>'                       >> ${htmlFile}
#echo '<br>'                            >> ${htmlFile}
#echo '<img src="'${madImg}'">'         >> ${htmlFile}
echo '<br>'                            >> ${htmlFile}
echo '<img src="'${biasImg}'">'        >> ${htmlFile}
echo '<br>'                            >> ${htmlFile}
echo '<img src="'${corrImg}'">'        >> ${htmlFile}
echo '<br>'                            >> ${htmlFile}
echo '<img src="'${SIImg}'">'          >> ${htmlFile}
echo '<br>'                            >> ${htmlFile}
echo '<img src="'${novImg}'">'         >> ${htmlFile}
echo '<img src="'${scatterImg}'">'     >> ${htmlFile}
echo '</td></tr></table>'              >> ${htmlFile}

echo '<p><small>'                      >> ${htmlFile}
echo 'Generated on '`date`             >> ${htmlFile}
echo '</small></p>'                    >> ${htmlFile}

cat tail.html                          >> ${htmlFile}

exit
