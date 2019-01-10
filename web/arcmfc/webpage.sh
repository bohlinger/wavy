#!/bin/bash


# Create graph plots based on validation results:

# image files
VALpath=/lustre/storeB/project/fou/om/ARCMFC/S3a/ValidationFigures
rmseImg=${VALpath}/ARCMFC_fig_val_ts_rmsd_lt012h_201901.png
madImg=${VALpath}/ARCMFC_fig_val_ts_mad_lt012h_201901.png
biasImg=${VALpath}/ARCMFC_fig_val_ts_bias_lt012h_201901.png
corrImg=${VALpath}/ARCMFC_fig_val_ts_corr_lt012h_201901.png
SIImg=${VALpath}/ARCMFC_fig_val_ts_SI_lt012h_201901.png
novImg=${VALpath}/ARCMFC_fig_val_ts_nov_lt012h_201901.png
scatterImg=${VALpath}/ARCMFC_fig_val_scatter_lt012h_201901.png

# Create html document:
htmlFile=index.html
bname=test

cp head1.html ${htmlFile}
echo '<title>'${bName}'</title>'       >> ${htmlFile}
cat head2.html                         >> ${htmlFile}

#echo '<h2>Results for '${bName}'</h2>' >> ${htmlFile}

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
