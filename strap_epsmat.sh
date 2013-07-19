grep "(" epsmat.out > epsinv
sed -i 's/(/ /g' epsinv
sed -i 's/)/ /g' epsinv
sed -i 's/,/ /g' epsinv