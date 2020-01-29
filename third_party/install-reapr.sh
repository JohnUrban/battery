tar -xzf Reapr_1.0.18-for-Battery.tar.gz

## Would need greadlink for this on Mac OS
BAMTOOLS_PATH=`readlink -f ${1}` ## e.g. ./bamtools

cd Reapr_1.0.18


bash install-battery-reapr.sh ${BAMTOOLS_PATH}
