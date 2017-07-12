cd /camb

make camb EQUATIONS=equations_ppf || exit $?
make clean || exit $?
make Release || exit $?
make libcamb || exit $?

cd pycamb

python setup.py install
python setup.py test || exit $?
cd ..

python python/CAMB_test_files.py testfiles --make_ini
cd testfiles
case "$TRAVIS_BRANCH" in
 devel) BRANCH="devel" ;;
    *) BRANCH="master" ;;
esac
git clone -b $BRANCH --depth=1 https://github.com/cmbant/CAMB_test_outputs.git
cd ..

python python/CAMB_test_files.py testfiles --diff_to CAMB_test_outputs/test_outputs --verbose
rc=$?

rm -Rf testfiles/CAMB_test_outputs
exit $rc
