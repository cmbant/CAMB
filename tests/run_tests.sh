cd /camb

gfortran -dumpversion

make camb EQUATIONS=equations_ppf || exit $?
make clean || exit $?
make Release || exit $?
make libcamb || exit $?

pushd pycamb
python setup.py install || exit $?
python setup.py test || exit $?

conda create -q -y -n py3-environment python=3 numpy scipy sympy six
source activate py3-environment 
mkdir -p fortran
find ../ -maxdepth 1 -type f | xargs cp -t fortran
mv ../lensing.f90 ../lensing.f90_tmp
python setup.py install || exit $?
python -c "import camb; print(camb.__version__)"
python setup.py test || exit $?
mv ../lensing.f90_tmp ../lensing.f90
popd

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
