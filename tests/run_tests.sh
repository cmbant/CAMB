
set -e

pushd pycamb

#copy clean files to fortran folder
mkdir -p fortran
find ../ -maxdepth 1 -type f | xargs cp -t fortran

source activate py3-environment
python --version
python setup.py install
python -c "import camb; print(camb.__version__)"
python setup.py test
pip uninstall -y camb
rm -Rf dist/*
rm -Rf build/*

source activate py2-environment
python --version
python setup.py install
python -c "import camb; print(camb.__version__)"
python setup.py test
pip uninstall -y camb
rm -Rf dist/*
rm -Rf build/*

popd

case "$TRAVIS_BRANCH" in
 devel*) 
       BRANCH="devel" 
       ;;
    *) 
       make camb EQUATIONS=equations_ppf
       BRANCH="master" 
       export CAMB_TESTS_NO_SOURCES=1
       export CAMB_TESTS_NO_DE=1
       ;;
esac


make clean
make Release

mkdir testfiles
python python/CAMB_test_files.py testfiles --make_ini

pushd testfiles
git clone -b $BRANCH --depth=1 https://github.com/cmbant/CAMB_test_outputs.git
popd

python python/CAMB_test_files.py testfiles --diff_to CAMB_test_outputs/test_outputs --verbose

rm -Rf testfiles/CAMB_test_outputs
