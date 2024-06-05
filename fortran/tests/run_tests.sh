set -e

gfortran --version
python --version
pip install -e .
python -c "import camb; print(camb.__version__)"
python -m unittest camb.tests.camb_test
rm -Rf HMcode_test_outputs
git clone https://github.com/alexander-mead/HMcode_test_outputs.git
python -m unittest camb.tests.hmcode_test
rm -Rf HMcode_test_outputs
pip uninstall -y camb
rm -Rf dist/*
rm -Rf build/*
rm -f camb/*.so

case "$TRAVIS_BRANCH" in
 devel*)
       BRANCH="master"
       ;;
    *)
       BRANCH="master"
       ;;
esac

pushd fortran

make clean
make

mkdir -p testfiles
python tests/CAMB_test_files.py testfiles --make_ini

pushd testfiles
echo "cloning test output branch:" $BRANCH
git clone -b $BRANCH --depth=1 https://github.com/cmbant/CAMB_test_outputs.git
popd

python tests/CAMB_test_files.py testfiles --diff_to CAMB_test_outputs/test_outputs --verbose

rm -Rf testfiles/CAMB_test_outputs

popd
