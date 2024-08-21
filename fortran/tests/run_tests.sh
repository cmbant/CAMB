set -e

python --version
gfortran --version

pushd fortran

make clean
make

mkdir -p testfiles
python tests/CAMB_test_files.py testfiles --make_ini

pushd testfiles
echo "cloning test output"
git clone --depth=1 https://github.com/cmbant/CAMB_test_outputs.git
popd

python tests/CAMB_test_files.py testfiles --diff_to CAMB_test_outputs/test_outputs --verbose

rm -Rf testfiles/CAMB_test_outputs

popd
