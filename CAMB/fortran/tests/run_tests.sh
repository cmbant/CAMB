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

if [[ $TRAVIS_REPO_SLUG == "cmbant/CAMB" && $PYPI_DIST == "true" && "$TRAVIS_PULL_REQUEST" == "false" ]]
then
 case "$TRAVIS_BRANCH" in
 devel*) export CAMB_PACKAGE_NAME=camb_devel ;;
    *) export CAMB_PACKAGE_NAME=camb
 esac
 python setup.py sdist
 pip install twine
 twine upload -r pypitest --repository-url https://test.pypi.org/legacy/ dist/* || true
#too much delay on test.pypi to reliably immediately test install
# mkdir -p test_dir
# pushd test_dir
# pip install --index-url https://test.pypi.org/simple/ $CAMB_PACKAGE_NAME
# python -c "import camb; print(camb.__version__)"
# python -m unittest tests.camb_test
# pip uninstall -y $CAMB_PACKAGE_NAME
# popd
fi

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
