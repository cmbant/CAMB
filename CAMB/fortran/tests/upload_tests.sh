pushd fortran

git clone -b ${TRAVIS_BRANCH}_travis --depth 1 "https://$GH_TOKEN@github.com/cmbant/CAMB_test_outputs" > /dev/null 2>&1

if [ $? -ne 0 ]; then
mkdir CAMB_test_outputs
cd CAMB_test_outputs
git init
git checkout -b ${TRAVIS_BRANCH}_travis
git remote add origin https://$GH_TOKEN@github.com/cmbant/CAMB_test_outputs
else
cd CAMB_test_outputs
fi

git config user.name "Travis CI"
git config user.email "none@travis"
git rm -r * > /dev/null 2>&1
cp -r ../testfiles/* .
git add --all .
git commit -m "Travis build $TRAVIS_BUILD_NUMBER for $TRAVIS_BRANCH commit $TRAVIS_COMMIT"
git push --quiet origin ${TRAVIS_BRANCH}_travis
#> /dev/null 2>&1

popd