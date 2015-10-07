git config user.name "Travis CI"
git config user.email "antony@cosmologist.info"
git clone -b ${TRAVIS_BRANCH}_travis --depth 1 "https://$GH_TOKEN@github.com/cmbant/CAMB_test_outputs" > /dev/null 2>&1
cd CAMB_test_outputs
git rm -r *
cp -r ../testfiles/* .
git add --all .
git commit -m "Travis build $TRAVIS_BUILD_NUMBER for ${TRAVIS_BRANCH} commit $TRAVIS_COMMIT"
git push --quiet origin ${TRAVIS_BRANCH}_travis > /dev/null 2>&1

#git init
#git add --all .
#git commit -m "Travis build $TRAVIS_BUILD_NUMBER for ${TRAVIS_BRANCH} commit $TRAVIS_COMMIT"
#git push --force --quiet "https://$GH_TOKEN@github.com/cmbant/CAMB_test_outputs" master:${TRAVIS_BRANCH}_travis > /dev/null 2>&1

