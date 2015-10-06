cd testfiles
rm -Rf CAMB_test_outputs
git init
git config user.name "Travis CI"
git config user.email "antony@cosmologist.info"
git add --all
git commit -m "Travis output for build $TRAVIS_BUILD_NUMBER"
git push --force --quiet "https://$GH_TOKEN@github.com/cmbant/CAMB_test_outputs" master:$(TRAVIS_BRANCH)_travis > /dev/null 2>&1

