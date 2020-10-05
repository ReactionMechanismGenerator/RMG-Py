#!/bin/bash

set -e # exit with nonzero exit code if anything fails

if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
  echo "RMG-tests should not deploy from pull requests"
  exit 0
elif [ "$TRAVIS_BRANCH" == "master" ]; then
  echo "RMG-tests should not deploy from the master branch"
  exit 0
elif [ "$TRAVIS_BRANCH" == "stable" ]; then
  echo "RMG-tests should not deploy from the stable branch"
  exit 0
else
  DEPLOY_BRANCH=$TRAVIS_BRANCH
fi

echo "TRAVIS_BUILD_DIR: $TRAVIS_BUILD_DIR"
echo "DEPLOY_BRANCH: $DEPLOY_BRANCH"

# URL for the official RMG-tests repository
REPO=https://${GH_TOKEN}@github.com/ReactionMechanismGenerator/RMG-tests.git

# create a temporary folder:
REPO_NAME=$(basename $REPO)
TARGET_DIR=$(mktemp -d /tmp/$REPO_NAME.XXXX)
REV=$(git rev-parse HEAD)

# clone RMG-tests repo in the newly created folder:
git clone ${REPO} ${TARGET_DIR}

# go inside the newly created folder:
cd $TARGET_DIR

# create a new branch in RMG-tests with the name equal to
# the branch name of the tested RMG-Py branch:
RMGTESTSBRANCH=rmgpy-$DEPLOY_BRANCH

git checkout -b $RMGTESTSBRANCH || true 
git checkout $RMGTESTSBRANCH

# create an empty commit with the SHA-ID of the 
# tested commit of the RMG-Py branch:\
DB_DEPLOY_BRANCH="Metal_Attributes"
git commit --allow-empty -m rmgpydb-$REV-${DB_DEPLOY_BRANCH}

# push to the branch to the RMG/RMG-tests repo:
git push -f $REPO $RMGTESTSBRANCH > /dev/null
