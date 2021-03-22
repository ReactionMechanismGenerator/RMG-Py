#!/bin/bash

set -e # exit with nonzero exit code if anything fails

git config --global user.name "RMG Bot"
git config --global user.email "rmg_dev@mit.edu"

BRANCH=${GITHUB_REF#refs/heads/}

if [ "$GITHUB_EVENT_NAME" != "push" ]; then
  echo "RMG-tests should not deploy from pull requests"
  exit 0
elif [ "$BRANCH" == "master" ]; then
  echo "RMG-tests should not deploy from the master branch"
  exit 0
elif [ "$BRANCH" == "stable" ]; then
  echo "RMG-tests should not deploy from the stable branch"
  exit 0
else
  DEPLOY_BRANCH=$BRANCH
fi

echo "GITHUB_WORKSPACE: $GITHUB_WORKSPACE"
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
# tested commit of the RMG-Py branch:
git commit --allow-empty -m rmgpy-$REV

# push to the branch to the RMG/RMG-tests repo:
git push -f $REPO $RMGTESTSBRANCH > /dev/null
