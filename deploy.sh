#!/bin/bash

set -e # exit with nonzero exit code if anything fails

echo 'Travis Build Dir: '$TRAVIS_BUILD_DIR

if $TRAVIS_PULL_REQUEST
  then
    DEPLOY_BRANCH=$TRAVIS_PULL_REQUEST
else
  DEPLOY_BRANCH=$TRAVIS_BRANCH
fi

# Deploy built site to this branch
echo Deploy branch: $DEPLOY_BRANCH

# SSH URL of the RMG/RMG-tests repo that is pushed to:
REPO=git@github.com:ReactionMechanismGenerator/RMG-tests.git

if [ -n "$TRAVIS_BUILD_ID" ]; then

  #
  # Set the following environment variables in the travis configuration (.travis.yml)
  #
  #   DEPLOY_BRANCH    - The only branch that Travis should deploy from
  #   GIT_NAME         - The Git user name
  #   GIT_EMAIL        - The Git user email
  #
  GIT_NAME="Travis Deploy"
  GIT_EMAIL="travisci@rmg.edu"
  echo DEPLOY_BRANCH: $DEPLOY_BRANCH
  echo GIT_NAME: $GIT_NAME
  echo GIT_EMAIL: $GIT_EMAIL
  if [ "$TRAVIS_BRANCH" != "$DEPLOY_BRANCH" ]; then
    echo "Travis should only deploy from the DEPLOY_BRANCH ($DEPLOY_BRANCH) branch"
    exit 0
  else
    if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
      echo "Travis should not deploy from pull requests"
      exit 0
    else

      # use the decrypted deploy SSH key as
      # the credentials to push to the RMG-tests repo:
      chmod 600 deploy_key
      eval `ssh-agent -s`
      ssh-add deploy_key
      git config --global user.name "$GIT_NAME"
      git config --global user.email "$GIT_EMAIL"
    fi
  fi
fi

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
git commit --allow-empty -m $REV

# push to the branch to the RMG/RMG-tests repo:
git push -f $REPO $RMGTESTSBRANCH > /dev/null