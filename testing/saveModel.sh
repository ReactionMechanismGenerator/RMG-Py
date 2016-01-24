#!/bin/bash

echo 'Travis Build Dir: '$TRAVIS_BUILD_DIR

if $TRAVIS_PULL_REQUEST
  then
    GIT_BRANCH=$TRAVIS_PULL_REQUEST #"$(git rev-parse --abbrev-ref HEAD)"
else
  GIT_BRANCH=$TRAVIS_BRANCH
fi

echo 'Checked out RMG-Py branch: '$GIT_BRANCH

echo 'Creating RMG-tests branch...'
cd $TRAVIS_BUILD_DIR/../RMG-tests
git checkout -b $GIT_BRANCH || true 
git checkout $GIT_BRANCH

echo 'Committing to RMG-tests branch with RMG-Py commit:  '$TRAVIS_COMMIT
git commit --allow-empty -m $TRAVIS_COMMIT


