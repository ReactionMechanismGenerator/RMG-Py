#!/bin/bash

# This script saves the CHEMKIN file and dictionary associated with an RMG run
# It commits it to a folder in RMG-tests automatically

# Types of models to save
eg1=minimal
# eg2, etc...

#git_hash="$(git rev-parse --short=0 HEAD)"

if $TRAVIS_PULL_REQUEST
  then
    GIT_BRANCH=$TRAVIS_PULL_REQUEST #"$(git rev-parse --abbrev-ref HEAD)"
else
  GIT_BRANCH=$TRAVIS_BRANCH
fi

echo $GIT_BRANCH


SOURCE_FOLDER=$TRAVIS_BUILD_DIR/RMG-Py/testing/$eg1/chemkin/
DESTINATION_FOLDER=$TRAVIS_BUILD_DIR/RMG-tests/rmg/$eg1/

mkdir -p $DESTINATION_FOLDER || true
cd $DESTINATION_FOLDER
git checkout -b $GIT_BRANCH || true 
git checkout $GIT_BRANCH


# Copy CHEMKIN and species dictionary files
cp $SOURCE_FOLDER/chem_annotated.inp $DESTINATION_FOLDER
cp $SOURCE_FOLDER/species_dictionary.txt $DESTINATION_FOLDER

git commit -m $TRAVIS_COMMIT


