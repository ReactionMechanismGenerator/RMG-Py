#!/bin/bash

set -e # exit with nonzero exit code if anything fails

branch=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
git push origin --delete $branch
