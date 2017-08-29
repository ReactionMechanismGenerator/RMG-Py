#!/bin/bash

target=$1

if [ -z ${TRAVIS_BUILD_DIR+x} ]; then TRAVIS_BUILD_DIR=$PWD; fi

echo 'Travis Build Dir: '$TRAVIS_BUILD_DIR

benchmarkmodel=$TRAVIS_BUILD_DIR/testing/benchmark/$target
testmodel=$TRAVIS_BUILD_DIR/testing/testmodel/$target
echo 'benchmark model folder: '$benchmarkmodel
echo 'Test model folder: '$testmodel

# check generated models:
# core:
python $RMG/scripts/checkModels.py $target $benchmarkmodel/chemkin/chem_annotated.inp $benchmarkmodel/chemkin/species_dictionary.txt $testmodel/chemkin/chem_annotated.inp $testmodel/chemkin/species_dictionary.txt

echo core for $target:
if grep "checkModels" $target.log | cut -f2- -d'=' > $target.core ; then
	cat $target.core
	rm $target.core
fi

# edge:
python $RMG/scripts/checkModels.py $target $benchmarkmodel/chemkin/chem_edge_annotated.inp $benchmarkmodel/chemkin/species_edge_dictionary.txt $testmodel/chemkin/chem_edge_annotated.inp $testmodel/chemkin/species_edge_dictionary.txt
echo edge for $target:
if grep "checkModels" $target.log | cut -f2- -d'=' > $target.edge ; then
	cat $target.edge
	rm $target.edge
fi

echo 'Execution time, Benchmark:'
grep "Execution time" $benchmarkmodel/RMG.log | tail -1
echo 'Execution time, Tested:'
grep "Execution time" $testmodel/RMG.log | tail -1

echo 'Memory used, Benchmark:'
grep "Memory used:" $benchmarkmodel/RMG.log | tail -1
echo 'Memory used, Tested:'
grep "Memory used:" $testmodel/RMG.log | tail -1

# regression testing
regr=examples/rmg/$target/regression_input.py
if [ -f "$regr" ];
then
	python $RMG/rmgpy/tools/regression.py $regr $benchmarkmodel/chemkin $testmodel/chemkin/
else
	echo "Regression input file not found. Not running a regression test."
fi