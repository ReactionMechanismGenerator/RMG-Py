.. _isotopes:

********
Isotopes
********

Describing isotopes in adjacency lists
--------------------------------------

Isotopic enrichment can be indicated in a molecular structure's adjacency list. 
The example below is methane with an isotopically labeled carbon of isotope 
number 13, which is indicated with ``i13``::

    1 C u0 p0 c0 i13 {2,S} {3,S} {4,S} {5,S}
    2 H u0 p0 c0 {1,S}
    3 H u0 p0 c0 {1,S}
    4 H u0 p0 c0 {1,S}
    5 H u0 p0 c0 {1,S}


Running the RMG isotopes algorithm
----------------------------------

The isotopes script is located in the folder scripts. To run the algorithm, 
ensure the RMG packages are loaded and type::

    python /path/to/rmg/scripts/isotopes.py /path/to/input/file.py

The input file is identical to a standard RMG input file and should contain the
conditions you want to run (unless you are inputting an already completed RMG
model). Without any options, the script will run the original RMG input file to
generate a model. Once the RMG job is finished, it will create new species for
all isotopologues of previously generated species and then generate all
reactions between the isotopologues.

Some arguments can be used to alter the behavior of the script. If you already
have a model (which includes atom mapping in RMG's format) which you would like
to add isotope labels to, you can use the option ``--original path/to/model/directory``
with the desired model files stored within with structure ``chemkin/chem_annotated.inp``
and ``chemkin/species_dictionary.txt``. With this option, the isotope script
will use the specified model instead of re-running an RMG job.

If you only desire the reactions contained in the specific RMG job,
you can add ``--useOriginalReactions`` in addition to ``--original``.
This will create a full set of isotopically labeled versions of the reactions
you input and avoid a time-consuming generate reactions proceedure.

The arguement ``--maximumIsotopicAtoms [integer]`` limits the number of enriched
atoms in any isotopologue in the model. This is beneficial for decreasing model 
size, runtime of model creation and runtime necessary for analysis.

Adding kinetic isotope effects which are described in this paper can be obtained
through the argument ``--kineticIsotopeEffect simple``. Currently this is the
only supported method, though this can be extended to other effects.

If you have a desired output folder, ``--output output_folder_name`` can direct
all output files to the specified folder. The default is to create a folder
named 'iso'.

There are some limitations in what can be used in isotope models. In general,
RMG Reaction libraries and other methods of kinetic estimation that do not
involve atom mapping to reaction recipes are not compatible (though they can be
functional if all isotopologues are included in the reaction library). The
algorithm also does not function with pressure dependent mechanisms generated
by RMG, and has only been tested for gas phase kinetics. This algorithm currently
only works for Carbon-13 enrichments.

Following the generation, a number of diagnostics check model accuracy.
Isotopologues are checked to ensure their symmetries are consistent.
Then, the reaction path degeneracy among reactions differing only in isotope
labeling is checked to ensure it is consistent with the symmetry values of reactions.
If one of these checks throws a warning, the model will likely exhibit non-natural
fluctuations in enrichment ten to one hundred times larger than from non-hydrogen
kinetic isotope effects.

Output from script
------------------

The isotope generation script will output two files inside the nested folders
``iso/chemkin``, unless ``--output`` is specified. The file
``species_dictionary.txt`` lists the structure of all isotopologue using the
RMG adjacency list structure. The other file of importance ``chem_annotated.inp``
is a chemkin input file containing elements, species, thermo, and reactions of
the entire system.
