"""
Optimized YAML input reader for RMG
Reads YAML format RMG input files and calls existing input.py functions
Preserves compatibility with legacy Python input files
"""

import yaml
import logging
import os
from pathlib import Path
from functools import lru_cache
from typing import Dict, Any, List, Union, Optional, Tuple

# import all existing functions from input.py
from rmgpy.rmg.input import (
    database, catalyst_properties, species, forbidden,
    simple_reactor, constant_V_ideal_gas_reactor, constant_TP_ideal_gas_reactor,
    liquid_cat_reactor, constant_T_V_liquid_reactor, liquid_reactor,
    surface_reactor, mb_sampled_reactor, simulator, solvation,
    model, quantum_mechanics, ml_estimator, pressure_dependence,
    options, generated_species_constraints, thermo_central_database,
    uncertainty, restart_from_seed, liquid_volumetric_mass_transfer_coefficient_power_law,
    smiles, inchi, adjacency_list, adjacency_list_group, smarts,
    fragment_adj, fragment_smiles, react
)

class YAMLInputReader:
    """
    Optimized YAML input reader with improved structure and performance
    (using mapping from dicts now)
    """
    
    # mapping of YAML keys to internal processor methods 
    # avoids long if/else chains and allows for scalable additions into input file and
    # corresponding processor methods here
    PROCESSORS = {
        'database': '_process_database',                    
        'catalystProperties': '_process_catalyst_properties', 
        'species': '_process_species',                      
        'forbidden': '_process_forbidden',                 
        'react': '_process_react',                          
        'simpleReactor': '_process_simple_reactor',         
        'constantVIdealGasReactor': '_process_constant_v_reactor',     
        'constantTPIdealGasReactor': '_process_constant_tp_reactor',   
        'liquidCatReactor': '_process_liquid_cat_reactor',  
        'constantTVLiquidReactor': '_process_constant_tv_liquid_reactor', 
        'liquidReactor': '_process_liquid_reactor',         
        'surfaceReactor': '_process_surface_reactor',       
        'mbSampledReactor': '_process_mb_sampled_reactor',  
        'solvation': '_process_solvation',                 
        'liquidVolumetricMassTransferCoefficientPowerLaw': '_process_volumetric_mass_transfer', 
        'simulator': '_process_simulator',                 
        'model': '_process_model',                          
        'quantumMechanics': '_process_quantum_mechanics',   
        'mlEstimator': '_process_ml_estimator',             
        'pressureDependence': '_process_pressure_dependence', 
        'generatedSpeciesConstraints': '_process_species_constraints', 
        'thermoCentralDatabase': '_process_thermo_central_database',
        'uncertainty': '_process_uncertainty',              
        'restartFromSeed': '_process_restart_from_seed',
        'options': '_process_options'                       
    }
    
    # mapping structure types to their corresponding conversion functions
    # already existing in input.py (the ones we imported up above)
    # allows flexible specification of molecular structures in different formats
    STRUCTURE_TYPES = {
        'SMILES': smiles,                           # SMILES string notation
        'InChI': inchi,                             # InChI string notation
        'adjacencyList': adjacency_list,            # RMG adjacency list format
        'adjacencyListGroup': adjacency_list_group, # RMG adjacency list for reaction groups
        'SMARTS': smarts,                           # SMARTS pattern notation
        'fragmentAdjacencyList': fragment_adj,      # Fragment adjacency list
        'fragmentSMILES': fragment_smiles           # Fragment SMILES notation
    }
    
    def __init__(self, path: Union[str, Path]):
        """
        Initialize YAML input reader with a file path

        Parameters
        ----------
        path : Union[str, Path]
            Path to the YAML input file
        """
        self.path = Path(path) # ensures path is a pathlib.Path obj for consistent path handling
        self.data = None # will store the parsed YAML data from input file
        self.species_dict = {} # store the species dictionary for reference lookup
        
    def read(self) -> None:
        """
        Read and parse YAML input file

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If YAML file has invalid syntax or structure
        IOError
            If file cannot be read
        """
        logging.info(f'Reading YAML input file "{self.path}"...') # log file being processed
        
        try:
            with open(self.path, 'r') as file:  
                content = file.read()            
                
            self.data = yaml.safe_load(content) # parse YAML content safely (prevents code execution)
            
            if not isinstance(self.data, dict): # validate that top-level structure is a dictionary
                raise ValueError("YAML file must contain a dictionary at the top level")
            
            logging.info(content) # log raw yaml text for debugging
            
        except yaml.YAMLError as e:
            # check for any YAML syntax errors
            raise ValueError(f"Invalid YAML syntax: {e}") 
        except IOError as e:
            # check to see if can access file or if file exists
            raise IOError(f"Failed to read file: {e}")  
    
    def process(self) -> None:
        """
        Process the loaded YAML data and call appropriate RMG input functions based
        on content of YAML input file

        Returns
        -------
        None

        Raises
        ------
        RuntimeError
            If no data has been loaded
        ValueError
            If error occurs while processing any section
        """
        if not self.data: # ensure data loaded before processing
            raise RuntimeError("No data loaded. Run read() first")
        
        # loop through all possible keys in order
        # replaced the previous method of using a bazillion if statements
        for key, processor_name in self.PROCESSORS.items(): # iterate through all processor mappings
            if key in self.data: # if the section exists in input file
                processor = getattr(self, processor_name) # call on corresponding processor func (dynamic method lookup)
                try:
                    processor(self.data[key]) # call on the processor function with the section data
                except Exception as e:
                    raise ValueError(f"Error processing {key}: {e}") # report processing errors with context
    
    def _process_structure(self, structure_data: Union[str, Dict], 
                          entity_type: str = "entity") -> Any:
        """
        Generic structure processor to reduce code duplication

        Parameters
        ----------
        structure_data : Union[str, Dict]
            Structure data (string for direct adjacency list or dict with type/value pairs)
        entity_type : str
            Type of entity for error messages ("species" or "forbidden")

        Returns
        -------
        Any
            Processed structure object (e.g., Molecule, Group)

        Raises
        ------
        ValueError
            If structure format is unknown or invalid
        """
        # default to adjacency list for species, adjacency list group for forbidden
        if isinstance(structure_data, str): # handle simple string format (most common case)
            # if just string, default to adjacency list by calling the imported conversion funcs
            if entity_type == "forbidden":
                return adjacency_list_group(structure_data) # forbidden structures use group format
            else:
                return adjacency_list(structure_data) # reg species use molecule format
                
        elif isinstance(structure_data, dict): # handle dictionary format
            # if is a dict, find the structure type and process it
            for key, func in self.STRUCTURE_TYPES.items(): # check each structure type
                if key in structure_data:
                    return func(structure_data[key]) # call on processor func for the found struc type
                    
            raise ValueError(f"Unknown structure format in {entity_type}") # no recognized format found
        else:
            raise ValueError(f"Invalid structure format for {entity_type}") # not string or dict
    
    ##############################################################
    # the general technique for these functions is just taking
    # the data from the YAML file and plugging them into the 
    # parameters of the imported pre-existing functions from the
    # original input handler
    ###############################################################

    def _process_database(self, db_data: Dict[str, Any]) -> None:
        """
        Process database input configuration

        Parameters
        ----------
        db_data : Dict[str, Any]
            Database configuration dictionary

        Returns
        -------
        None
        """
        # convert reaction libraries
        # handle reaction libraries which can be simple strings or dicts with seed flags
        reaction_libraries = [
            (lib.get('name'), lib.get('seed', False)) if isinstance(lib, dict) else lib # convert dict to tuple
            for lib in db_data.get('reactionLibraries', []) # auto default to empty list if not specified
        ]
        
        # call original database function 
        database(
            thermoLibraries=db_data.get('thermoLibraries'),
            transportLibraries=db_data.get('transportLibraries'),
            reactionLibraries=reaction_libraries,
            frequenciesLibraries=db_data.get('frequenciesLibraries'),
            seedMechanisms=db_data.get('seedMechanisms'),
            kineticsFamilies=db_data.get('kineticsFamilies', 'default'),
            kineticsDepositories=db_data.get('kineticsDepositories', 'default'),
            kineticsEstimator=db_data.get('kineticsEstimator', 'rate rules'),
            adsorptionGroups=db_data.get('adsorptionGroups', 'adsorptionPt111')
        )
    
    def _process_catalyst_properties(self, cat_data: Dict[str, Any]) -> None:
        """
        Process catalyst properties configuration

        Parameters
        ----------
        cat_data : Dict[str, Any]
            Catalyst properties configuration dictionary

        Returns
        -------
        None
        """
        # call original catalyst_properties function
        catalyst_properties(
            bindingEnergies=cat_data.get('bindingEnergies'),
            surfaceSiteDensity=cat_data.get('surfaceSiteDensity'),
            metal=cat_data.get('metal'),
            coverageDependence=cat_data.get('coverageDependence', False)
        )
    
    def _process_species(self, spec_list: List[Dict[str, Any]]) -> None:
        """
        Process species definitions

        Parameters
        ----------
        spec_list : List[Dict[str, Any]]
            List of species configuration dictionaries

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If species structure is missing or invalid
        """
        for spec in spec_list: # process each species in input file
            if 'structure' not in spec: # check if structure provided
                raise ValueError(f"No structure provided for species {spec.get('label', 'unknown')}")
            
            structure = self._process_structure(spec['structure'], "species") # convert structure to RMG format
            
            # call original species function
            species(
                label=spec['label'],
                structure=structure,
                reactive=spec.get('reactive', True),
                cut=spec.get('cut', False),
                size_threshold=spec.get('sizeThreshold')
            )
    
    def _process_forbidden(self, forb_list: List[Dict[str, Any]]) -> None:
        """
        Process forbidden structures

        Parameters
        ----------
        forb_list : List[Dict[str, Any]]
            List of forbidden structure configuration dictionaries

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If forbidden structure is missing or invalid
        """
        for forb in forb_list: # process each forb structure
            if 'structure' not in forb: # check if structure provided
                raise ValueError(f"No structure provided for forbidden {forb.get('label', 'unknown')}")
            
            structure = self._process_structure(forb['structure'], "forbidden") # convert to group format
            
            # call original forb function
            forbidden(
                label=forb['label'],
                structure=structure
            )
    
    def _process_react(self, react_data: Any) -> None:
        """
        Process react specifications

        Parameters
        ----------
        react_data : Any
            React configuration data

        Returns
        -------
        None
        """
        # pass to original react function
        react(react_data)
    
    def _process_reactor_common(self, reactor_func, reactor_data: Union[Dict, List],
                               field_mapping: Dict[str, str]) -> None:
        """
        Common reactor processor to reduce code duplication

        Parameters
        ----------
        reactor_func : callable
            The reactor function to call
        reactor_data : Union[Dict, List]
            Reactor configuration data (single dict or list of dicts)
        field_mapping : Dict[str, str]
            Mapping of YAML field names to function parameter names

        Returns
        -------
        None
        """
        reactors = reactor_data if isinstance(reactor_data, list) else [reactor_data] # normalize to list format
        
        for reactor in reactors: # process reactor(s) in input file
            kwargs = {} # keyword args dictionary for function call (basically parameter names of func)
            
            for yaml_key, param_name in field_mapping.items(): # map YAML keys to provided function parameters
                if yaml_key in reactor: # only process keys in input file
                    value = reactor[yaml_key]
                    
                    # special handling for sensitivity
                    if yaml_key == 'sensitivity' and value is not None: # sensitivity can be string or list
                        value = [value] if isinstance(value, str) else value # normalize to list format
                    
                    # convert quantities
                    elif any(keyword in yaml_key.lower() for keyword in # check if field requires unit conversion
                            ['temperature', 'pressure', 'time', 'volume', 'rate', 
                             'coefficient', 'viscosity', 'potential', 'distance']):
                        value = self._convert_quantity(value) # convert to (value, units) tuple format
                    
                    # convert concentration dictionaries
                    elif 'concentration' in yaml_key.lower() and isinstance(value, dict):
                        value = self._convert_concentration_dict(value)  # process concentration mappings
                    
                    kwargs[param_name] = value # store converted value with function parameter name
            
            reactor_func(**kwargs) # unpack the kwargs into params so func can handle
            # call on the func itself with the kwargs as params
    
    def _process_simple_reactor(self, reactor_data: Union[Dict, List]) -> None:
        """
        Process simple reactor configuration

        Parameters
        ----------
        reactor_data : Union[Dict, List]
            Simple reactor configuration data

        Returns
        -------
        None
        """
        # define mapping from YAML keys to function param for simple reactor
        # this basically what those kwargs above are 
        field_mapping = {
            'temperature': 'temperature',
            'pressure': 'pressure',
            'initialMoleFractions': 'initialMoleFractions',
            'nSims': 'nSims',
            'terminationConversion': 'terminationConversion',
            'terminationTime': 'terminationTime',
            'terminationRateRatio': 'terminationRateRatio',
            'balanceSpecies': 'balanceSpecies',
            'sensitivity': 'sensitivity',
            'sensitivityThreshold': 'sensitivityThreshold',
            'sensitivityTemperature': 'sensitivityTemperature',
            'sensitivityPressure': 'sensitivityPressure',
            'sensitivityMoleFractions': 'sensitivityMoleFractions',
            'constantSpecies': 'constantSpecies'
        }
        # use common processing
        self._process_reactor_common(simple_reactor, reactor_data, field_mapping)
    
    def _process_constant_v_reactor(self, reactor_data: Union[Dict, List]) -> None:
        """
        Process constant V ideal gas reactor configuration

        Parameters
        ----------
        reactor_data : Union[Dict, List]
            Constant V reactor configuration data

        Returns
        -------
        None
        """
        # mapping for constant V reactor parameters
        field_mapping = {
            'temperature': 'temperature',
            'pressure': 'pressure',
            'initialMoleFractions': 'initialMoleFractions',
            'terminationConversion': 'terminationConversion',
            'terminationTime': 'terminationTime',
            'terminationRateRatio': 'terminationRateRatio',
            'balanceSpecies': 'balanceSpecies'
        }
        self._process_reactor_common(constant_V_ideal_gas_reactor, reactor_data, field_mapping)
    
    def _process_constant_tp_reactor(self, reactor_data: Union[Dict, List]) -> None:
        """
        Process constant T,P ideal gas reactor configuration

        Parameters
        ----------
        reactor_data : Union[Dict, List]
            Constant T,P reactor configuration data

        Returns
        -------
        None
        """
        # mapping for constant T and P reactor
        field_mapping = {
            'temperature': 'temperature',
            'pressure': 'pressure',
            'initialMoleFractions': 'initialMoleFractions',
            'terminationConversion': 'terminationConversion',
            'terminationTime': 'terminationTime',
            'terminationRateRatio': 'terminationRateRatio',
            'balanceSpecies': 'balanceSpecies'
        }
        self._process_reactor_common(constant_TP_ideal_gas_reactor, reactor_data, field_mapping)
    
    def _process_liquid_cat_reactor(self, reactor_data: Union[Dict, List]) -> None:
        """
        Process liquid catalyst reactor configuration

        Parameters
        ----------
        reactor_data : Union[Dict, List]
            Liquid catalyst reactor configuration data

        Returns
        -------
        None
        """
        # mapping for liquid-phase catalytic reactor parameters
        field_mapping = {
            'temperature': 'temperature',
            'initialConcentrations': 'initialConcentrations',
            'initialSurfaceCoverages': 'initialSurfaceCoverages',
            'surfaceVolumeRatio': 'surfaceVolumeRatio',
            'distance': 'distance',
            'viscosity': 'viscosity',
            'surfPotential': 'surfPotential',
            'liqPotential': 'liqPotential',
            'terminationConversion': 'terminationConversion',
            'terminationTime': 'terminationTime',
            'terminationRateRatio': 'terminationRateRatio',
            'constantSpecies': 'constantSpecies'
        }
        self._process_reactor_common(liquid_cat_reactor, reactor_data, field_mapping)
    
    def _process_constant_tv_liquid_reactor(self, reactor_data: Union[Dict, List]) -> None:
        """
        Process constant T,V liquid reactor configuration

        Parameters
        ----------
        reactor_data : Union[Dict, List]
            Constant T,V liquid reactor configuration data

        Returns
        -------
        None
        """
        # mapping for constant T and V liquid reactor
        field_mapping = {
            'temperature': 'temperature',
            'initialConcentrations': 'initialConcentrations',
            'liquidVolume': 'liquidVolume',
            'residenceTime': 'residenceTime',
            'inletVolumetricFlowRate': 'inletVolumetricFlowRate',
            'outletVolumetricFlowRate': 'outletVolumetricFlowRate',
            'inletConcentrations': 'inletConcentrations',
            'vaporPressure': 'vaporPressure',
            'vaporMoleFractions': 'vaporMoleFractions',
            'terminationConversion': 'terminationConversion',
            'terminationTime': 'terminationTime',
            'terminationRateRatio': 'terminationRateRatio',
            'constantSpecies': 'constantSpecies'
        }
        self._process_reactor_common(constant_T_V_liquid_reactor, reactor_data, field_mapping)
    
    def _process_liquid_reactor(self, reactor_data: Union[Dict, List]) -> None:
        """
        Process liquid reactor configuration

        Parameters
        ----------
        reactor_data : Union[Dict, List]
            Liquid reactor configuration data

        Returns
        -------
        None
        """
        # mapping for general liquid reactor parameters
        field_mapping = {
            'temperature': 'temperature',
            'initialConcentrations': 'initialConcentrations',
            'terminationConversion': 'terminationConversion',
            'nSims': 'nSims',
            'terminationTime': 'terminationTime',
            'terminationRateRatio': 'terminationRateRatio',
            'sensitivity': 'sensitivity',
            'sensitivityThreshold': 'sensitivityThreshold',
            'sensitivityTemperature': 'sensitivityTemperature',
            'sensitivityConcentrations': 'sensitivityConcentrations',
            'constantSpecies': 'constantSpecies'
        }
        self._process_reactor_common(liquid_reactor, reactor_data, field_mapping)
    
    def _process_surface_reactor(self, reactor_data: Union[Dict, List]) -> None:
        """
        Process surface reactor configuration

        Parameters
        ----------
        reactor_data : Union[Dict, List]
            Surface reactor configuration data

        Returns
        -------
        None
        """
        # mapping for surface catalysis reactor parameters
        field_mapping = {
            'temperature': 'temperature',
            'initialPressure': 'initialPressure',
            'initialGasMoleFractions': 'initialGasMoleFractions',
            'initialSurfaceCoverages': 'initialSurfaceCoverages',
            'surfaceVolumeRatio': 'surfaceVolumeRatio',
            'nSims': 'nSims',
            'terminationConversion': 'terminationConversion',
            'terminationTime': 'terminationTime',
            'terminationRateRatio': 'terminationRateRatio',
            'sensitivity': 'sensitivity',
            'sensitivityThreshold': 'sensitivityThreshold'
        }
        self._process_reactor_common(surface_reactor, reactor_data, field_mapping)
    
    def _process_mb_sampled_reactor(self, reactor_data: Union[Dict, List]) -> None:
        """
        Process MB sampled reactor configuration

        Parameters
        ----------
        reactor_data : Union[Dict, List]
            MB sampled reactor configuration data

        Returns
        -------
        None
        """
        # mapping for Maxwell-Boltzmann sampled reactor (for T fluctuations)
        field_mapping = {
            'temperature': 'temperature',
            'pressure': 'pressure',
            'initialMoleFractions': 'initialMoleFractions',
            'mbsamplingRate': 'mbsamplingRate',
            'terminationConversion': 'terminationConversion',
            'terminationTime': 'terminationTime',
            'sensitivity': 'sensitivity',
            'sensitivityThreshold': 'sensitivityThreshold',
            'constantSpecies': 'constantSpecies'
        }
        self._process_reactor_common(mb_sampled_reactor, reactor_data, field_mapping)
                
    def _process_solvation(self, solv_data: Dict[str, Any]) -> None:
        """
        Process solvation settings

        Parameters
        ----------
        solv_data : Dict[str, Any]
            Solvation configuration dictionary

        Returns
        -------
        None

        Raises
        ------
        ImportError
            If SolventData cannot be imported when needed
        """
        # handle SolventData if provided
        solvent_data = None # initialize as None
        if 'solventData' in solv_data: # check if custom solvent data is provided in input file
            try:
                from rmgpy.data.solvation import SolventData # import solvent data class
            except ImportError:
                raise ImportError(
                    "SolventData could not be imported. Make sure RMG's solvation module is installed."
                )
            
            sd = solv_data['solventData'] # extract solvent data dictionary from input file
            # if custom solvent then specify Abraham-Mintz values
            solvent_data = SolventData( # create SolventData object with all params for custom sovlents
                # Abraham values 
                s_g=sd.get('s_g'),
                b_g=sd.get('b_g'),
                e_g=sd.get('e_g'),
                l_g=sd.get('l_g'),
                a_g=sd.get('a_g'),
                c_g=sd.get('c_g'),
                
                # solvent descriptors for enthalpy effects
                s_h=sd.get('s_h'),
                b_h=sd.get('b_h'),
                e_h=sd.get('e_h'),
                l_h=sd.get('l_h'),
                a_h=sd.get('a_h'),
                c_h=sd.get('c_h'),
                 # viscosity correlation coefficients
                A=sd.get('A'),
                B=sd.get('B'),
                C=sd.get('C'),
                D=sd.get('D'),
                E=sd.get('E'),
                 # additional solvent properties
                alpha=sd.get('alpha'),
                beta=sd.get('beta'),
                # dielectric constant
                eps=sd.get('eps'),
                # name lol
                name=sd.get('name')
            )
            
        # call original solvation function
        solvation( 
            solvent=solv_data['solvent'],
            solventData=solvent_data
        )
     
    def _process_volumetric_mass_transfer(self, lmt_data: Dict[str, Any]) -> None:
        """
        Process liquid volumetric mass transfer coefficient settings

        Parameters
        ----------
        lmt_data : Dict[str, Any]
            Mass transfer coefficient configuration dictionary

        Returns
        -------
        None
        """
        # call original function
        liquid_volumetric_mass_transfer_coefficient_power_law(
            prefactor=self._convert_quantity(lmt_data.get('prefactor', (0, "1/s"))),
            diffusionCoefficientPower=lmt_data.get('diffusionCoefficientPower', 0),
            solventViscosityPower=lmt_data.get('solventViscosityPower', 0),
            solventDensityPower=lmt_data.get('solventDensityPower', 0)
        )
        
    def _process_simulator(self, sim_data: Dict[str, Any]) -> None:
        """
        Process simulator settings

        Parameters
        ----------
        sim_data : Dict[str, Any]
            Simulator configuration dictionary

        Returns
        -------
        None
        """
        # call original simulator function
        simulator(
            atol=float(sim_data.get('atol', 1e-16)),
            rtol=float(sim_data.get('rtol', 1e-8)),
            sens_atol=float(sim_data.get('sens_atol', 1e-6)),
            sens_rtol=float(sim_data.get('sens_rtol', 1e-4))
        )
    
    def _process_model(self, model_data: Dict[str, Any]) -> None:
        """
        Process model settings

        Parameters
        ----------
        model_data : Dict[str, Any]
            Model configuration dictionary

        Returns
        -------
        None
        """
        # call og model func
        model(
            toleranceMoveToCore=model_data.get('toleranceMoveToCore'),
            toleranceRadMoveToCore=model_data.get('toleranceRadMoveToCore', float('inf')),
            toleranceMoveEdgeReactionToCore=model_data.get('toleranceMoveEdgeReactionToCore', float('inf')),
            toleranceKeepInEdge=model_data.get('toleranceKeepInEdge', 0.0),
            toleranceInterruptSimulation=model_data.get('toleranceInterruptSimulation', 1.0),
            toleranceMoveEdgeReactionToSurface=model_data.get('toleranceMoveEdgeReactionToSurface', float('inf')),
            toleranceMoveSurfaceSpeciesToCore=model_data.get('toleranceMoveSurfaceSpeciesToCore', float('inf')),
            toleranceMoveSurfaceReactionToCore=model_data.get('toleranceMoveSurfaceReactionToCore', float('inf')),
            toleranceMoveEdgeReactionToSurfaceInterrupt=model_data.get('toleranceMoveEdgeReactionToSurfaceInterrupt'),
            toleranceMoveEdgeReactionToCoreInterrupt=model_data.get('toleranceMoveEdgeReactionToCoreInterrupt'),
            maximumEdgeSpecies=model_data.get('maximumEdgeSpecies', 1000000),
            minCoreSizeForPrune=model_data.get('minCoreSizeForPrune', 50),
            minSpeciesExistIterationsForPrune=model_data.get('minSpeciesExistIterationsForPrune', 2),
            filterReactions=model_data.get('filterReactions', False),
            filterThreshold=model_data.get('filterThreshold', 1e8),
            ignoreOverallFluxCriterion=model_data.get('ignoreOverallFluxCriterion', False),
            maxNumSpecies=model_data.get('maxNumSpecies'),
            maxNumObjsPerIter=model_data.get('maxNumObjsPerIter', 1),
            terminateAtMaxObjects=model_data.get('terminateAtMaxObjects', False),
            toleranceThermoKeepSpeciesInEdge=model_data.get('toleranceThermoKeepSpeciesInEdge', float('inf')),
            dynamicsTimeScale=self._convert_quantity(model_data.get('dynamicsTimeScale', (0.0, 'sec'))),
            toleranceBranchReactionToCore=model_data.get('toleranceBranchReactionToCore', 0.0),
            branchingIndex=model_data.get('branchingIndex', 0.5),
            branchingRatioMax=model_data.get('branchingRatioMax', 1.0),
            toleranceTransitoryDict=model_data.get('toleranceTransitoryDict', {}),
            transitoryStepPeriod=model_data.get('transitoryStepPeriod', 20),
            toleranceReactionToCoreDeadendRadical=model_data.get('toleranceReactionToCoreDeadendRadical', 0.0)
        )
        
    def _process_quantum_mechanics(self, qm_data: Dict[str, Any]) -> None:
        """
        Process quantum mechanics settings

        Parameters
        ----------
        qm_data : Dict[str, Any]
            Quantum mechanics configuration dictionary

        Returns
        -------
        None
        """
        # call og qm func
        quantum_mechanics(
            software=qm_data['software'],
            method=qm_data['method'],
            fileStore=qm_data.get('fileStore'),
            scratchDirectory=qm_data.get('scratchDirectory'),
            onlyCyclics=qm_data.get('onlyCyclics', False),
            maxRadicalNumber=qm_data.get('maxRadicalNumber', 0)
        )
        
    def _process_ml_estimator(self, ml_data: Dict[str, Any]) -> None:
        """
        Process ML estimator settings

        Parameters
        ----------
        ml_data : Dict[str, Any]
            ML estimator configuration dictionary

        Returns
        -------
        None
        """
        # call og ML estimator func
        ml_estimator(
            thermo=ml_data.get('thermo', True),
            name=ml_data.get('name', 'main'),
            minHeavyAtoms=ml_data.get('minHeavyAtoms', 1),
            maxHeavyAtoms=ml_data.get('maxHeavyAtoms'),
            minCarbonAtoms=ml_data.get('minCarbonAtoms', 0),
            maxCarbonAtoms=ml_data.get('maxCarbonAtoms'),
            minOxygenAtoms=ml_data.get('minOxygenAtoms', 0),
            maxOxygenAtoms=ml_data.get('maxOxygenAtoms'),
            minNitrogenAtoms=ml_data.get('minNitrogenAtoms', 0),
            maxNitrogenAtoms=ml_data.get('maxNitrogenAtoms'),
            onlyCyclics=ml_data.get('onlyCyclics', False),
            onlyHeterocyclics=ml_data.get('onlyHeterocyclics', False),
            minCycleOverlap=ml_data.get('minCycleOverlap', 0),
            H298UncertaintyCutoff=self._convert_quantity(ml_data.get('H298UncertaintyCutoff', (3.0, 'kcal/mol'))),
            S298UncertaintyCutoff=self._convert_quantity(ml_data.get('S298UncertaintyCutoff', (2.0, 'cal/(mol*K)'))),
            CpUncertaintyCutoff=self._convert_quantity(ml_data.get('CpUncertaintyCutoff', (2.0, 'cal/(mol*K)')))
        )
        
    def _process_pressure_dependence(self, pd_data: Dict[str, Any]) -> None:
        """
        Process pressure dependence settings

        Parameters
        ----------
        pd_data : Dict[str, Any]
            Pressure dependence configuration dictionary

        Returns
        -------
        None
        """
        # process T, can be dict with min/max/count/units or list
        temps = pd_data['temperatures'] # get T specification
        if isinstance(temps, dict): # handle range specification
            temperatures = [temps['min'], temps['max'], temps['units'], temps['count']] # convert to list
        else:
            temperatures = temps # use list as is
            
        # process P, can be dict with min/max/count/units or list lkie T
        press = pd_data['pressures'] # get P specification
        if isinstance(press, dict):  # handle range specification
            pressures = [press['min'], press['max'], press['units'], press['count']] # convert to list 
        else:
            pressures = press # use list as is
            
        # process interpolation, can be list or tuple
        interp = pd_data.get('interpolation') # het interpolation method
        if isinstance(interp, list) and len(interp) > 1:  # handle list format
            # convert list format [method, param1, param2] to tuple
            interpolation = tuple(interp) # convert to tuple for function
        else:
            interpolation = interp  # use as is
            
        # call og pressure dependence func
        pressure_dependence(
            method=pd_data['method'],
            temperatures=temperatures,
            pressures=pressures,
            maximumGrainSize=self._convert_quantity(pd_data.get('maximumGrainSize', 0.0)),
            minimumNumberOfGrains=pd_data.get('minimumNumberOfGrains', 0),
            interpolation=interpolation,
            maximumAtoms=pd_data.get('maximumAtoms')
        )
        
    def _process_species_constraints(self, constraints: Dict[str, Any]) -> None:
        """
        Process generated species constraints

        Parameters
        ----------
        constraints : Dict[str, Any]
            Species constraints configuration dictionary

        Returns
        -------
        None
        """
        # create copy to avoid modifying original
        constraints_copy = constraints.copy() # more safety to avoid modifying input data
        
        # handle the special 'allowed' field
        if 'allowed' in constraints_copy: # process allowed specs list
            allowed_list = constraints_copy['allowed']
            # convert special string values to expected format
            processed_allowed = []  
            for item in allowed_list: # check each allowed item
                if item == 'input species':
                    processed_allowed.append('input species') # keep as string
                elif item == 'seed mechanisms':
                    processed_allowed.append('seed mechanisms')   
                elif item == 'reaction libraries':
                    processed_allowed.append('reaction libraries')
                else:
                    processed_allowed.append(item) # keep others as is
            constraints_copy['allowed'] = processed_allowed  # update w processed list
            
        # call og func
        generated_species_constraints(**constraints_copy) # unpack constraints like w kwargs above
        
    def _process_thermo_central_database(self, tcd_data: Dict[str, Any]) -> None:
        """
        Process thermo central database settings

        Parameters
        ----------
        tcd_data : Dict[str, Any]
            Thermo central database configuration dictionary

        Returns
        -------
        None
        """
        # call og thermo database func
        thermo_central_database(
            host=tcd_data['host'],
            port=tcd_data['port'],
            username=tcd_data['username'],
            password=tcd_data['password'],
            application=tcd_data['application']
        )
        
    def _process_uncertainty(self, unc_data: Dict[str, Any]) -> None:
        """
        Process uncertainty analysis settings

        Parameters
        ----------
        unc_data : Dict[str, Any]
            Uncertainty analysis configuration dictionary

        Returns
        -------
        None
        """
        # call og uncertainty func
        uncertainty(
            localAnalysis=unc_data.get('localAnalysis', False),
            globalAnalysis=unc_data.get('globalAnalysis', False),
            uncorrelated=unc_data.get('uncorrelated', True),
            correlated=unc_data.get('correlated', True),
            localNumber=unc_data.get('localNumber', 10),
            globalNumber=unc_data.get('globalNumber', 5),
            terminationTime=self._convert_quantity(unc_data.get('terminationTime')),
            pceRunTime=unc_data.get('pceRunTime', 1800),
            pceErrorTol=unc_data.get('pceErrorTol'),
            pceMaxEvals=unc_data.get('pceMaxEvals'),
            logx=unc_data.get('logx', True)
        )
        
    def _process_restart_from_seed(self, restart_data: Dict[str, Any]) -> None:
        """
        Process restart from seed settings

        Parameters
        ----------
        restart_data : Dict[str, Any]
            Restart from seed configuration dictionary

        Returns
        -------
        None
        """
        # call og restart func
        restart_from_seed(
            path=restart_data.get('path'),
            coreSeed=restart_data.get('coreSeed'),
            edgeSeed=restart_data.get('edgeSeed'),
            filters=restart_data.get('filters'),
            speciesMap=restart_data.get('speciesMap')
        )
        
    def _process_options(self, opt_data: Dict[str, Any]) -> None:
        """
        Process general RMG options

        Parameters
        ----------
        opt_data : Dict[str, Any]
            General options configuration dictionary

        Returns
        -------
        None
        """
        # call og options func
        options(
            name=opt_data.get('name', 'Seed'),
            generateSeedEachIteration=opt_data.get('generateSeedEachIteration', True),
            saveSeedToDatabase=opt_data.get('saveSeedToDatabase', False),
            units=opt_data.get('units', 'si'),
            saveRestartPeriod=opt_data.get('saveRestartPeriod'),
            generateOutputHTML=opt_data.get('generateOutputHTML', False),
            generatePlots=opt_data.get('generatePlots', False),
            saveSimulationProfiles=opt_data.get('saveSimulationProfiles', False),
            verboseComments=opt_data.get('verboseComments', False),
            saveEdgeSpecies=opt_data.get('saveEdgeSpecies', False),
            keepIrreversible=opt_data.get('keepIrreversible', False),
            trimolecularProductReversible=opt_data.get('trimolecularProductReversible', True),
            wallTime=opt_data.get('wallTime', '00:00:00:00'),
            saveSeedModulus=opt_data.get('saveSeedModulus', -1)
        )
        
    def _convert_quantity(self, value: Union[Dict, List, Tuple, str, float, None]) -> Union[Tuple, Any]:
        """
        Convert YAML quantity representation to tuple format expected by RMG functions

        Parameters
        ----------
        value : Union[Dict, List, Tuple, str, float, None]
            Quantity value in various formats

        Returns
        -------
        Union[Tuple, Any]
            Converted quantity as (value, units) tuple or original value

        """
        if value is None:
            return None
            
        if isinstance(value, dict): # handle dict
            if 'value' in value and 'units' in value: # standard {value: X, units: Y} format
                return (value['value'], value['units']) # convert to tuple
            # also handle case where quantity IS the dict value
            elif len(value) == 1: # single key value pair format
                # e.g., {0.5: 'kcal/mol'} format
                val, unit = next(iter(value.items())) # extract key-value pair
                return (val, unit)  # return as tuple
            else:
                # return the dict as is if it dont match expected formats
                return value # no conversion possible
        elif isinstance(value, (list, tuple)): # handle list/tuple
            if len(value) == 2: # standard (value, units) format
                return tuple(value)  # make tuple format
            elif len(value) == 4: # for T/P ranges in P dependence
                return value # keep as is (for range specs)
            else:
                # other list formats
                return value  # no conversion
        elif isinstance(value, str): # handle string
            # for cases like wallTime which is just a string
            return value  # return string as is
        else:
            # for single numeric values, return as-is
            # func being called will handle unit defaults if needed
            return value  # no conversion
            
    def _convert_concentration_dict(self, conc_dict: Dict[str, Any]) -> Dict[str, Tuple]:
        """
        Convert concentration dictionary with quantity values

        Parameters
        ----------
        conc_dict : Dict[str, Any]
            Dictionary with species names as keys and quantities as values

        Returns
        -------
        Dict[str, Tuple]
            Dictionary with converted quantity tuples
        """
        if not conc_dict:
            return {}
            
        result = {}
        for species, conc in conc_dict.items(): # process each spec concentration
            result[species] = self._convert_quantity(conc) # conv concentration to tuple
        return result  # return processed dict
     

###################################
# the actual reader function itself
###################################

def read_yaml_input_file(path: Union[str, Path], rmg0) -> None:
    """
    Read an RMG YAML input file and process it using existing input.py functions

    Parameters
    ----------
    path : Union[str, Path]
        Path to the YAML input file
    rmg0 : RMG
        RMG object to populate with input data

    Returns
    -------
    None
    """
    # import necessary modules for processing
    from rmgpy.rmg.input import set_global_rmg
    from rmgpy.rmg.model import CoreEdgeReactionModel
    
    # set up global RMG object
    set_global_rmg(rmg0)
    rmg0.reaction_model = CoreEdgeReactionModel()
    rmg0.initial_species = []
    rmg0.reaction_systems = []
    
    # clear the global species_dict
    from rmgpy.rmg import input as rmg_input
    rmg_input.species_dict = {} # clear global spec dict
    rmg_input.mol_to_frag = {} # clear molecular fragment
    
    # set spec constraints default
    rmg0.species_constraints = {'explicitlyAllowedMolecules': []}  # initialize with empty allowed list
    
    # process YAML input file
    reader = YAMLInputReader(path)
    reader.read() # read and parse YAML input file
    reader.process() # process parsed data and call RMG functions
    
    # post-processing (similar to original read_input_file)
    for reaction_system in rmg0.reaction_systems: # process each reactor system
        if hasattr(reaction_system, 'convert_initial_keys_to_species_objects'): # check for conversion method
            reaction_system.convert_initial_keys_to_species_objects(rmg_input.species_dict) # convert spec keys to objects
    
    if rmg0.quantum_mechanics:  # if quantum mechanics is enabled
        rmg0.quantum_mechanics.set_default_output_directory(rmg0.output_directory) # set qm output directory
        rmg0.quantum_mechanics.initialize() # initialize qm calculations
    
    logging.info('') # log empty line for spacing


def read_input_file_wrapper(path: Union[str, Path], rmg0) -> None:
    """
    Read an RMG input file (either Python or YAML format) and process it

    Parameters
    ----------
    path : Union[str, Path]
        Path to the input file (.py or .yaml/.yml)
    rmg0 : RMG
        RMG object to populate with input data

    Returns
    -------
    None

    Raises
    ------
    IOError
        If the input file cannot be found
    ValueError
        If the file format is unsupported
    """
    import os
    from pathlib import Path
    
    # get the file extension
    file_path = Path(path)
    extension = file_path.suffix.lower() # extract extension
    
    # check if file exists
    if not file_path.exists(): # validate file existence
        raise IOError(f'The input file "{path}" could not be found.')
    
    # route to appropriate reader based on extension
    # if .py (og input file) then run og input reader func
    if extension == '.py':
        from rmgpy.rmg.input import read_input_file as read_python_input_file
        logging.info(f'Detected Python input file format (.py)') # log file type
        read_python_input_file(path, rmg0)
    # if .yaml or .yml, use the new func above to read
    elif extension in ['.yaml', '.yml']:
        logging.info(f'Detected YAML input file format ({extension})') # log file type
        read_yaml_input_file(path, rmg0)
    else:
        raise ValueError(
            f'Unsupported input file format "{extension}". '
            f'RMG supports .py and .yaml/.yml input files.'
        )