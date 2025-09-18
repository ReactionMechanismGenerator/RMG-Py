"""
YAML input reader for RMG
this module reads YAML format RMG input files
and calls the existing input.py functions

preserves ability to still use legacy python input files
"""

import yaml
import logging
import os
from pathlib import Path

# import ALL the existing functions from original input.py
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
    def __init__(self, path):
        """
        initialize the YAML input reader with a file path
        
        :param path: path to the YAML input file
        """
        self.path = Path(path)
        self.data = None
        self.species_dict = {}

    def read(self):
        """
        read and parse YAML input file
        """
        logging.info(f'Reading YAML input file  "{self.path}"...')
        # let user know in termal that file is being read

        with open(self.path, 'r') as file:
            #store content to log onto terminal
            content = file.read()

            self.data = yaml.safe_load(file)

        # check if input file is empty
        if not isinstance(self.data, dict):
            raise ValueError("Yo ur missing a dictionary bro YAML file needs a dictionary")
        
        # log contents of file into terminal
        logging.info(content)

    def process(self):
        """
        processes the info from YAML file and calls on the preexisting functions in input.py
        """
        if not self.data:
            raise RuntimeError("Yo no data loaded run read() first")
        
        # now process the data in sections (same order as for .py inputs)
        # lowkey a really monkey way of doing this w a ton of if statements but I cant really do 
        # switch cases since i need to check if EVERY PIECE OF INFO IS VALID 
        # ⬇ libraries and species information ⬇
        if 'database' in self.data:
            self._process_database(self.data['database'])
        if 'catalystProperties' in self.data:
            self._process_catalyst_properties(self.data['catalystProperties'])
        if 'species' in self.data:
            self._process_species(self.data['species'])
        if 'forbidden' in self.data:
            self._process_forbidden(self.data['forbidden'])
        if 'react' in self.data:
            self._process_react(self.data['react'])
        # process each reactor type individually 
        # each reactor can appear at the top level
        # ⬇ reactor information ⬇
        if 'simpleReactor' in self.data:
            self._process_simple_reactor(self.data['simpleReactor'])
        if 'constantVIdealGasReactor' in self.data:
            self._process_constant_v_reactor(self.data['constantVIdealGasReactor'])
        if 'constantTPIdealGasReactor' in self.data:
            self._process_constant_tp_reactor(self.data['constantTPIdealGasReactor'])
        if 'liquidCatReactor' in self.data:
            self._process_liquid_cat_reactor(self.data['liquidCatReactor'])
        if 'constantTVLiquidReactor' in self.data:
            self._process_constant_tv_liquid_reactor(self.data['constantTVLiquidReactor'])
        if 'liquidReactor' in self.data:
            self._process_liquid_reactor(self.data['liquidReactor'])
        if 'surfaceReactor' in self.data:
            self._process_surface_reactor(self.data['surfaceReactor'])
        if 'mbSampledReactor' in self.data:
            self._process_mb_sampled_reactor(self.data['mbSampledReactor'])
        # ⬇ other info/options ⬇ 
        if 'solvation' in self.data:
            self._process_solvation(self.data['solvation'])
        if 'liquidVolumetricMassTransferCoefficientPowerLaw' in self.data:
            self._process_volumetric_mass_transfer(self.data['liquidVolumetricMassTransferCoefficientPowerLaw'])
        if 'simulator' in self.data:
            self._process_simulator(self.data['simulator'])
        if 'model' in self.data:
            self._process_simulator(self.data['model'])
        if 'quantumMechanics' in self.data:
            self._process_quantum_mechanics(self.data['quantumMechanics'])
        if 'mlEstimator' in self.data:
            self._process_ml_estimator(self.data['mlEstimator'])
        if 'pressureDependence' in self.data:
            self._process_pressure_dependence(self.data['pressureDependence'])
        if 'generatedSpeciesConstraints' in self.data:
            self._process_species_constraints(self.data['generatedSpeciesConstraints'])
        if 'thermoCentralDatabase' in self.data:
            self._process_thermo_central_database(self.data['thermoCentralDatabase'])
        if 'uncertainty' in self.data:
            self._process_uncertainty(self.data['uncertainty'])
        if 'restartFromSeed' in self.data:
            self._process_restart_from_seed(self.data['restartFromSeed'])
        if 'options' in self.data:
            self._process_options(self.data['options'])
        

    # FOR PROCESSOR FUNCTIONS DO THIS:
    # FOR DATABASE PROCESSOR, JUST PASS THE INFO INTO THE OG DB FUNC
    def _process_database(self, db_data):
        """
        process database input
        """
        reaction_libraries = []
        if 'reactionLibraries' in db_data:
            for lib in db_data['reactionLibraries']:
                if isinstance(lib, str):
                    reaction_libraries.append(lib)
                elif isinstance(lib, dict):
                    # convert dict format to tuple format
                    name = lib.get('name')
                    seed = lib.get('seed', False) 
                    # if no seed bool set, default to False
                    reaction_libraries.append((name, seed))
        database(
            thermoLibraries = db_data.get('thermoLibraries'),
            transportLibraries = db_data.get('transportLibraries'),
            reactionLibraries = db_data.get('reactionLibraries'),
            frequenciesLibraries = db_data.get('frequenciesLibraries'),
            seedMechanisms = db_data.get('seedMechanisms'),
            kineticsFamilies = db_data.get('kineticsFamilies', 'default'),
            kineticsDepositories = db_data.get('kineticsDepositories', 'default'),
            kineticsEstimator = db_data.get('kineticsEstimator', 'rate rules'),
            adsorptionGroups = db_data.get('adsorptionGroups', 'adsorptionPt111')
        )
    # FOR SPECIES PROCESSORS, CALL ON THE EXISTING FUNCS WITH DATA FROM YAML FILE
    def _process_catalyst_properties(self, cat_data):
        """
        process catalyst input
        """
        catalyst_properties(
            bindingEnergies = cat_data.get('bindingEnergies'),
            surfaceSiteDensity = cat_data.get('surfaceSiteDensity'),
            metal = cat_data.get('metal'),
            coverageDependence = cat_data.get('coverageDependence', False) 
            # ^ if no coverage dependence bool is set, default to False like in input.py
        )
    def _process_species(self, spec_list):
        """
        process species definition
        """
        for spec in spec_list:
        # handle structure based on nested format or explicit type
            if 'structure' in spec:
                struc_data = spec['structure']

                if isinstance(struc_data, str):
                    structure = adjacency_list(struc_data)
                    # assumes explicit adjacency list if is only string
                elif isinstance(struc_data, dict):
                    # if not, check if its a dict w a key | [name]: [value]
                    if 'SMILES' in struc_data:
                        structure = smiles(struc_data['SMILES'])
                    elif 'InChI' in struc_data:
                        structure = inchi(struc_data['InChI'])
                    elif 'adjacencyList' in struc_data:
                        structure = adjacency_list(struc_data['adjacencyList'])
                    elif 'fragmentAdjacencyList' in struc_data:
                        structure = fragment_adj(struc_data['fragmentAdjacencyList'])
                    elif 'fragmentSMILES' in struc_data:
                        structure = fragment_smiles(struc_data['fragmentSMILES'])
                    else:
                        raise ValueError(f"Unknown structure format in species {spec.get('label', 'unknown')}")
                else: 
                    raise ValueError(f"Invalid structure format for forbidden {spec.get('label', 'unknown')}")
            else:
                raise ValueError(f"No structure provided for species {spec.get('label', 'unknown')}")

            species(
                label = spec['label'],
                structure = structure,
                reactive = spec.get('reactive', True),
                cut = spec.get('cut', False),
                size_threshold = spec.get('sizeThreshold')
            )
    def _process_forbidden(self, forb_list):
        """
        process forbidden structures
        similar method as did species
        """
        for forb in forb_list:
            if 'structure' in forb:
                struc_data = forb['structure']

                if isinstance(struc_data, str):
                    # assume adjacency list group
                    structure = adjacency_list_group(struc_data)
                elif isinstance(struc_data, dict):
                    if 'SMILES' in struc_data:
                        structure = smiles(struc_data['SMILES'])
                    elif 'SMARTS' in struc_data:
                        structure = smarts(struc_data['SMARTS'])
                    elif 'adjacencyList' in struc_data:
                        structure = adjacency_list(struc_data['adjacencyList'])
                    elif 'adjacencyListGroup' in struc_data:
                        structure = adjacency_list_group(struc_data['adjacencyListGroup'])
                    else:
                        raise ValueError(f"Unknown structure format in forbidden {forb.get('label', 'unknown')}")
                else:
                    raise ValueError(f"Invalid structure format for forbidden {forb.get('label', 'unknown')}")
            else:
                raise ValueError(f"No structure provided for species {forb.get('label', 'unknown')}")
            forbidden(
                label=forb['label'],
                structure=structure
            )



    # FOR REACTOR PROCESSOR, ALSO DO THE SAME 

    def _process_react(self, react_data):
        """
        process react specifications
        """
        react(react_data)
        
    def _process_simple_reactor(self, reactor_data):
        """
        process simple reactor configuration
        """
        # Handle both single reactor and list of reactors
        reactors = reactor_data if isinstance(reactor_data, list) else [reactor_data]
        
        for reactor in reactors:
            # Handle sensitivity which could be None, a string, or a list
            sensitivity = reactor.get('sensitivity')
            if sensitivity is None:
                sensitivity = None
            elif isinstance(sensitivity, str):
                sensitivity = [sensitivity]
            else:
                sensitivity = sensitivity
                
            simple_reactor(
                temperature=self._convert_quantity(reactor['temperature']),
                pressure=self._convert_quantity(reactor['pressure']),
                initialMoleFractions=reactor['initialMoleFractions'],
                nSims=reactor.get('nSims', 6),
                terminationConversion=reactor.get('terminationConversion'),
                terminationTime=self._convert_quantity(reactor.get('terminationTime')),
                terminationRateRatio=reactor.get('terminationRateRatio'),
                balanceSpecies=reactor.get('balanceSpecies'),
                sensitivity=sensitivity,
                sensitivityThreshold=reactor.get('sensitivityThreshold', 1e-3),
                sensitivityTemperature=self._convert_quantity(reactor.get('sensitivityTemperature')),
                sensitivityPressure=self._convert_quantity(reactor.get('sensitivityPressure')),
                sensitivityMoleFractions=reactor.get('sensitivityMoleFractions'),
                constantSpecies=reactor.get('constantSpecies')
            )
    
    def _process_constant_v_reactor(self, reactor_data):
        """
        process constant V ideal gas reactor configuration
        """
        reactors = reactor_data if isinstance(reactor_data, list) else [reactor_data]
        
        for reactor in reactors:
            constant_V_ideal_gas_reactor(
                temperature=self._convert_quantity(reactor['temperature']),
                pressure=self._convert_quantity(reactor['pressure']),
                initialMoleFractions=reactor['initialMoleFractions'],
                terminationConversion=reactor.get('terminationConversion'),
                terminationTime=self._convert_quantity(reactor.get('terminationTime')),
                terminationRateRatio=reactor.get('terminationRateRatio'),
                balanceSpecies=reactor.get('balanceSpecies')
            )
    
    def _process_constant_tp_reactor(self, reactor_data):
        """
        process constant T,P ideal gas reactor configuration
        """
        reactors = reactor_data if isinstance(reactor_data, list) else [reactor_data]
        
        for reactor in reactors:
            constant_TP_ideal_gas_reactor(
                temperature=self._convert_quantity(reactor['temperature']),
                pressure=self._convert_quantity(reactor['pressure']),
                initialMoleFractions=reactor['initialMoleFractions'],
                terminationConversion=reactor.get('terminationConversion'),
                terminationTime=self._convert_quantity(reactor.get('terminationTime')),
                terminationRateRatio=reactor.get('terminationRateRatio'),
                balanceSpecies=reactor.get('balanceSpecies')
            )
    
    def _process_liquid_cat_reactor(self, reactor_data):
        """
        process liquid catalyst reactor configuration
        """
        reactors = reactor_data if isinstance(reactor_data, list) else [reactor_data]
        
        for reactor in reactors:
            liquid_cat_reactor(
                temperature=self._convert_quantity(reactor['temperature']),
                initialConcentrations=self._convert_concentration_dict(reactor['initialConcentrations']),
                initialSurfaceCoverages=reactor['initialSurfaceCoverages'],
                surfaceVolumeRatio=self._convert_quantity(reactor['surfaceVolumeRatio']),
                distance=self._convert_quantity(reactor.get('distance')),
                viscosity=self._convert_quantity(reactor.get('viscosity')),
                surfPotential=self._convert_quantity(reactor.get('surfPotential')),
                liqPotential=self._convert_quantity(reactor.get('liqPotential')),
                terminationConversion=reactor.get('terminationConversion'),
                terminationTime=self._convert_quantity(reactor.get('terminationTime')),
                terminationRateRatio=reactor.get('terminationRateRatio'),
                constantSpecies=reactor.get('constantSpecies', [])
            )
    
    def _process_constant_tv_liquid_reactor(self, reactor_data):
        """
        process constant T,V liquid reactor configuration
        """
        reactors = reactor_data if isinstance(reactor_data, list) else [reactor_data]
        
        for reactor in reactors:
            constant_T_V_liquid_reactor(
                temperature=self._convert_quantity(reactor['temperature']),
                initialConcentrations=self._convert_concentration_dict(reactor['initialConcentrations']),
                liquidVolume=self._convert_quantity(reactor.get('liquidVolume')),
                residenceTime=self._convert_quantity(reactor.get('residenceTime')),
                inletVolumetricFlowRate=self._convert_quantity(reactor.get('inletVolumetricFlowRate')),
                outletVolumetricFlowRate=self._convert_quantity(reactor.get('outletVolumetricFlowRate')),
                inletConcentrations=self._convert_concentration_dict(reactor.get('inletConcentrations', {})),
                vaporPressure=self._convert_quantity(reactor.get('vaporPressure')),
                vaporMoleFractions=reactor.get('vaporMoleFractions'),
                terminationConversion=reactor.get('terminationConversion'),
                terminationTime=self._convert_quantity(reactor.get('terminationTime')),
                terminationRateRatio=reactor.get('terminationRateRatio'),
                constantSpecies=reactor.get('constantSpecies', [])
            )
    
    def _process_liquid_reactor(self, reactor_data):
        """
        process liquid reactor configuration
        """
        reactors = reactor_data if isinstance(reactor_data, list) else [reactor_data]
        
        for reactor in reactors:
            # Handle sensitivity which could be None, a string, or a list
            sensitivity = reactor.get('sensitivity')
            if sensitivity is None:
                sensitivity = None
            elif isinstance(sensitivity, str):
                sensitivity = [sensitivity]
            else:
                sensitivity = sensitivity
                
            liquid_reactor(
                temperature=self._convert_quantity(reactor['temperature']),
                initialConcentrations=self._convert_concentration_dict(reactor['initialConcentrations']),
                terminationConversion=reactor.get('terminationConversion'),
                nSims=reactor.get('nSims', 4),
                terminationTime=self._convert_quantity(reactor.get('terminationTime')),
                terminationRateRatio=reactor.get('terminationRateRatio'),
                sensitivity=sensitivity,
                sensitivityThreshold=reactor.get('sensitivityThreshold', 1e-3),
                sensitivityTemperature=self._convert_quantity(reactor.get('sensitivityTemperature')),
                sensitivityConcentrations=self._convert_concentration_dict(reactor.get('sensitivityConcentrations', {})),
                constantSpecies=reactor.get('constantSpecies')
            )
    
    def _process_surface_reactor(self, reactor_data):
        """
        process surface reactor configuration
        """
        reactors = reactor_data if isinstance(reactor_data, list) else [reactor_data]
        
        for reactor in reactors:
            # Handle sensitivity which could be None, a string, or a list
            sensitivity = reactor.get('sensitivity')
            if sensitivity is None:
                sensitivity = None
            elif isinstance(sensitivity, str):
                sensitivity = [sensitivity]
            else:
                sensitivity = sensitivity
                
            surface_reactor(
                temperature=self._convert_quantity(reactor['temperature']),
                initialPressure=self._convert_quantity(reactor['initialPressure']),
                initialGasMoleFractions=reactor['initialGasMoleFractions'],
                initialSurfaceCoverages=reactor['initialSurfaceCoverages'],
                surfaceVolumeRatio=self._convert_quantity(reactor['surfaceVolumeRatio']),
                nSims=reactor.get('nSims', 4),
                terminationConversion=reactor.get('terminationConversion'),
                terminationTime=self._convert_quantity(reactor.get('terminationTime')),
                terminationRateRatio=reactor.get('terminationRateRatio'),
                sensitivity=sensitivity,
                sensitivityThreshold=reactor.get('sensitivityThreshold', 1e-3)
            )
    
    def _process_mb_sampled_reactor(self, reactor_data):
        """
        process MB sampled reactor configuration
        """
        reactors = reactor_data if isinstance(reactor_data, list) else [reactor_data]
        
        for reactor in reactors:
            # Handle sensitivity which could be None, a string, or a list
            sensitivity = reactor.get('sensitivity')
            if sensitivity is None:
                sensitivity = None
            elif isinstance(sensitivity, str):
                sensitivity = [sensitivity]
            else:
                sensitivity = sensitivity
                
            mb_sampled_reactor(
                temperature=self._convert_quantity(reactor['temperature']),
                pressure=self._convert_quantity(reactor['pressure']),
                initialMoleFractions=reactor['initialMoleFractions'],
                mbsamplingRate=self._convert_quantity(reactor['mbsamplingRate']),
                terminationConversion=reactor.get('terminationConversion'),
                terminationTime=self._convert_quantity(reactor.get('terminationTime')),
                sensitivity=sensitivity,
                sensitivityThreshold=reactor.get('sensitivityThreshold', 1e-3),
                constantSpecies=reactor.get('constantSpecies')
            )
                
    def _process_solvation(self, solv_data):
        """
        process solvation settings
        """
        # Handle SolventData if provided
        
        # (FINISH THIS IDK HOW TO DO IT RN)
        
    def _process_liquid_mass_transfer(self, lmt_data):
        """
        process liquid volumetric mass transfer coefficient
        """
        liquid_volumetric_mass_transfer_coefficient_power_law(
            prefactor=self._convert_quantity(lmt_data.get('prefactor', (0, "1/s"))),
            diffusionCoefficientPower=lmt_data.get('diffusionCoefficientPower', 0),
            solventViscosityPower=lmt_data.get('solventViscosityPower', 0),
            solventDensityPower=lmt_data.get('solventDensityPower', 0)
        )
        
    def _process_simulator(self, sim_data):
        """
        process simulator settings
        """
        simulator(
            atol=sim_data.get('atol', 1e-16),
            rtol=sim_data.get('rtol', 1e-8),
            sens_atol=sim_data.get('sens_atol', 1e-6),
            sens_rtol=sim_data.get('sens_rtol', 1e-4)
        )
    
    def _process_model(self, model_data):
        """
        process model settings
        """
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
        
    def _process_quantum_mechanics(self, qm_data):
        """
        process quantum mechanics settings
        """
        quantum_mechanics(
            software=qm_data['software'],
            method=qm_data['method'],
            fileStore=qm_data.get('fileStore'),
            scratchDirectory=qm_data.get('scratchDirectory'),
            onlyCyclics=qm_data.get('onlyCyclics', False),
            maxRadicalNumber=qm_data.get('maxRadicalNumber', 0)
        )
        
    def _process_ml_estimator(self, ml_data):
        """
        process ML estimator settings
        """
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
        
    def _process_pressure_dependence(self, pd_data):
        """
        process pressure dependence settings
        """
        # Process temperatures - can be dict with min/max/count/units or list
        temps = pd_data['temperatures']
        if isinstance(temps, dict):
            temperatures = [temps['min'], temps['max'], temps['units'], temps['count']]
        else:
            temperatures = temps
            
        # Process pressures - can be dict with min/max/count/units or list
        press = pd_data['pressures']
        if isinstance(press, dict):
            pressures = [press['min'], press['max'], press['units'], press['count']]
        else:
            pressures = press
            
        # Process interpolation - can be list or tuple
        interp = pd_data.get('interpolation')
        if isinstance(interp, list) and len(interp) > 1:
            # Convert list format [method, param1, param2] to tuple
            interpolation = tuple(interp)
        else:
            interpolation = interp
            
        pressure_dependence(
            method=pd_data['method'],
            temperatures=temperatures,
            pressures=pressures,
            maximumGrainSize=self._convert_quantity(pd_data.get('maximumGrainSize', 0.0)),
            minimumNumberOfGrains=pd_data.get('minimumNumberOfGrains', 0),
            interpolation=interpolation,
            maximumAtoms=pd_data.get('maximumAtoms')
        )
        
    def _process_species_constraints(self, constraints):
        """
        process generated species constraints
        """
        # Create a copy to avoid modifying the original
        constraints_copy = constraints.copy()
        
        # Handle the special 'allowed' field
        if 'allowed' in constraints_copy:
            allowed_list = constraints_copy['allowed']
            # Convert special string values to their expected format
            processed_allowed = []
            for item in allowed_list:
                if item == 'input species':
                    processed_allowed.append('input species')
                elif item == 'seed mechanisms':
                    processed_allowed.append('seed mechanisms')
                elif item == 'reaction libraries':
                    processed_allowed.append('reaction libraries')
                else:
                    processed_allowed.append(item)
            constraints_copy['allowed'] = processed_allowed
            
        generated_species_constraints(**constraints_copy)
        
    def _process_thermo_central_database(self, tcd_data):
        """
        process thermo central database settings
        """
        thermo_central_database(
            host=tcd_data['host'],
            port=tcd_data['port'],
            username=tcd_data['username'],
            password=tcd_data['password'],
            application=tcd_data['application']
        )
        
    def _process_uncertainty(self, unc_data):
        """
        process uncertainty settings
        """
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
        
    def _process_restart_from_seed(self, restart_data):
        """
        process restart from seed settings
        """
        restart_from_seed(
            path=restart_data.get('path'),
            coreSeed=restart_data.get('coreSeed'),
            edgeSeed=restart_data.get('edgeSeed'),
            filters=restart_data.get('filters'),
            speciesMap=restart_data.get('speciesMap')
        )
        
    def _process_options(self, opt_data):
        """
        process general options
        """
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
            wallTime=opt_data.get('wallTime', '00:00:00:00'),  # This is a string, not a quantity
            saveSeedModulus=opt_data.get('saveSeedModulus', -1)
        )
        
    def _convert_quantity(self, value):
        """
        convert YAML quantity representation to tuple format expected by functions
        
        :param value: Either a dict with 'value' and 'units', a list/tuple, 
                      a single number/string, or None
        :return: Tuple (value, units), the original value, or None
        """
        if value is None:
            return None
            
        if isinstance(value, dict):
            if 'value' in value and 'units' in value:
                return (value['value'], value['units'])
            # Also handle the case where the quantity is directly the dict value
            elif len(value) == 1:
                # e.g., {0.5: 'kcal/mol'} format
                val, unit = next(iter(value.items()))
                return (val, unit)
            else:
                # Return the dict as-is if it doesn't match expected formats
                return value
        elif isinstance(value, (list, tuple)):
            if len(value) == 2:
                # Standard (value, units) format
                return tuple(value)
            elif len(value) == 4:
                # For temperature/pressure ranges in pressure dependence
                return value
            else:
                # Other list formats
                return value
        elif isinstance(value, str):
            # For cases like wallTime which is just a string
            return value
        else:
            # For single numeric values, return as-is
            # The function being called will handle unit defaults if needed
            return value
            
    def _convert_concentration_dict(self, conc_dict):
        """
        convert concentration dictionary with quantity values
        
        :param conc_dict: Dictionary with species as keys and quantities as values
        :return: Dictionary with converted quantities
        """
        if not conc_dict:
            return {}
            
        result = {}
        for species, conc in conc_dict.items():
            result[species] = self._convert_quantity(conc)
        return result
     
# Actual reader function itself now
def read_yaml_input_file(path, rmg0):
    """
    read an RMG YAML input file and process it using the existing input.py functions.
    
    :param path: Path to the YAML input file
    :param rmg0: RMG object to populate
    """
    # Import necessary modules for processing
    from rmgpy.rmg.input import set_global_rmg
    from rmgpy.rmg.model import CoreEdgeReactionModel
    
    # Set up the global RMG object
    set_global_rmg(rmg0)
    rmg0.reaction_model = CoreEdgeReactionModel()
    rmg0.initial_species = []
    rmg0.reaction_systems = []
    
    # Clear the global species_dict
    from rmgpy.rmg import input as rmg_input
    rmg_input.species_dict = {}
    rmg_input.mol_to_frag = {}
    
    # Set species constraints default
    rmg0.species_constraints = {'explicitlyAllowedMolecules': []}
    
    # Process YAML file
    reader = YAMLInputReader(path)
    reader.read()
    reader.process()
    
    # Post-processing (similar to original read_input_file)
    for reaction_system in rmg0.reaction_systems:
        if hasattr(reaction_system, 'convert_initial_keys_to_species_objects'):
            reaction_system.convert_initial_keys_to_species_objects(rmg_input.species_dict)
    
    if rmg0.quantum_mechanics:
        rmg0.quantum_mechanics.set_default_output_directory(rmg0.output_directory)
        rmg0.quantum_mechanics.initialize()
    
    logging.info('')

def read_input_file_wrapper(path, rmg0):
    """
    read an RMG input file (either Python or YAML format) and process it.
    
    this function automatically detects the file format based on the extension
    and calls the appropriate reader.
    
    :param path: Path to the input file (.py or .yaml/.yml)
    :param rmg0: RMG object to populate
    """
    import os
    from pathlib import Path
    
    # Get the file extension
    file_path = Path(path)
    extension = file_path.suffix.lower()
    
    # Check if file exists
    if not file_path.exists():
        raise IOError(f'The input file "{path}" could not be found.')
    
    # Route to appropriate reader based on extension
    if extension == '.py':
        # Use the original Python input file reader
        from rmgpy.rmg.input import read_input_file as read_python_input_file
        logging.info(f'Detected Python input file format (.py)')
        read_python_input_file(path, rmg0)
    elif extension in ['.yaml', '.yml']:
        # Use the YAML input file reader
        logging.info(f'Detected YAML input file format ({extension})')
        read_yaml_input_file(path, rmg0)
    else:
        raise ValueError(
            f'Unsupported input file format "{extension}". '
            f'RMG supports .py and .yaml/.yml input files.'
        )
            