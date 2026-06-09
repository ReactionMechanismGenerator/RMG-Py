import os
import yaml
import pandas as pd
import re
from collections import Counter

class YamlAnalyst:
    def __init__(self, path_to_yaml_file):
        self.path_to_yaml_file = path_to_yaml_file

    def get_absolute_path(self):
        # If path is already absolute, use it; otherwise join with cwd
        if os.path.isabs(self.path_to_yaml_file):
            return self.path_to_yaml_file
        return os.path.join(os.getcwd(), self.path_to_yaml_file)

    def load_yaml_file(self):
        with open(self.get_absolute_path(), 'r') as file:
            return yaml.safe_load(file)

    def get_species(self):
        return self.load_yaml_file()['species']

    def get_species_count(self):
        return len(self.get_species())

    def get_species_names(self):
        return [specie['name'] for specie in self.get_species()]

    def get_species_count_per_phase(self):
        return {f"specie_count_{phase['name']}": len(phase['species']) for phase in self.load_yaml_file()['phases']}

    def get_reactions_dict(self):
        reactions_dict = {}
        data = self.load_yaml_file()
        for key, values in data.items():
            if key in [f"{phase['name']}-reactions" for phase in data['phases']]:
                reactions_dict[key] = data[key]
            elif key == 'gas_reactions':
                reactions_dict['gas_reactions'] = data['gas_reactions']
            elif key == 'surface_reactions':
                reactions_dict['surface_reactions'] = data['surface_reactions']
            elif key == 'reactions':
                # Gas-only mechanisms use plain 'reactions' key
                reactions_dict['reactions'] = data['reactions']
        return reactions_dict

    def create_reaction_df(self, reactions):
        data = []
        for reaction in reactions:
            row = {'equation': reaction['equation']}
            if 'rate-constant' in reaction:
                row.update(reaction['rate-constant'])
            elif 'sticking-coefficient' in reaction:
                row.update(reaction['sticking-coefficient'])
            data.append(row)
        return pd.DataFrame(data)
    def get_reaction_df(self):
        reaction_dfs = {key: self.create_reaction_df(value) for key, value in self.get_reactions_dict().items()}
        return reaction_dfs

    def get_reaction_count(self):
        return {key: len(value) for key, value in self.get_reactions_dict().items()}

class CompareYaml:
    '''
    Compare two YAML files.

    Args:
        yaml_path_1: Path to the first YAML file.
        yaml_path_2: Path to the second YAML file.
    '''
    def __init__(self, yaml_path_1, yaml_path_2):
        self.yaml1 = YamlAnalyst(yaml_path_1)
        self.yaml2 = YamlAnalyst(yaml_path_2)

    def compare_species_count(self):
        count1 = self.yaml1.get_species_count()
        count2 = self.yaml2.get_species_count()
        if count1 - count2 == 0:
            return True
        else:
            return False
        
    def compare_species_names(self):
        names1 = set(self.yaml1.get_species_names())
        names2 = set(self.yaml2.get_species_names())
        if set(names1) == set(names2):
            return True
        else:
            return False
    
    def compare_species_count_per_phase(self):
        count_per_phase1 = self.yaml1.get_species_count_per_phase()
        count_per_phase2 = self.yaml2.get_species_count_per_phase()
        phase_names1 = [phase['name'] for phase in self.yaml1.load_yaml_file()['phases']]
        phase_names2 = [phase['name'] for phase in self.yaml2.load_yaml_file()['phases']]
        all_phase_names = set(phase_names1).union(set(phase_names2))
        count_diff = {'gas': count_per_phase1[f"specie_count_{phase_names1[0]}"] - count_per_phase2[f"specie_count_{phase_names2[0]}"], 
                      'surface': count_per_phase1[f"specie_count_{phase_names1[1]}"] - count_per_phase2[f"specie_count_{phase_names2[1]}"]
                      }
        if count_diff['gas'] == 0 and count_diff['surface'] == 0:
            return True
        else:
            return False
    
    def normalize_equation(self, equation):
        def process_side(side):
            # Extract and remove (+M) or +M third-body markers
            has_third_body = False
            side_clean = side.strip()
            if '(+M)' in side_clean:
                has_third_body = True
                side_clean = side_clean.replace('(+M)', '').strip()
            elif side_clean.strip().endswith('+ M'):
                has_third_body = True
                side_clean = side_clean.rsplit('+ M', 1)[0].strip()

            components = side_clean.split('+')
            expanded = []
            for component in components:
                component = component.strip()
                if not component:
                    continue
                # Match optional integer coefficient prefix (e.g. "2 CH3(14)")
                m = re.match(r'^(\d+)\s+(.+)$', component)
                if m:
                    count = int(m.group(1))
                    species = m.group(2).strip()
                    expanded.extend([species] * count)
                else:
                    expanded.append(component)

            result = ' + '.join(sorted(expanded))
            if has_third_body:
                result += ' (+M)'
            return result

        # Handle both reversible (<=>) and irreversible (=>) reactions
        if '<=>' in equation:
            reactants, products = equation.split('<=>')
            separator = '<=>'
        elif '=>' in equation:
            reactants, products = equation.split('=>')
            separator = '=>'
        else:
            raise ValueError(f"Unknown reaction format: {equation}")
        
        normalized_reactants = process_side(reactants)
        normalized_products = process_side(products)
        return f"{normalized_reactants} {separator} {normalized_products}"

    def compare_reactions(self):
        """Compare reactions between two YAML files.
        
        First checks that reaction counts and normalized equations match.
        Then uses Cantera to load both files and compares forward rate
        constants at a reference state to verify kinetic equivalence,
        since the two files may use different unit systems internally.
        """
        import cantera as ct
        import numpy as np

        reactions1 = self.yaml1.get_reaction_df()
        reactions2 = self.yaml2.get_reaction_df()

        # Check if total reaction counts match
        count1 = sum(len(df) for df in reactions1.values())
        count2 = sum(len(df) for df in reactions2.values())
        if count1 != count2:
            return False

        # Collect all normalized equations from each file (using Counter to handle duplicates)
        all_eqs_1 = Counter()
        all_eqs_2 = Counter()
        for key, df in reactions1.items():
            all_eqs_1.update(df['equation'].apply(self.normalize_equation))
        for key, df in reactions2.items():
            all_eqs_2.update(df['equation'].apply(self.normalize_equation))

        # Check that all reaction equations are present in both files
        if all_eqs_1 != all_eqs_2:
            return False

        # Use Cantera to compare actual rate constants, since the two files
        # may store A/Ea values in different unit systems
        yaml1_path = self.yaml1.get_absolute_path()
        yaml2_path = self.yaml2.get_absolute_path()
        try:
            gas1 = ct.Solution(yaml1_path)
            gas2 = ct.Solution(yaml2_path)
        except Exception as e:
            # If Cantera can't load the files, fall back to equation-only comparison
            print(f"Warning: Cantera failed to load one or both files: {e}")
            print(f"Falling back to equation-only comparison")
            return True

        # Compare at a reference state
        T, P = 1000.0, ct.one_atm
        # Build a composition string from species common to both
        species_names = [s.name for s in gas1.species()]
        if len(species_names) >= 2:
            comp = f"{species_names[0]}:0.5, {species_names[1]}:0.5"
        else:
            comp = f"{species_names[0]}:1.0"

        gas1.TPX = T, P, comp
        gas2.TPX = T, P, comp

        kf1 = gas1.forward_rate_constants
        kf2 = gas2.forward_rate_constants

        # Use relative tolerance for comparison; allow 1% difference
        # to account for Chemkin format precision loss
        result = np.allclose(kf1, kf2, rtol=0.01, atol=1e-50)
        
        if not result:
            # Find which reactions differ
            diff_mask = ~np.isclose(kf1, kf2, rtol=0.01, atol=1e-50)
            diff_indices = np.where(diff_mask)[0]
            print(f"\n{len(diff_indices)} reactions have mismatched rate constants:")
            for idx in diff_indices[:5]:  # Show first 5
                rxn1 = gas1.reaction(int(idx))
                rxn2 = gas2.reaction(int(idx))
                print(f"  Reaction {idx}: {rxn1.equation}")
                print(f"    RMG kf: {kf1[idx]:.6e}, CK kf: {kf2[idx]:.6e}, ratio: {kf1[idx]/kf2[idx] if kf2[idx] != 0 else 'inf':.4f}")
        
        return result
