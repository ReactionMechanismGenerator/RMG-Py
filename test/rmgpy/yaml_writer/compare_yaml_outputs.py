import os
import yaml
import pandas as pd
import re

class YamlAnalyst:
    def __init__(self, path_to_yaml, yaml_file):
        self.path_to_yaml = path_to_yaml
        self.yaml_file = yaml_file

    def get_absolute_path(self):
        path = os.path.join(self.path_to_yaml, self.yaml_file)
        # If path is already absolute, use it; otherwise join with cwd
        if os.path.isabs(path):
            return path
        return os.path.join(os.getcwd(), path)

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
    Takes a dictionary with keys yaml1 and yaml2, and values a 
    list of the directory and file name of the yaml files.

    e.g.

    yaml_files = {
        'yaml1': [yaml1_file_directory, file1.yaml],
        'yaml2': [yaml2_file_directory, file2.yaml]
    }
    '''
    def __init__(self, yaml_files):
        self.yaml1 = YamlAnalyst(yaml_files['yaml1'][0], yaml_files['yaml1'][1])
        self.yaml2 = YamlAnalyst(yaml_files['yaml2'][0], yaml_files['yaml2'][1])

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
            components = side.split(' + ')
            normalized_components = []
            for component in components:
                # Remove any prefix integers/coefficients
                normalized_component = re.sub(r'^\d*\s*', '', component).strip()
                normalized_components.append(normalized_component)
            return ' + '.join(sorted(set(normalized_components)))

        reactants, products = equation.split('<=>')
        normalized_reactants = process_side(reactants)
        normalized_products = process_side(products)
        return f"{normalized_reactants} <=> {normalized_products}"

    def compare_reactions(self):
        reactions1 = self.yaml1.get_reaction_df()
        reactions2 = self.yaml2.get_reaction_df()
        comparison_results = {}

        for key1, df1 in reactions1.items():
            df1['normalized_equation'] = df1['equation'].apply(self.normalize_equation)
            for key2, df2 in reactions2.items():
                df2['normalized_equation'] = df2['equation'].apply(self.normalize_equation)
                merged_df = pd.merge(df1, df2, on='normalized_equation', suffixes=('_1', '_2'), how='inner')
                if not merged_df.empty:
                    merged_df['A_diff'] = merged_df['A_1'].round(2) - merged_df['A_2'].round(2)
                    merged_df['b_diff'] = merged_df['b_1'].round(2) - merged_df['b_2'].round(2)
                    merged_df['Ea_diff'] = merged_df['Ea_1'].round(2) - merged_df['Ea_2'].round(2)
                    comparison_results[f'{key1}_{key2}'] = merged_df[['normalized_equation', 'A_diff', 'b_diff', 'Ea_diff']]
        for key, df in comparison_results.items():
            if not (df['A_diff'].eq(0).all() and df['b_diff'].eq(0).all() and df['Ea_diff'].eq(0).all()):
                return False
        return True
