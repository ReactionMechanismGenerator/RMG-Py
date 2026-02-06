from compare_yaml_outputs import *
import os
import pytest

@pytest.fixture(scope="module")
def compare_manager():
    '''Create instance of a Compare Yaml before each test'''
    test_data_folder='test/rmgpy/test_data/yaml_writer_data/chemkin/'
    yaml_path_1 = os.path.join(test_data_folder, 'chem0047-gas.yaml')
    yaml_path_2 = os.path.join(test_data_folder, 'chem47.yaml')
    return CompareYaml(yaml_path_1, yaml_path_2)

def test_compare_number_of_species(compare_manager):
    assert compare_manager.compare_species_count() == True

def test_compare_species_names(compare_manager):
    assert compare_manager.compare_species_names() == True

def test_compare_species_count_per_phase(compare_manager):
    assert compare_manager.compare_species_count_per_phase() == True

def test_compare_reactions(compare_manager):
    assert compare_manager.compare_reactions() == True