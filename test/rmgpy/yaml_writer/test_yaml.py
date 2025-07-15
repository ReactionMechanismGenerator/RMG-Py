from compare_yaml_outputs import *
import pytest

@pytest.fixture
def compare_manager():
    '''Create instance of a Compare Yaml before each test'''
    yaml_files = {
        'yaml1': ['RMG_yaml_writer_addition/RMG-Py/test/rmgpy/test_data/yaml_writer_data/chemkin/', 'chem0047-gas.yaml'],
        'yaml2': ['RMG_yaml_writer_addition/RMG-Py/test/rmgpy/test_data/yaml_writer_data/cantera/', 'chem47.yaml']
    }
    return CompareYaml(yaml_files)

def test_compare_number_of_species(compare_manager):
    assert compare_manager.compare_species_count() == True

def test_compare_species_names(compare_manager):
    assert compare_manager.compare_species_names() == True

def test_compare_species_count_per_phase(compare_manager):
    assert compare_manager.compare_species_count_per_phase() == True

def test_compare_reactions(compare_manager):
    assert compare_manager.compare_reactions() == True