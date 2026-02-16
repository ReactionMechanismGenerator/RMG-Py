#!/usr/bin/env python3
"""
Script to compare two Cantera YAML files.

This script loads two Cantera mechanism YAML files and compares them for
structural and numerical differences. It compares:
- YAML metadata (generator, date, units, elements, phases structure)
- Species in each phase (names, ordering, and thermodynamic data)
- Reactions in each phase (equations, ordering, and kinetic data)

Usage:
    python compare_cantera_yaml.py <file1.yaml> <file2.yaml>
    python compare_cantera_yaml.py <file1.yaml> <file1.yaml>  # Should show no differences
"""

import sys
import argparse
import logging
from pathlib import Path
from typing import Any, List, Tuple, Dict
from itertools import chain

import yaml
import numpy as np
import cantera as ct


class CanteraModel:
    """
    Represents a Cantera mechanism model loaded from a YAML file.
    
    This class loads both the raw YAML data structure and creates Cantera
    Solution objects for each phase defined in the file.
    
    Attributes
    ----------
    file_path : str
        Path to the Cantera YAML file.
    yaml_data : dict
        Raw YAML data loaded from the file.
    phases : dict
        Dictionary of phase_name -> ct.Solution objects.
    """
    
    def __init__(self, file_path: str):
        """
        Initialize a CanteraModel by loading the YAML file and all its phases.
        
        Parameters
        ----------
        file_path : str
            Path to the Cantera YAML file.
            
        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        Exception
            If loading the YAML or phases fails.
        """
        self.file_path = file_path
        self.yaml_data = None
        self.phases = {}
        
        # Load the YAML structure
        self.load_yaml()
        
        # Extract phase names and load each phase
        phase_names = self._extract_phase_names()
        self.load_phases(phase_names)
    
    def load_yaml(self) -> dict:
        """
        Load the YAML file into a dictionary structure using yaml.safe_load().
        
        Returns
        -------
        dict
            The YAML data structure.
            
        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        yaml.YAMLError
            If YAML parsing fails.
        """
        with open(self.file_path, 'r') as f:
            self.yaml_data = yaml.safe_load(f)
        return self.yaml_data
    
    def _extract_phase_names(self) -> List[str]:
        """
        Extract phase names from the loaded YAML data.
        
        Returns
        -------
        list of str
            List of phase names found in the 'phases' section.
            
        Raises
        ------
        ValueError
            If 'phases' key is not found or is empty.
        """
        if self.yaml_data is None:
            raise ValueError("YAML data not loaded. Call load_yaml() first.")
        
        if 'phases' not in self.yaml_data:
            raise ValueError(f"No 'phases' key found in {self.file_path}")
        
        phases_list = self.yaml_data['phases']
        if not phases_list:
            raise ValueError(f"'phases' list is empty in {self.file_path}")
        
        phase_names = []
        for phase in phases_list:
            if 'name' not in phase:
                raise ValueError(f"Phase definition missing 'name' key: {phase}")
            phase_names.append(phase['name'])
        
        return phase_names
    
    def load_phases(self, phase_names: List[str]):
        """
        Load Cantera Solution objects for the specified phases.
        
        Parameters
        ----------
        phase_names : list of str
            List of phase names to load from the YAML file.
            
        Raises
        ------
        RuntimeError
            If loading any phase fails (with original exception chained).
        """
        for phase_name in phase_names:
            try:
                solution = ct.Solution(self.file_path, name=phase_name)
                self.phases[phase_name] = solution
            except Exception as e:
                raise RuntimeError(f"Failed to load phase '{phase_name}' from {self.file_path}") from e



def is_numeric(value: Any) -> bool:
    """Check if a value is numeric (int or float)."""
    return isinstance(value, (int, float, np.number)) and not isinstance(value, bool)


def compare_values(val1: Any, val2: Any, path: str, atol: float = 1e-12, 
                  rtol: float = 1e-3) -> List[str]:
    """Compare two values and return a list of differences.
    
    Parameters
    ----------
    val1 : Any
        First value to compare.
    val2 : Any
        Second value to compare.
    path : str
        Path in the dictionary (for error reporting).
    atol : float
        Absolute tolerance for numerical comparisons.
    rtol : float
        Relative tolerance for numerical comparisons.
        
    Returns
    -------
    list of str
        List of difference descriptions, empty if no differences found.
    """
    differences = []

    ### SPECIAL CASES
    # Special handling for 'elements' path - normalize to title case and sort
    if path.split('.')[-1] == 'elements':
        if isinstance(val1[0], str) and isinstance(val2[0], str):
            val1_normalized = sorted([v.title() for v in val1])
            val2_normalized = sorted([v.title() for v in val2])
            if val1_normalized != val2_normalized:
                differences.append(f"Elements list mismatch at {path}: {val1_normalized} vs {val2_normalized}")
            return differences
        if isinstance(val1[0], dict) and isinstance(val2[0], dict):
            for d in chain(val1, val2):
                d['symbol'] = d['symbol'].title()
            val1 = sorted([d for d in val1], key=lambda x: x['symbol'])
            val2 = sorted([d for d in val2], key=lambda x: x['symbol'])

    if path.endswith('Troe.T1') or path.endswith('Troe.T2') or path.endswith('Troe.T3'):
        rtol = 5e-3 # Relax tolerance due to rounding.
    
    ### END OF SPECIAL CASES
    
    # Type checking
    if type(val1) != type(val2):
        differences.append(f"Type mismatch at {path}: {type(val1).__name__} vs {type(val2).__name__}")
        return differences
    
    # Handle dictionaries recursively
    if isinstance(val1, dict):
        for key in set(list(val1.keys()) + list(val2.keys())):
            if key not in val1:
                if key == 'reference-pressure':
                    continue  # Ignore missing 'reference-pressure' key SPECIAL CASE
                differences.append(f"Missing key in first file at {path}.{key}")
            elif key not in val2:
                differences.append(f"Missing key in second file at {path}.{key}")
            else:
                new_path = f"{path}.{key}" if path else key
                differences.extend(compare_values(val1[key], val2[key], new_path, atol, rtol))
    
    # Handle lists recursively
    elif isinstance(val1, list):
        if len(val1) != len(val2):
            differences.append(f"List length mismatch at {path}: {len(val1)} vs {len(val2)}")
            # Compare up to the shorter length
            min_len = min(len(val1), len(val2))
        else:
            min_len = len(val1)
        
        for i in range(min_len):
            new_path = f"{path}[{i}]"
            differences.extend(compare_values(val1[i], val2[i], new_path, atol, rtol))
    
    # Handle numeric values with tolerance
    elif is_numeric(val1) and is_numeric(val2):
        # Use numpy.allclose for comparison
        if not np.isclose(val1, val2, atol=atol, rtol=rtol):
            # Compute the difference for reporting
            abs_diff = abs(val1 - val2)
            rel_diff = abs(abs_diff / val2) if val2 != 0 else float('inf')
            differences.append(f"Numerical difference at {path}: {val1} vs {val2} "
                             f"(abs_diff={abs_diff:.2e}, rel_diff={rel_diff:.2e})")
    
    # Handle strings and other comparable types
    elif val1 != val2:
        differences.append(f"Value mismatch at {path}: {val1!r} vs {val2!r}")
    
    return differences


def compare_yaml_files(file1: str, file2: str, atol: float = 1e-12,
                      rtol: float = 1e-3) -> List[str]:
    """Compare two Cantera YAML files.
    
    Compares both the raw YAML structure and the species/reactions
    from loaded Cantera phases.
    
    Parameters
    ----------
    file1 : str
        Path to the first YAML file.
    file2 : str
        Path to the second YAML file.
    atol : float
        Absolute tolerance for numerical comparisons.
    rtol : float
        Relative tolerance for numerical comparisons.
        
    Returns
    -------
    list of str
        List of difference descriptions (empty if files are equivalent).
    """
    differences = []
    
    logging.info("Loading %s...", file1)
    model1 = CanteraModel(file1)
    
    logging.info("Loading %s...", file2)
    model2 = CanteraModel(file2)
    
    # Compare YAML metadata (everything except species and reactions details)
    logging.info("Comparing YAML metadata...")
    yaml_meta1 = _extract_yaml_metadata(model1.yaml_data)
    yaml_meta2 = _extract_yaml_metadata(model2.yaml_data)

    for ym in (yaml_meta1, yaml_meta2):
        ym.pop('cantera-version', None)
        ym.pop('input-files', None)
        ym.pop('date', None)

    differences.extend(compare_values(yaml_meta1, yaml_meta2, "metadata", atol, rtol))
    
    # Compare phases sequentially by order
    logging.info("Comparing phases...")
    phase_list1 = model1.yaml_data.get('phases', [])
    phase_list2 = model2.yaml_data.get('phases', [])
    
    if len(phase_list1) != len(phase_list2):
        differences.append(f"Number of phases differs: {len(phase_list1)} vs {len(phase_list2)}")
    
    # Compare each phase by order
    for i in range(max(len(phase_list1), len(phase_list2))):
        if i >= len(phase_list1):
            differences.append(f"Phase {i}: missing in first file (second file has '{phase_list2[i]['name']}')")
            continue
        if i >= len(phase_list2):
            differences.append(f"Phase {i}: missing in second file (first file has '{phase_list1[i]['name']}')")
            continue
        
        phase1_name = phase_list1[i]['name']
        phase2_name = phase_list2[i]['name']
        
        if phase1_name != phase2_name:
            differences.append(f"Phase {i}: name differs: '{phase1_name}' vs '{phase2_name}'")
        
        # Compare species and reactions for this phase
        phase1 = model1.phases[phase1_name]
        phase2 = model2.phases[phase2_name]
        
        logging.info("  Comparing species in phase '%s'...", phase1_name)
        differences.extend(_compare_species(phase1, phase2, phase1_name, atol, rtol))
        
        logging.info("  Comparing reactions in phase '%s'...", phase1_name)
        differences.extend(_compare_reactions(phase1, phase2, phase1_name, atol, rtol))
    
    return differences


def compare_yaml_files_and_report(file1: str, file2: str, atol: float = 1e-12,
                                   rtol: float = 1e-3, output: str = None) -> bool:
    """Compare two Cantera YAML files and report results via logging.
    
    Performs a comparison between two Cantera YAML files and logs the results.
    If an output file path is provided, also writes the report to that file.
    
    Parameters
    ----------
    file1 : str
        Path to the first YAML file.
    file2 : str
        Path to the second YAML file.
    atol : float
        Absolute tolerance for numerical comparisons.
    rtol : float
        Relative tolerance for numerical comparisons.
    output : str, optional
        Path to an output file where the comparison report will be written.
        If None, only logs to the standard logger.
        
    Returns
    -------
    bool
        True if files are equivalent (no differences), False otherwise.
    """
    file_handler = None
    root_logger = logging.getLogger()
    
    try:
        # Set up optional file logging if output path provided
        if output:
            file_handler = logging.FileHandler(output)
            file_handler.setFormatter(logging.Formatter("%(message)s"))
            root_logger.addHandler(file_handler)
        
        # Check if file2 exists
        if not Path(file2).exists():
            logging.warning("Cantera YAML comparison skipped; file not found at %s", file2)
            return False
        
        # Perform the comparison
        differences = compare_yaml_files(file1, file2, atol, rtol)
        
        # Log and report results
        if differences:
            logging.warning(
                "Cantera YAML comparison found %d difference(s) between %s and %s",
                len(differences),
                file1,
                file2,
            )
            for diff in differences:
                logging.warning("  %s", diff)
            return False
        else:
            logging.info(
                "Cantera YAML comparison passed: %s matches %s",
                file1,
                file2,
            )
            return True
    
    except Exception:
        logging.exception("Cantera YAML comparison failed for %s vs %s", file1, file2)
        return False
    
    finally:
        # Clean up the file handler
        if file_handler:
            file_handler.close()
            root_logger.removeHandler(file_handler)


def _extract_yaml_metadata(yaml_data: dict) -> dict:
    """Extract metadata from YAML (excluding detailed species/reactions)."""
    metadata = yaml_data.copy()
    # Remove the detailed species and reactions lists since we'll compare those separately
    metadata.pop('species', None)
    metadata.pop('reactions', None)
    reaction_blocks = []
    for phase in yaml_data.get('phases', []):
        reactions = phase.get('reactions', [])
        if reactions in ('declared-species', 'all', 'none'):
            continue
        reaction_blocks.extend(reactions)
    for block in reaction_blocks:
        if block not in metadata:
            raise ValueError(f"Phase mentioned reactions block '{block}' not found in top-level YAML keys")
        metadata.pop(block)
    return metadata


def _compare_species(phase1, phase2, phase_name: str, atol: float, rtol: float) -> List[str]:
    """Compare species in two Cantera phases."""
    differences = []
    
    species1 = phase1.species()
    species2 = phase2.species()
    
    if len(species1) != len(species2):
        differences.append(f"Phase '{phase_name}': number of species differs: {len(species1)} vs {len(species2)}")
    
    # Build name-to-species mapping for comparison
    species1_map = {sp.name: sp for sp in species1}
    species2_map = {sp.name: sp for sp in species2}
    
    # Check for missing species
    names1 = set(species1_map.keys())
    names2 = set(species2_map.keys())
    
    for name in names1 - names2:
        differences.append(f"Phase '{phase_name}': species '{name}' only in first file")
    for name in names2 - names1:
        differences.append(f"Phase '{phase_name}': species '{name}' only in second file")
    
    # Compare ordering
    common_species = names1 & names2
    for i in range(min(len(species1), len(species2))):
        if species1[i].name != species2[i].name:
            if species1[i].name in common_species and species2[i].name in common_species:
                differences.append(f"Phase '{phase_name}': species order differs at index {i}: "
                                 f"'{species1[i].name}' vs '{species2[i].name}'")
            break  # Only report first ordering difference
    
    # Compare input_data for matching species
    for name in sorted(common_species):
        sp1 = species1_map[name]
        sp2 = species2_map[name]
        path = f"phase.{phase_name}.species.{name}"
        differences.extend(compare_values(sp1.input_data, sp2.input_data, path, atol, rtol))
    
    return differences


def _compare_reactions(phase1, phase2, phase_name: str, atol: float, rtol: float) -> List[str]:
    """Compare reactions in two Cantera phases."""
    differences = []
    
    reactions1 = phase1.reactions()
    reactions2 = phase2.reactions()
    
    if len(reactions1) != len(reactions2):
        differences.append(f"Phase '{phase_name}': number of reactions differs: {len(reactions1)} vs {len(reactions2)}")
    
    # Compare reactions by index (assuming they should be in the same order)
    for i in range(min(len(reactions1), len(reactions2))):
        rxn1 = reactions1[i]
        rxn2 = reactions2[i]
        
        # Check if equations match
        eq1 = rxn1.equation
        eq2 = rxn2.equation
        
        if eq1 != eq2:
            differences.append(f"Phase '{phase_name}': reaction {i} equation differs: '{eq1}' vs '{eq2}'")
            # Still compare input_data even if equations differ
        
        path = f"phase.{phase_name}.reaction[{i}]"
        differences.extend(compare_values(rxn1.input_data, rxn2.input_data, path, atol, rtol))
    
    return differences


def main():
    """Main entry point for the comparison script."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    parser = argparse.ArgumentParser(
        description="Compare two Cantera YAML mechanism files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python compare_cantera_yaml.py file1.yaml file2.yaml
  python compare_cantera_yaml.py file1.yaml file1.yaml  # Should show no differences
  python compare_cantera_yaml.py --abs-tol 1e-6 file1.yaml file2.yaml
        """
    )
    parser.add_argument("file1", help="First Cantera YAML file")
    parser.add_argument("file2", help="Second Cantera YAML file")
    parser.add_argument("--abs-tol", type=float, default=1e-11,
                       help="Absolute tolerance for numerical comparisons (default: 1e-9)")
    parser.add_argument("--rel-tol", type=float, default=1e-3,
                       help="Relative tolerance for numerical comparisons (default: 1e-9)")
    
    args = parser.parse_args()
    
    # Verify files exist
    for file_path in [args.file1, args.file2]:
        if not Path(file_path).exists():
            logging.error("File not found: %s", file_path)
            sys.exit(1)
    
    try:
        differences = compare_yaml_files(
            args.file1, args.file2, args.abs_tol, args.rel_tol
        )
        
        if len(differences) == 0:
            logging.info("Files are equivalent (within specified tolerances)")
            sys.exit(0)
        else:
            logging.warning("Files differ. Found %d difference(s):", len(differences))
            for i, diff in enumerate(differences, 1):
                logging.info("%3d. %s", i, diff)
            sys.exit(1)
            
    except Exception as e:
        logging.exception("Error: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    if len(sys.argv) == 1:
        logging.info("No arguments provided. Using default test files for demonstration.")
        # sys.argv.extend([
        #     "test/rmgpy/test_data/yaml_writer_data/chemkin/from_main_test.yaml",
        #     "test/rmgpy/test_data/yaml_writer_data/cantera/from_main_test.yaml"
        # ])

        sys.argv.extend([
            "/Users/rwest/Code/RMG-Py/testing/eg0/cantera_from_ck/chem.yaml",
            "/Users/rwest/Code/RMG-Py/testing/eg0/cantera2/chem.yaml"
        ])

    main()
