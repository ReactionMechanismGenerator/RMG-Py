#!/usr/bin/env python3
"""
Script to compare two Cantera YAML files.

This script loads two Cantera mechanism YAML files using cantera.Solution,
extracts their input_data dictionaries, and compares them for structural and
numerical differences.

Usage:
    python compare_cantera_yaml.py <file1.yaml> <file2.yaml>
    python compare_cantera_yaml.py <file1.yaml> <file1.yaml>  # Should show no differences
"""

import sys
import argparse
from pathlib import Path
from typing import Any, List, Tuple

import numpy as np
import cantera as ct


def load_cantera_input_data(yaml_file: str) -> dict:
    """Load a Cantera YAML file and return its input_data dictionary.
    
    Parameters
    ----------
    yaml_file : str
        Path to the Cantera YAML file.
        
    Returns
    -------
    dict
        The input_data dictionary from ct.Solution.
        
    Raises
    ------
    Exception
        If the file cannot be loaded.
    """
    try:
        solution = ct.Solution(yaml_file)
        return solution.input_data
    except Exception as e:
        raise Exception(f"Failed to load {yaml_file}: {e}")


def is_numeric(value: Any) -> bool:
    """Check if a value is numeric (int or float)."""
    return isinstance(value, (int, float, np.number)) and not isinstance(value, bool)


def compare_values(val1: Any, val2: Any, path: str, atol: float = 1e-9, 
                  rtol: float = 1e-9) -> List[str]:
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
    print(path)  # Debug: print the current path being compared
    # Special handling for 'elements' path - normalize to title case and sort
    if path == 'elements' and isinstance(val1, list) and isinstance(val2, list):
        val1_normalized = sorted([str(v).title() for v in val1])
        val2_normalized = sorted([str(v).title() for v in val2])
        if val1_normalized != val2_normalized:
            differences.append(f"Elements list mismatch at {path}: {val1_normalized} vs {val2_normalized}")
        return differences
    
    # Type checking
    if type(val1) != type(val2):
        differences.append(f"Type mismatch at {path}: {type(val1).__name__} vs {type(val2).__name__}")
        return differences
    
    # Handle dictionaries recursively
    if isinstance(val1, dict):
        for key in set(list(val1.keys()) + list(val2.keys())):
            if key not in val1:
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


def compare_yaml_files(file1: str, file2: str, atol: float = 1e-9,
                      rtol: float = 1e-9) -> List[str]:
    """Compare two Cantera YAML files.
    
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
    print(f"Loading {file1}...")
    data1 = load_cantera_input_data(file1)
    
    print(f"Loading {file2}...")
    data2 = load_cantera_input_data(file2)
    
    print("Comparing files...")
    return compare_values(data1, data2, "", atol, rtol)


def main():
    """Main entry point for the comparison script."""
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
    parser.add_argument("--abs-tol", type=float, default=1e-9,
                       help="Absolute tolerance for numerical comparisons (default: 1e-9)")
    parser.add_argument("--rel-tol", type=float, default=1e-9,
                       help="Relative tolerance for numerical comparisons (default: 1e-9)")
    
    args = parser.parse_args()
    
    # Verify files exist
    for file_path in [args.file1, args.file2]:
        if not Path(file_path).exists():
            print(f"Error: File not found: {file_path}", file=sys.stderr)
            sys.exit(1)
    
    try:
        differences = compare_yaml_files(
            args.file1, args.file2, args.abs_tol, args.rel_tol
        )
        
        print("\n" + "="*70)
        if len(differences) == 0:
            print("✓ Files are equivalent (within specified tolerances)")
            sys.exit(0)
        else:
            print(f"✗ Files differ. Found {len(differences)} difference(s):\n")
            for i, diff in enumerate(differences, 1):
                print(f"{i:3d}. {diff}")
            sys.exit(1)
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
