#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

from rmgpy.util import strip_yaml_notes


class UtilTest:
    def strip_yaml_notes(self, tmp_path, source_text):
        source_path = tmp_path / "chem_annotated.yaml"
        destination_path = tmp_path / "chem.yaml"
        source_path.write_text(source_text)

        strip_yaml_notes(source_path, destination_path)

        return destination_path.read_text()

    def test_strip_yaml_notes_removes_block_style_note(self, tmp_path):
        source = """species:
- name: CH4
  composition: {C: 1, H: 4}
  note: RMG-generated species note
- name: O2
  composition: {O: 2}
"""
        expected = """species:
- name: CH4
  composition: {C: 1, H: 4}
- name: O2
  composition: {O: 2}
"""

        assert self.strip_yaml_notes(tmp_path, source) == expected

    def test_strip_yaml_notes_removes_block_style_multiline_note(self, tmp_path):
        source = """reactions:
- equation: CH4 + O2 <=> CH3 + HO2
  rate-constant: {A: 1.0e+06, b: 0.0, Ea: 10000.0}
  note: |
    Estimated by RMG.
    Includes database comments.
- equation: H + O2 <=> O + OH
  rate-constant: {A: 1.0e+14, b: 0.0, Ea: 15000.0}
"""
        expected = """reactions:
- equation: CH4 + O2 <=> CH3 + HO2
  rate-constant: {A: 1.0e+06, b: 0.0, Ea: 10000.0}
- equation: H + O2 <=> O + OH
  rate-constant: {A: 1.0e+14, b: 0.0, Ea: 15000.0}
"""

        assert self.strip_yaml_notes(tmp_path, source) == expected

    def test_strip_yaml_notes_removes_single_line_flow_note(self, tmp_path):
        source = """species:
- name: Ar
  transport: {model: gas, note: RMG transport, geometry: atom, diameter: 3.33, well-depth: 136.5}
"""
        expected = """species:
- name: Ar
  transport: {model: gas, geometry: atom, diameter: 3.33, well-depth: 136.5}
"""

        assert self.strip_yaml_notes(tmp_path, source) == expected

    def test_strip_yaml_notes_removes_wrapped_flow_note(self, tmp_path):
        source = """species:
- name: CH4
  transport: {model: gas, geometry: nonlinear, diameter: 3.746,
    note: RMG transport note,
    well-depth: 141.4}
"""
        expected = """species:
- name: CH4
  transport: {model: gas, geometry: nonlinear, diameter: 3.746,
    well-depth: 141.4}
"""

        assert self.strip_yaml_notes(tmp_path, source) == expected

    def test_strip_yaml_notes_removes_wrapped_flow_multiline_note(self, tmp_path):
        source = """species:
- name: CH4
  transport: {model: gas, geometry: nonlinear, diameter: 3.746,
    note: RMG transport note
      with wrapped detail,
    well-depth: 141.4}
"""
        expected = """species:
- name: CH4
  transport: {model: gas, geometry: nonlinear, diameter: 3.746,
    well-depth: 141.4}
"""

        assert self.strip_yaml_notes(tmp_path, source) == expected
