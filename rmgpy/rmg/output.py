#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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

"""
This module contains functionality for saving the output of RMG jobs to output
files.
"""

import logging
import os.path
import re
import textwrap

from rmgpy.chemkin import get_species_identifier
from rmgpy.exceptions import OutputError
from rmgpy.util import make_output_subdirectory


################################################################################

def save_output_html(path, reaction_model, part_core_edge='core'):
    """
    Save the current set of  species and reactions of `reactionModel` to
    an HTML file `path` on disk. As part of this process, drawings of all 
    species are created in the species folder (if they don't already exist)
    using the :mod:`rmgpy.molecule.draw` module. The :mod:`jinja`
    package is used to generate the HTML; if this package is not found, no
    HTML will be generated (but the program will carry on).
    """

    from rmgpy.rmg.model import PDepReaction

    from rmgpy.molecule.draw import MoleculeDrawer

    try:
        import jinja2
    except ImportError:
        logging.warning("jinja2 package not found; HTML output will not be saved.")
        return

    path = os.path.abspath(path)
    dirname = os.path.dirname(path)

    # Prepare parameters to pass to jinja template
    title = 'RMG Output'

    if part_core_edge == 'core':
        species = reaction_model.core.species[:] + reaction_model.output_species_list
    elif part_core_edge == 'edge':
        species = reaction_model.edge.species[:] + reaction_model.output_species_list

    if not os.path.isdir(os.path.join(dirname, 'species')):
        os.makedirs(os.path.join(dirname, 'species'))

    re_index_search = re.compile(r'\((\d+)\)$').search

    for spec in species:
        # if the species dictionary came from an RMG-Java job, make them prettier
        # We use the presence of a trailing index on the label to discern this
        # (A single open parenthesis is not enough (e.g. when using SMILES strings as labels!)
        match = re_index_search(spec.label)
        if match:
            spec.index = int(match.group(0)[1:-1])
            spec.label = spec.label[0:match.start()]
        # Draw molecules if necessary
        fstr = os.path.join(dirname, 'species', '{0}.png'.format(spec))
        if not os.path.exists(fstr):
            try:
                MoleculeDrawer().draw(spec.molecule[0], 'png', fstr)
            except IndexError:
                logging.error("{0} species could not be drawn because it did not contain a molecular structure. "
                              "Please recheck your files.".format(get_species_identifier(spec)))
                raise
        # spec.thermo.comment=
        # Text wrap the thermo comments
    # We want to keep species sorted in the original order in which they were added to the RMG core.
    # Rather than ordered by index
    #    species.sort(key=lambda x: x.index)

    if part_core_edge == 'core':
        reactions = [rxn for rxn in reaction_model.core.reactions] + reaction_model.output_reaction_list
    elif part_core_edge == 'edge':
        reactions = [rxn for rxn in reaction_model.edge.reactions] + reaction_model.output_reaction_list

    # We want to keep reactions sorted in original order in which they were added to core
    # rather than ordered by index
    # reactions.sort(key=lambda x: x.index)

    family_count = {}
    for rxn in reactions:

        if isinstance(rxn, PDepReaction):
            family = "PDepNetwork"
        else:
            family = rxn.get_source()
        if family in family_count:
            family_count[family] += 1
        else:
            family_count[family] = 1
    families = list(family_count.keys())
    families.sort()

    ## jinja2 filters etc.
    to_remove_from_css_names = re.compile(r'[/.\-+,]')

    def csssafe(input):
        """Replace unsafe CSS class name characters with an underscore."""
        return to_remove_from_css_names.sub('_', input)

    environment = jinja2.Environment()
    environment.filters['csssafe'] = csssafe

    # Make HTML file
    template = environment.from_string(
"""<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN">
<html lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html;charset=utf-8" >
<title>{{ title }}</title>
<style type="text/css">
    body {
        font-family: sans-serif;
    }
    a {
        color: #993333;
        text-decoration: none;
    }
    a:visited {
        color: #993333;
    }
    a:hover {
        text-decoration: underline;
    }

    
    table.speciesList, table.reactionList {
        border-collapse: collapse;
        align: center;
    }

        
    table.speciesList th, table.reactionList th {
        text-align: left;
        vertical-align: top;
    }
    table.reaction {
        border-top: 1px solid #808080;        
        padding: 10px;
    }
    td.reactants {
        text-align: right;
    }
    td.products {
        text-align: left;
    }
    td.reactionArrow {
        text-align: center;
        font-size: 16px;
    }
    td.species img, td.reactants img, td.products img {
        vertical-align: middle;
    }
    
    img.surface_species {
        vertical-align: bottom;
    }
    
    tr.species{
        border-top: 1px solid #808080;        
    }

    tr.rxnStart{
        border-top: 1px solid #808080;        
    }
    
    td, .speciesList th{
        padding: 10px;
        }
    
    td.index {
    width: 50px;
    }
    
    .thermo td, .thermo th{
    padding: 2px;
    }
    
    .kinetics td,  .kinetics th {
    padding: 2px;
    }
    
    tr.kinetics {
        font-size: small;
    }
    .KineticsData {
        # border: 1px solid gray;
    }
    .KineticsData th {
        width: 15em;
        word-wrap: none;
    }
    .KineticsData td {
        width: 3em;
    }
    
    .chemkin, .KineticsData_repr {
       white-space: pre-wrap;
       font-size: x-small;
       font-family: "Andale Mono", monospace;
    }
    
    .energy {
        font-size: small;
    }
    
    .thermoComment {
       white-space: pre-wrap;
       font-size: small;
       font-family: "Andale Mono", monospace;
    }
    .hide_kinetics .kinetics{
       display: none !important;
    }
    .hide_chemkin .chemkin{
       display: none !important;
    }
    .hide_reaction{
       display: none !important;
    }
    .hide_energy .energy{
       display: none !important;
    }
    .hide_thermoComment .thermoComment{
       display: none !important;
    }

           
</style>
<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.4.1/jquery.min.js"></script>
<script type="text/javascript" src="../../../external/jquery.min.js"></script>
<script type="text/javascript">
function updateFamily(family) {
    if (family.checked) {
        $("."+family.value).show();
    } else {
        $("."+family.value).hide();
    }
}
function updateDetails(type) {
    if (type.checked) {
        $("."+type.value).show();
    } else {
        $("."+type.value).hide();
    }
}
function checkAllFamilies() {
    $("#familySelector").find("[name='family']").each(function() { this.checked = true; updateFamily(this); });
    return false;
}
function uncheckAllFamilies() {
    $("#familySelector").find("[name='family']").each(function() { this.checked = false; updateFamily(this); });
    return false;
}
function checkAllDetails() {
    $("#familySelector").find("[name='detail']").each(function() { this.checked = true; updateDetails(this); });
    return false;
}
function uncheckAllDetails() {
    $("#familySelector").find("[name='detail']").each(function() { this.checked = false; updateDetails(this); });
    return false;
}

function updateFamilyDetails() {
    $("#familySelector").find("[name='family']").each(function() { 
        updateDetails(this); });
    return false;
}
function updateReactionDetails() {
    $("#familySelector").find("[name='detail']").each(function() { 
        updateDetails(this); });
    return false;
}

function updateThermoDetails(type) {
    if (type.checked) {
        $(".thermoComment").removeClass("hide_"+type.value);
    } else {
        $(".thermoComment").addClass("hide_"+type.value);
    }
}

function uncheckThermoDetails() {
    $("#thermoSelector").find("[name='detail']").each(function() { this.checked = false; updateThermoDetails(this); });
    return false;
}

function resetReactionFilter() {
    $.each($(".reaction"), function() {
        $(this).removeClass("hide_reaction");
    });
}

function submitReactionFilter(){
    resetReactionFilter();
    _r1 = $("#reactant1").val().toLowerCase();
    _r2 = $("#reactant2").val().toLowerCase();
    _p1 = $("#product1").val().toLowerCase();
    _p2 = $("#product2").val().toLowerCase();
    $.each($(".reaction"), function() {
        _rxnRow = this;
        _matched = false;
        _rxn_spc_list = [""];
        _reactants = _rxnRow.getElementsByClassName('reactants');
        $.each($(_reactants).find("a"), function() {
            _a = this;
            $.each($(_a).find("img"), function() {
                _spec = this;
                _rxn_spc_list.push(_spec.getAttribute("alt").toLowerCase());
          
            });
        });
        _products = _rxnRow.getElementsByClassName('products');
        $.each($(_products).find("a"), function() {
            _a = this;
            $.each($(_a).find("img"), function() {
                _spec = this;
                _rxn_spc_list.push(_spec.getAttribute("alt").toLowerCase());
            
            });
        });

        if(_rxn_spc_list.indexOf(_r1) != -1 && _rxn_spc_list.indexOf(_r2) != -1 &&_rxn_spc_list.indexOf(_p1) != -1 && _rxn_spc_list.indexOf(_p2) != -1){
            _matched = true
        }
        if(_matched == true)
        $(_rxnRow).removeClass("hide_reaction");
        else
        $(_rxnRow).addClass("hide_reaction");

     });
}

$(document).ready(function() {
    uncheckThermoDetails();
    checkAllFamilies();
    uncheckAllDetails();
});
</script>
</head>

<body>

<h1>{{ title }}</h1>

<h2>Species ({{ species|length }})</h2>

<form id='thermoSelector' action="">
<input type="checkbox" id="thermoComment" name="detail" value="thermoComment" onclick="updateThermoDetails(this);" checked="false"><label for="thermoComment"><b>Show Thermo Details</b></label><br>
</form>


<table class="speciesList" hide_thermoComment>
<tr><th>Index</th><th>Thermo<br> H298 (kcal/mol), S298 (cal/mol*K), Cp (cal/mol*K)</th><th>Structure</th><th>Label</th><th>SMILES</th><th>MW<br> (g/mol)</th></tr>
{% for spec in species %}

<tr class="species">
    <td class="index" valign="top">
    {{ spec.index }}.</td>
    
    
 <td class="thermo" valign="top">
 
{% if spec.thermo %}
        <table class="thermo" align="left">
            <tr>
                <th>H298</th>
                <th>S298</th>
                <th>Cp300</th>
                <th>Cp500</th>
                <th>Cp1000</th>
                <th>Cp1500</th>
            </tr>
            <tr>
                <td>
                {% if spec.thermo.Tmin.value_si <= 298 %}                    
                {{ "%.2f"|format(spec.thermo.get_enthalpy(298) / 4184) }}
                {% endif %} </td>
                <td>{% if spec.thermo.Tmin.value_si <= 298 %}
                {{ "%.2f"|format(spec.thermo.get_entropy(298) / 4.184) }}
                {% endif %}</td>
                <td>{{ "%.2f"|format(spec.thermo.get_heat_capacity(300) / 4.184) }}</td>
                <td>{{ "%.2f"|format(spec.thermo.get_heat_capacity(500) / 4.184) }}</td>
                <td>{{ "%.2f"|format(spec.thermo.get_heat_capacity(1000) / 4.184) }}</td>
                <td>{{ "%.2f"|format(spec.thermo.get_heat_capacity(1500) / 4.184) }}</td>
            </tr>
<tr><td colspan="6" class="thermoComment">
<div id="thermoComment" class="thermoComment">{{textwrap.fill(spec.thermo.comment,80).replace('\n','<br>')}}</div>
</td></tr>
        </table>
    
  {% endif %}

 </td>
    
    <td class="structure" valign="top"><a href={{ spec.molecule[0].get_url() }}><img src="species/{{ spec|replace('#','%23') }}.png" alt="{{ get_species_identifier(spec) }}" title="{{ get_species_identifier(spec) }}"></a></td>
    <td class="label" valign="top">{{ get_species_identifier(spec) }}</td>
    <td class="SMILES" valign="top">{{ spec.molecule[0].to_smiles() }}</td>
    
  <td class="MW" valign="top">{{ "%.2f"|format(spec.molecule[0].get_molecular_weight() * 1000) }}</td>
    
</tr>
{% endfor %}
</table>

<h2>Reactions ({{ reactions|length }})</h2>

<form id='familySelector' action="">
<h4>Reaction families:</h4>
{% for family in families %}    <input type="checkbox" id="{{ family|csssafe }}" name="family" value="{{ family|csssafe }}" checked="checked" onclick="updateFamily(this);"><label for="{{ family|csssafe }}">{{ family }} ({{ family_count[family] }} rxn{{ 's' if family_count[family] != 1 }})</label><br>
{% endfor %}
<a href="javascript:checkAllFamilies();" onclick="checkAllFamilies()">check all</a> &nbsp; &nbsp; <a href="javascript:uncheckAllFamilies();" onclick="uncheckAllFamilies();">uncheck all</a><br>

<h4>Reaction Details:</h4>
<input type="checkbox" id="kinetics" name="detail" value="kinetics" onclick="updateDetails(this);"><label for="kinetics">Kinetics</label><br>
<input type="checkbox" id="energy" name="detail" value="energy" onclick="updateDetails(this);"><label for="energy">Heats of Reaction</label><br>
<input type="checkbox" id="chemkin" name="detail" value="chemkin" onclick="updateDetails(this);"><label for="chemkin">Chemkin strings</label><br>
<a href="javascript:checkAllDetails();" onclick="checkAllDetails()">check all</a> &nbsp; &nbsp; <a href="javascript:uncheckAllDetails();" onclick="uncheckAllDetails();">uncheck all</a>
</form>

<h4>Reaction Filter:</h4>

<form id="reactionFilter">
  Reactant 1: <input type="text" id="reactant1" value=""> &nbsp;
  Reactant 2: <input type="text" id="reactant2" value=""> &nbsp;
  Product 1: <input type="text" id="product1" value=""> &nbsp;
  Product 2: <input type="text" id="product2" value=""><br><br>

  
</form>

<input type="button" onclick="submitReactionFilter()" value="Search"> &nbsp;
<input type="button" onclick="resetReactionFilter()" value="Clear Filter"> 

<h4>Reaction List:</h4>

<table class="reactionList">
<thead>
<tr><th>Index</th><th colspan="3" style="text-align: center;">Reaction</th><th>Family</th></tr>
</thead>
{% for rxn in reactions %}
<tbody class="reaction">
<tr class="{{ rxn.get_source()|csssafe }} rxnStart">
    <td class="index"><a href="{{ rxn.get_url() }}" title="Search on RMG website" class="searchlink">{{ rxn.index }}.</a></td>
    <td class="reactants">{% for reactant in rxn.reactants %}<a href="{{ reactant.molecule[0].get_url() }}"><img src="species/{{ reactant|replace('#','%23') }}.png" alt="{{ get_species_identifier(reactant) }}" title="{{ get_species_identifier(reactant) }}, MW = {{ "%.2f g/mol"|format(reactant.molecule[0].get_molecular_weight() * 1000) }}" {% if reactant.contains_surface_site() %}class="surface_species" {% endif %}></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
    <td class="reactionArrow">{% if rxn.reversible %}&hArr;{% else %}&rarr;{% endif %}</td>
    <td class="products">{% for product in rxn.products %}<a href="{{ product.molecule[0].get_url() }}"><img src="species/{{ product|replace('#','%23') }}.png" alt="{{ get_species_identifier(product) }}" title="{{ get_species_identifier(product) }}, MW = {{ "%.2f g/mol"|format(product.molecule[0].get_molecular_weight() * 1000) }}" {% if product.contains_surface_site() %}class="surface_species" {% endif %}></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
    <td class="family">{{ rxn.get_source() }}</td>
</tr>
<tr class="kinetics {{ rxn.get_source()|csssafe }} hide_kinetics">
    <td></td>
    <td colspan="4">{{ rxn.kinetics.to_html() }}</td>
</tr>
<tr class="energy {{ rxn.get_source()|csssafe }} hide_energy">
    <td></td>
    <td colspan="3"><b>H298 (kcal/mol)</b> = {{ '%0.2f'| format(rxn.get_enthalpy_of_reaction(298)/4184) }}
    <br><b>S298 (cal/mol*K)</b> = {{ '%0.2f'| format(rxn.get_entropy_of_reaction(298)/4.184) }}
    <br><b>G298 (kcal/mol)</b> = {{ '%0.2f'| format(rxn.get_free_energy_of_reaction(298)/4184) }}</td>
    <td></td>
</tr>
<tr class="chemkin {{ rxn.get_source()|csssafe }} hide_chemkin">
    <td></td>
    <td colspan="4">{{ rxn.to_chemkin(species) }}</td>
</tr>
</tbody>
{% endfor %}

</table>

</body>

</html>
""")

    f = open(path, 'w')
    f.write(template.render(title=title, species=species, reactions=reactions, families=families,
                            family_count=family_count, get_species_identifier=get_species_identifier,
                            textwrap=textwrap))
    f.close()


def save_diff_html(path, common_species_list, species_list1, species_list2, common_reactions, unique_reactions1,
                   unique_reactions2):
    """
    This function outputs the species and reactions on an HTML page
    for the comparison of two RMG models.
    """
    from rmgpy.rmg.model import PDepReaction
    from rmgpy.kinetics import MultiArrhenius, MultiPDepArrhenius

    from rmgpy.molecule.draw import MoleculeDrawer
    try:
        import jinja2
    except ImportError:
        logging.warning("jinja package not found; HTML output will not be saved.")
        return

    path = os.path.abspath(path)
    dirname = os.path.dirname(path)

    # Prepare parameters to pass to jinja template
    title = 'RMG Model Comparison'

    species_list = [spec1 for spec1, spec2 in common_species_list] + \
                   [spec2 for spec1, spec2 in common_species_list] + species_list1 + species_list2
    re_index = re.compile(r'\((\d+)\)$')

    if not os.path.isdir(os.path.join(dirname, 'species1')):
        os.makedirs(os.path.join(dirname, 'species1'))

    if not os.path.isdir(os.path.join(dirname, 'species2')):
        os.makedirs(os.path.join(dirname, 'species2'))

    for spec1, spec2 in common_species_list:
        # if the species dictionary came from an RMG-Java job, make them prettier
        # We use the presence of a trailing index on the label to discern this
        # (A single open parenthesis is not enough (e.g. when using SMILES strings as labels!)
        match1 = re_index.search(spec1.label)
        if match1:
            spec1.index = int(match1.group(0)[1:-1])
            spec1.label = spec1.label[0:match1.start()]

        match2 = re_index.search(spec2.label)
        if match2:
            spec2.index = int(match2.group(0)[1:-1])
            spec2.label = spec2.label[0:match2.start()]

            # Draw molecules if necessary
        fstr = os.path.join(dirname, 'species1', '{0}.png'.format(spec1))
        if not os.path.exists(fstr):
            try:
                MoleculeDrawer().draw(spec1.molecule[0], 'png', fstr)
            except IndexError:
                raise OutputError('{0} species could not be drawn because it did not contain a molecular structure. '
                                  'Please recheck your files.'.format(get_species_identifier(spec1)))

        fstr = os.path.join(dirname, 'species2', '{0}.png'.format(spec2))
        if not os.path.exists(fstr):
            try:
                MoleculeDrawer().draw(spec2.molecule[0], 'png', fstr)
            except IndexError:
                raise OutputError(
                    '{0} species could not be drawn because it did not contain a molecular structure. Please recheck '
                    'your files.'.format(get_species_identifier(spec2)))

    for spec in species_list1:
        match = re_index.search(spec.label)
        if match:
            spec.index = int(match.group(0)[1:-1])
            spec.label = spec.label[0:match.start()]
        # Draw molecules if necessary
        fstr = os.path.join(dirname, 'species1', '{0}.png'.format(spec))
        if not os.path.exists(fstr):
            try:
                MoleculeDrawer().draw(spec.molecule[0], 'png', fstr)
            except IndexError:
                raise OutputError('{0} species could not be drawn because it did not contain a molecular structure. '
                                  'Please recheck your files.'.format(get_species_identifier(spec)))

    for spec in species_list2:
        match = re_index.search(spec.label)
        if match:
            spec.index = int(match.group(0)[1:-1])
            spec.label = spec.label[0:match.start()]
        # Draw molecules if necessary
        fstr = os.path.join(dirname, 'species2', '{0}.png'.format(spec))
        if not os.path.exists(fstr):
            try:
                MoleculeDrawer().draw(spec.molecule[0], 'png', fstr)
            except IndexError:
                raise OutputError('{0} species could not be drawn because it did not contain a molecular structure. '
                                  'Please recheck your files.'.format(get_species_identifier(spec)))

    # Add pictures for species that may not have different thermo but are in reactions with different kinetics
    all_rxns = [rxnTuple[0] for rxnTuple in common_reactions] + unique_reactions1 + unique_reactions2
    all_species = []
    for rxn in all_rxns:
        for prod in rxn.products:
            all_species.append(prod)
        for rxt in rxn.reactants:
            all_species.append(rxt)
    all_species = set(all_species)

    for spec in all_species:
        match = re_index.search(spec.label)
        if match:
            spec.index = int(match.group(0)[1:-1])
            spec.label = spec.label[0:match.start()]
        # Draw molecules if necessary
        fstr = os.path.join(dirname, 'species2', '{0}.png'.format(spec))
        if not os.path.exists(fstr):
            try:
                MoleculeDrawer().draw(spec.molecule[0], 'png', fstr)
            except IndexError:
                raise OutputError('{0} species could not be drawn because it did not contain a molecular structure. '
                                  'Please recheck your files.'.format(get_species_identifier(spec)))

    family_count1 = {}
    family_count2 = {}
    for rxn1, rxn2 in common_reactions:
        if isinstance(rxn2.kinetics, (MultiArrhenius, MultiPDepArrhenius)):
            rxn2.duplicate = True
        if isinstance(rxn1, PDepReaction):
            family = "PDepNetwork"
        else:
            family = rxn1.get_source()
        if family in family_count1:
            family_count1[family] += 1
            family_count2[family] += 1
        else:
            family_count1[family] = 1
            family_count2[family] = 1

    for rxn in unique_reactions1:
        if isinstance(rxn, PDepReaction):
            family = "PDepNetwork"
        else:
            family = rxn.get_source()
        if family in family_count1:
            family_count1[family] += 1
        else:
            family_count1[family] = 1

    for rxn in unique_reactions2:
        if isinstance(rxn, PDepReaction):
            family = "PDepNetwork"
        else:
            family = rxn.get_source()
        if family in family_count2:
            family_count2[family] += 1
        else:
            family_count2[family] = 1

    families1 = list(family_count1.keys())
    families2 = list(family_count2.keys())
    families1.sort()
    families2.sort()

    # jinja2 filters etc.
    to_remove_from_css_names = re.compile(r'[/.\-+,]')

    def csssafe(input):
        """Replace unsafe CSS class name characters with an underscore."""
        return to_remove_from_css_names.sub('_', input)

    environment = jinja2.Environment()
    environment.filters['csssafe'] = csssafe

    # Make HTML file
    template = environment.from_string(
"""<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN">
<html lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html;charset=utf-8" >
    <title>{{ title }}</title>
    <style type="text/css">
     body {
        font-family: sans-serif;
    }
    a {
        color: #993333;
        text-decoration: none;
    }
    a:visited {
        color: #993333;
    }
    a:hover {
        text-decoration: underline;
    }

    
    table.speciesList, table.reactionList {
        border-collapse: collapse;
        align: center;
    }

        
    table.speciesList th, table.reactionList th {
        text-align: left;
        vertical-align: top;
    }
    
    table.speciesList th, {
        font-size: small;
    }
    
    tr.reaction {
        border-top: 1px solid #808080;        
        padding: 10px;
    }
    td.reactants {
        text-align: right;
    }
    td.products {
        text-align: left;
    }
    td.reactionArrow {
        text-align: center;
        font-size: 16px;
    }
    td.species img, td.reactants img, td.products img {
        vertical-align: middle;
    }
    
    tr.species{
        border-top: 1px solid #808080;        
    }
    tr.commonSpecies{
        border-top: 1px solid #808080;        
    }
    
    td, .speciesList th{
        padding: 10px;
        }
    
    td.index {
    width: 50px;
    }
    
    .thermo td, .thermo th{
    padding: 2px;
    }
    
    .kinetics td,  .kinetics th {
    padding: 2px;
    }
    
    tr.kinetics {
        font-size: small;
    }
    .KineticsData {
        # border: 1px solid gray;
    }
    .KineticsData th {
        width: 15em;
        word-wrap: none;
    }
    .KineticsData td {
        width: 3em;
    }
    
    .energy {
        font-size: small;
        }
    .chemkin, .KineticsData_repr {
       white-space: pre-wrap;
       font-size: x-small;
       font-family: "Andale Mono", monospace;
    }
    .thermoComment {
       white-space: pre-wrap;
       font-size: small;
       font-family: "Andale Mono", monospace;
    }
    .hide_kinetics .kinetics{
       display: none !important;
    }
    .hide_chemkin .chemkin{
       display: none !important;
    }   
    .hide_energy .energy{
       display: none !important;
    }
    .hide_thermoComment .thermoComment{
       display: none !important;
    }

           
</style>
<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.4.1/jquery.min.js"></script>
<script type="text/javascript" src="../../../external/jquery.min.js"></script>
<script type="text/javascript">
function updateFamily(family) {
    if (family.checked) {
        $("."+family.value).show();
    } else {
        $("."+family.value).hide();
    }
}
function updateDetails(type) {
    if (type.checked) {
        $(".reactionList").removeClass("hide_"+type.value);
    } else {
        $(".reactionList").addClass("hide_"+type.value);
    }
}
function checkAllFamilies() {
    $("#familySelector").find("[name='family']").each(function() { this.checked = true; updateFamily(this); });
    return false;
}
function uncheckAllFamilies() {
    $("#familySelector").find("[name='family']").each(function() { this.checked = false; updateFamily(this); });
    return false;
}
function checkAllDetails() {
    $("#familySelector").find("[name='detail']").each(function() { this.checked = true; updateDetails(this); });
    return false;
}
function uncheckAllDetails() {
    $("#familySelector").find("[name='detail']").each(function() { this.checked = false; updateDetails(this); });
    return false;
}



function updateThermoDetails(type) {
    if (type.checked) {
        $(".thermoComment").removeClass("hide_"+type.value);
    } else {
        $(".thermoComment").addClass("hide_"+type.value);
    }
}

function uncheckThermoDetails() {
    $("#thermoSelector").find("[name='detail']").each(function() { this.checked = false; updateThermoDetails(this); });
    return false;
}

$(document).ready(function() {
    uncheckThermoDetails();
    checkAllFamilies();
    uncheckAllDetails();
});
    </script>
</head>

<body>

<h1>{{ title }}</h1>

<div align = "center">
<form id='thermoSelector' action="">
<input type="checkbox" id="thermoComment" name="detail" value="thermoComment" onclick="updateThermoDetails(this);" checked="false"><label for="thermoComment"><b>Show Thermo Details</b></label><br>
</form></div>

<h2 align="center">Common Species ({{ common_species|length }})</h2>

<table class="speciesList" align="center" hide_thermoComment>
    <tr><td align="center"><h3>Model 1</h3></td><td align="center"><h3>Model 2</h3></td></tr>
    
    {% for spec1, spec2 in common_species %}
    <tr class="commonSpecies">
        <td width="100%" colspan="2">
            <table align="center">
                <tr><th>Structure</th><th>SMILES</th><th>MW (g/mol)</th></tr>
                <tr>
                    <td>{{ spec1.molecule[0].to_smiles() }}</td>
                    <td class="structure" align="center"><a href="{{spec1.molecule[0].get_url()}}"><img src="species1/{{ spec1|replace('#','%23') }}.png"></a></td>
                    <td>{{ "%.2f"|format(spec1.molecule[0].get_molecular_weight() * 1000) }}</td>
                </tr>
            </table>
        </td>
    </tr>
    <tr>
    <td width="50%">
    <table width="80%">
    <tr><td width="20%" valign="top">{{ spec1.index }}. </td>
        <td width="20%"  valign="top">{{get_species_identifier(spec1)}}</td>
        <td width="80%"  valign="top">
        {% if spec1.thermo %}
            <table width="80%"  class="thermo" valign="top"> 
                <tr>
                    <th>H298</th>
                    <th>S298</th>
                    <th>Cp300</th>
                    <th>Cp500</th>
                    <th>Cp1000</th>
                    <th>Cp1500</th>
                </tr>
                <tr>
                    <td>{% if spec1.thermo.Tmin.value_si <= 298 %}
                    {{ "%.2f"|format(spec1.thermo.get_enthalpy(298) / 4184) }}
                    {% endif %}</td>
                    <td>
                    {% if spec1.thermo.Tmin.value_si <= 298 %}
                    {{ "%.2f"|format(spec1.thermo.get_entropy(298) / 4.184) }}
                    {% endif %}</td>
                    <td>{{ "%.2f"|format(spec1.thermo.get_heat_capacity(300) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec1.thermo.get_heat_capacity(500) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec1.thermo.get_heat_capacity(1000) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec1.thermo.get_heat_capacity(1500) / 4.184) }}</td>
                </tr>
                <tr><td colspan="6" class="thermoComment">
<div id="thermoComment" class="thermoComment">{{textwrap.fill(spec1.thermo.comment,80).replace('\n','<br>')}}</div>
</td></tr>
            </table>
            
            {% endif %}
        </td></tr>
        </table>
        </td>
        
        <td width="50%">
        <table width="80%" class="thermo" valign="top">
    <tr><td width="20%"  valign="top">{{ spec2.index }}. </td>
        <td width="20%"  valign="top">{{get_species_identifier(spec2)}}</td>
        <td width="80%"  valign="top">
        
        {% if spec2.thermo %}
            <table width="100%">
                <tr>
                    <th>H298</th>
                    <th>S298</th>
                    <th>Cp300</th>
                    <th>Cp500</th>
                    <th>Cp1000</th>
                    <th>Cp1500</th>
                </tr>
                <tr>
                    <td>{% if spec2.thermo.Tmin.value_si <= 298 %}
                    {{ "%.2f"|format(spec2.thermo.get_enthalpy(298) / 4184) }}
                    {% endif %}</td>
                    <td>{% if spec2.thermo.Tmin.value_si <= 298 %}
                    {{ "%.2f"|format(spec2.thermo.get_entropy(298) / 4.184) }}
                    {% endif %}</td>
                    <td>{{ "%.2f"|format(spec2.thermo.get_heat_capacity(300) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec2.thermo.get_heat_capacity(500) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec2.thermo.get_heat_capacity(1000) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec2.thermo.get_heat_capacity(1500) / 4.184) }}</td>
                </tr>
                                <tr><td colspan="6" class="thermoComment">
<div id="thermoComment" class="thermoComment">{{textwrap.fill(spec2.thermo.comment,80).replace('\n','<br>')}}</div>
</td></tr>
            </table>
            {% endif %}
        </td>
    </tr></table>
    </td></tr>
    
    {% if spec1.thermo and spec2.thermo %}
    {% if spec1.thermo.is_identical_to(spec2.thermo) %}
    <tr width=100%>
         <td colspan="2" valign="top" width=50%><div align="center"><font color="blue">IDENTICAL THERMO WAS FOUND FOR THIS SPECIES.</font></div>
    </tr>
    {% elif spec1.thermo.is_similar_to(spec2.thermo) %}
    <tr width=100%>
         <td colspan="2" valign="top" width=50%><div align="center"><font color="green">SIMILAR THERMO WAS FOUND FOR THIS SPECIES.</font></div>
    </tr>
    {% else %}
     <tr width=100%>
         <td colspan="2" valign="top" width=50%><div align="center"><font color="red">DIFFERENT THERMO WAS FOUND FOR THIS SPECIES.</font></div>
    </tr>
    {% endif%}{% endif %}
    
    
    {% endfor %}
</table>

<table width=100%>
<tr colspan="2">
<td width=50% valign="top">

<h2>Model 1: Unique Species ({{ species_list1|length }})</h2>
<table class="speciesList" width="80%" hide_thermoComment>
    <tr><th>Index</th><th>Structure</th><th>Label</th><th>SMILES</th><th>MW (g/mol)</th></tr>
    {% for spec in species_list1 %}
    <tr class="species">
        <td class="index">
        {{ spec.index }}.</td>
        <td class="structure"><a href="{{ spec.molecule[0].get_url() }}"><img src="species1/{{ spec|replace('#','%23') }}.png" alt="{{ get_species_identifier(spec) }}" title="{{ get_species_identifier(spec) }}"></a></td>
        <td class="label">{{ get_species_identifier(spec) }}</td>
        <td>{{spec.molecule[0].to_smiles()}}</td>
        <td>{{ "%.2f"|format(spec.molecule[0].get_molecular_weight() * 1000) }}</td>
    </tr>
    <tr><td colspan="5">
    
    {% if spec.thermo %}
            <table width="80%"  class="thermo" valign="top">
                <tr>
                    <th>H298</th>
                    <th>S298</th>
                    <th>Cp300</th>
                    <th>Cp500</th>
                    <th>Cp1000</th>
                    <th>Cp1500</th>
                </tr>
                <tr>
                    <td>{% if spec.thermo.Tmin.value_si <= 298 %}
                    {{ "%.2f"|format(spec.thermo.get_enthalpy(298) / 4184) }}
                    {% endif %}</td>
                    <td>
                    {% if spec.thermo.Tmin.value_si <= 298 %}
                    {{ "%.2f"|format(spec.thermo.get_entropy(298) / 4.184) }}
                    {% endif %}</td>
                    <td>{{ "%.2f"|format(spec.thermo.get_heat_capacity(300) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec.thermo.get_heat_capacity(500) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec.thermo.get_heat_capacity(1000) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec.thermo.get_heat_capacity(1500) / 4.184) }}</td>
                </tr>
                                <tr><td colspan="6" class="thermoComment">
<div id="thermoComment" class="thermoComment">{{textwrap.fill(spec.thermo.comment,80).replace('\n','<br>')}}</div>
</td></tr>
            </table>
            {% endif %}
    
    </td></tr>
    
    {% endfor %}
</table>
</td>
<td width=50% valign="top">
<h2>Model 2: Unique Species ({{ species_list2|length }})</h2>
<table class="speciesList" width="80%" hide_thermoComment>
    <tr><th>Index</th><th>Structure</th><th>Label</th><th>SMILES</th><th>MW (g/mol)</th></tr>
    {% for spec in species_list2 %}
    <tr class="species">
        <td class="index">
        {{ spec.index }}.</td>
        <td class="structure"><a href="{{ spec.molecule[0].get_url() }}"><img src="species2/{{ spec|replace('#','%23') }}.png" alt="{{ get_species_identifier(spec) }}" title="{{ get_species_identifier(spec) }}"></a></td>
        <td class="label">{{ get_species_identifier(spec) }}</td>
        <td>{{spec.molecule[0].to_smiles()}}</td>
        <td>{{ "%.2f"|format(spec.molecule[0].get_molecular_weight() * 1000) }}</td>
    </tr>
    
        <tr><td colspan="5">
    
    {% if spec.thermo %}
            <table width="80%"  class="thermo" valign="top">
                <tr>
                    <th>H298</th>
                    <th>S298</th>
                    <th>Cp300</th>
                    <th>Cp500</th>
                    <th>Cp1000</th>
                    <th>Cp1500</th>
                </tr>
                <tr>
                    <td>{% if spec.thermo.Tmin.value_si <= 298 %}
                    {{ "%.2f"|format(spec.thermo.get_enthalpy(298) / 4184) }}
                    {% endif %}</td>
                    <td>
                    {% if spec.thermo.Tmin.value_si <= 298 %}
                    {{ "%.2f"|format(spec.thermo.get_entropy(298) / 4.184) }}
                    {% endif %}</td>
                    <td>{{ "%.2f"|format(spec.thermo.get_heat_capacity(300) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec.thermo.get_heat_capacity(500) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec.thermo.get_heat_capacity(1000) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec.thermo.get_heat_capacity(1500) / 4.184) }}</td>
                </tr>
                                <tr><td colspan="6" class="thermoComment">
<div id="thermoComment" class="thermoComment">{{textwrap.fill(spec.thermo.comment,80).replace('\n','<br>')}}</div>
</td></tr>
            </table>
            {% endif %}
    
    </td></tr>
    
    {% endfor %}
</table>
</td></tr>



<tr><td colspan="2" align="center">
<form id='familySelector' action="">
    <h4>Reaction families:</h4>
{% for family in families_union %}    <input type="checkbox" id="{{ family|csssafe }}" name="family" value="{{ family|csssafe }}" checked="checked" onclick="updateFamily(this);"><label for="{{ family|csssafe }}">{{ family }}</label><br>
{% endfor %}
    <a href="javascript:checkAllFamilies();" onclick="checkAllFamilies()">check all</a> &nbsp; &nbsp; <a href="javascript:uncheckAllFamilies();" onclick="uncheckAllFamilies();">uncheck all</a><br>

    <h4>Reaction Details:</h4>
    <input type="checkbox" id="kinetics" name="detail" value="kinetics" onclick="updateDetails(this);"><label for="kinetics">Kinetics</label><br>
    <input type="checkbox" id="energy" name="detail" value="energy" onclick="updateDetails(this);"><label for="energy">Heats of Reaction</label><br>
    <input type="checkbox" id="chemkin" name="detail" value="chemkin" onclick="updateDetails(this);"><label for="chemkin">Chemkin strings</label><br>
    <a href="javascript:checkAllDetails();" onclick="checkAllDetails()">check all</a> &nbsp; &nbsp; <a href="javascript:uncheckAllDetails();" onclick="uncheckAllDetails();">uncheck all</a>
</form>



</td></tr>


<tr colspan="2">
<td width=50% valign="top">
<h2>Model 1 Reactions ({{ common_reactions|length + unique_reactions1|length}})</h2>
</td>

<td width=50% valign="top">
<h2>Model 2 Reactions ({{ common_reactions|length + unique_reactions2|length}})</h2>
</td>
</tr>


<tr colspan="2"><td width=100% align="center" colspan="2">
<h2>Common Reactions ({{ common_reactions|length}})</h2></td></tr>


<tr colspan="2"><td width=100% colspan="2">

<table class="reactionList" hide_kinetics hide_chemkin cellpadding="10" align="center">
    <tr colspan="4" width=100%><th>Index.</th><th>Family</th><th>Index.</th><th>Family</th></tr>

    {% for rxn1, rxn2 in common_reactions %}


<tr class="reaction  {{ rxn1.get_source()|csssafe }}">

<td width=100% colspan="4" align="center">


<table width=100%>
<tr>
<td width=100% colspan="4">
<table align="center">
<tr>
    <td class="reactants" align="right">{% for reactant in rxn1.reactants %}<a href="{{reactant.molecule[0].get_url() }}"><img src="species1/{{ reactant|replace('#','%23') }}.png" alt="{{ reactant }}" title="{{ reactant }}, MW = {{ "%.2f"|format(reactant.molecule[0].get_molecular_weight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
    <td class="reactionArrow" align="center">{% if rxn1.reversible %}&hArr;{% else %}&rarr;{% endif %}</td>
    <td class="products" align="left">{% for product in rxn1.products %}<a href="{{product.molecule[0].get_url()}}"><img src="species1/{{ product|replace('#','%23') }}.png" alt="{{ product }}" title="{{ product }}, MW = {{ "%.2f"|format(product.molecule[0].get_molecular_weight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
</tr>
</table>
</td>
</tr>


{% if rxn1.kinetics.is_identical_to(rxn2.kinetics) %}

 <tr width=100%>
     <td colspan="4" valign="top" width=50%><div align="center"><font color="blue">IDENTICAL KINETICS WERE FOUND FOR THIS REACTION.</font></div>

</tr>
{% elif rxn1.kinetics.is_similar_to(rxn2.kinetics) %}

 <tr width=100%>
     <td colspan="4" valign="top" width=50%><div align="center"><font color="green">SIMILAR KINETICS WERE FOUND FOR THIS REACTION.</font></div>

</tr>

{% else %}
     <tr width=100%>
         <td colspan="4" valign="top" width=50%><div align="center"><font color="red">DIFFERENT KINETICS WERE FOUND FOR THIS REACTION.</font></div>
    </tr>
{% endif%}


<tr width=100%>
     <td class="index" width=10%><a href="{{ rxn1.get_url() }}" title="Search on RMG website" class="searchlink">{{ rxn1.index }}.</a></td>
     <td class="family" width=40%>{{ rxn1.get_source() }}</td>

     <td class="index" width=10%><a href="{{ rxn2.get_url() }}" title="Search on RMG website" class="searchlink">{{ rxn2.index }}.</a></td>
     <td class="family" width=40%>{{ rxn2.get_source() }}</td>
 </tr>

<tr "width=100%" class="kinetics">{% if not rxn1.is_isomorphic(rxn2, either_direction=False) %} 
<td colspan="2" width=50%></td>
<td colspan="2" width=50%>* Reaction was found in reverse 

{% if not rxn2.duplicate %}

<P><b>Fitted Reverse Kinetics:</b>
{% if not rxn2.kinetics.is_pressure_dependent() %}
{{rxn2.generate_reverse_rate_coefficient(surface_site_density=2.483e-05).to_html() }}
{% else %} Pressure dependent
{% endif %}
{% endif %}

<P><b>Original Kinetics:</b>

{% endif %}</td>
</tr>

<tr width=100% class="kinetics">
     <td colspan="2" valign="top" width=50%>
     
     {{ rxn1.kinetics.to_html() }}</td>
     <td colspan="2" valign="top" width=50%>
     {{ rxn2.kinetics.to_html() }}</td>
</tr>

<tr width=100% class="energy">

    <td colspan="2" valign="top" width=50%>
    <b>H298 (kcal/mol)</b> = {{ '%0.2f'| format(rxn1.get_enthalpy_of_reaction(298)/4184) }}
    <br><b>S298 (cal/mol*K)</b> = {{ '%0.2f'| format(rxn1.get_entropy_of_reaction(298)/4.184) }}
    <br><b>G298 (kcal/mol)</b> = {{ '%0.2f'| format(rxn1.get_free_energy_of_reaction(298)/4184) }}
    </td>
    
    <td colspan="2" valign="top" width=50%><b>H298 (kcal/mol)</b> = {{ '%0.2f'| format(rxn2.get_enthalpy_of_reaction(298)/4184) }}
    <br><b>S298 (cal/mol*K)</b> = {{ '%0.2f'| format(rxn2.get_entropy_of_reaction(298)/4.184) }}
    <br><b>G298 (kcal/mol)</b> = {{ '%0.2f'| format(rxn2.get_free_energy_of_reaction(298)/4184) }}
    </td>

</tr>

<tr width=100% class="chemkin">
    <td colspan="2" valign="top" width=50%><font size="1pt" face="courier">{{ rxn1.to_chemkin(species_list) }}</font></td>
    <td colspan="2" valign="top" width=50%><font size="1pt" face="courier">{{ rxn2.to_chemkin(species_list) }}</font></td>
</tr>


</td></tr></table>
</td></tr>
{% endfor %}

</table>


<tr>
<td width=50% valign="top">
<h2>Model 1: Unique Reactions ({{ unique_reactions1|length}})</h2>
<br>
<table class="reactionList" hide_kinetics hide_chemkin >
    <tr><th>Index</th><th colspan="3" style="text-align: center;">Reaction</th><th>Family</th></tr>
    {% for rxn in unique_reactions1 %}
    <tr class="reaction {{ rxn.get_source()|csssafe }}">
        <td class="index"><a href="{{ rxn.get_url() }}" title="Search on RMG website" class="searchlink">{{ rxn.index }}.</a></td>
        <td class="reactants">{% for reactant in rxn.reactants %}<a href="{{ reactant.molecule[0].get_url() }}"><img src="species1/{{ reactant|replace('#','%23') }}.png" alt="{{ reactant }}" title="{{ reactant }}, MW = {{ "%.2f"|format(reactant.molecule[0].get_molecular_weight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="reactionArrow">{% if rxn.reversible %}&hArr;{% else %}&rarr;{% endif %}</td>
        <td class="products">{% for product in rxn.products %}<a href="{{ product.molecule[0].get_url() }}"><img src="species1/{{ product|replace('#','%23') }}.png" alt="{{ product }}" title="{{ product }}, MW = {{ "%.2f"|format(product.molecule[0].get_molecular_weight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="family">{{ rxn.get_source() }}</td>
    </tr>
    <tr class="kinetics {{ rxn.get_source()|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.kinetics.to_html() }}</td>
    </tr>
    <tr class="energy {{ rxn.get_source()|csssafe }} hide_energy">
    <td></td>
    <td colspan="3"><b>H298 (kcal/mol)</b> = {{ '%0.2f'| format(rxn.get_enthalpy_of_reaction(298)/4184) }}
    <br><b>S298 (cal/mol*K)</b> = {{ '%0.2f'| format(rxn.get_entropy_of_reaction(298)/4.184) }}
    <br><b>G298 (kcal/mol)</b> = {{ '%0.2f'| format(rxn.get_free_energy_of_reaction(298)/4184) }}</td>
    <td></td>
</tr>

    <tr class="chemkin {{ rxn.get_source()|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.to_chemkin(species) }}</td>
    </tr>
    {% endfor %}
    </table>

</td>

<td width=50% valign="top">
<h2>Model 2: Unique Reactions ({{ unique_reactions2|length}})</h2>
<br>
<table class="reactionList" hide_kinetics hide_chemkin>
    <tr><th>Index</th><th colspan="3" style="text-align: center;">Reaction</th><th>Family</th></tr>
    {% for rxn in unique_reactions2 %}
    <tr class="reaction {{ rxn.get_source()|csssafe }}">
        <td class="index"><a href="{{ rxn.get_url() }}" title="Search on RMG website" class="searchlink">{{ rxn.index }}.</a></td>
        <td class="reactants">{% for reactant in rxn.reactants %}<a href="{{ reactant.molecule[0].get_url() }}"><img src="species2/{{ reactant|replace('#','%23') }}.png" alt="{{ reactant }}" title="{{ reactant }}, MW = {{ "%.2f"|format(reactant.molecule[0].get_molecular_weight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="reactionArrow">{% if rxn.reversible %}&hArr;{% else %}&rarr;{% endif %}</td>
        <td class="products">{% for product in rxn.products %}<a href="{{ product.molecule[0].get_url() }}"><img src="species2/{{ product|replace('#','%23') }}.png" alt="{{ product }}" title="{{ product }}, MW = {{ "%.2f"|format(product.molecule[0].get_molecular_weight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="family">{{ rxn.get_source() }}</td>
    </tr>
    <tr class="kinetics {{ rxn.get_source()|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.kinetics.to_html() }}</td>
    </tr>
    <tr class="energy {{ rxn.get_source()|csssafe }} hide_energy">
    <td></td>
    <td colspan="3"><b>H298 (kcal/mol)</b> = {{ '%0.2f'| format(rxn.get_enthalpy_of_reaction(298)/4184) }}
    <br><b>S298 (cal/mol*K)</b> = {{ '%0.2f'| format(rxn.get_entropy_of_reaction(298)/4.184) }}
    <br><b>G298 (kcal/mol)</b> = {{ '%0.2f'| format(rxn.get_free_energy_of_reaction(298)/4184) }}</td>
    <td></td>
</tr>
    <tr class="chemkin {{ rxn.get_source()|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.to_chemkin(species) }}</td>
    </tr>
    {% endfor %}
    </table>

</td>

</table>



</body>

</html>
""")
    f = open(path, 'w')
    f.write(template.render(title=title, common_species=common_species_list, species_list1=species_list1,
                            species_list2=species_list2,
                            common_reactions=common_reactions, unique_reactions1=unique_reactions1,
                            unique_reactions2=unique_reactions2,
                            families1=families1, families2=families2, family_count1=family_count1,
                            family_count2=family_count2, families_union=set(families1 + families2),
                            species_list=species_list,
                            get_species_identifier=get_species_identifier, textwrap=textwrap))
    f.close()


def save_output(rmg):
    """
    Save the current reaction model to a pretty HTML file.
    """
    logging.info('Saving current model core to HTML file...')
    save_output_html(os.path.join(rmg.output_directory, 'output.html'), rmg.reaction_model, 'core')

    if rmg.save_edge_species:
        logging.info('Saving current model edge to HTML file...')
        save_output_html(os.path.join(rmg.output_directory, 'output_edge.html'), rmg.reaction_model, 'edge')


class OutputHTMLWriter(object):
    """
    This class listens to a RMG subject
    and writes a HTML file with the current state of the RMG model,
    to the 'species' subfolder.

    A new instance of the class can be appended to a subject as follows:
    
    rmg = ...
    listener = OutputHTMLWriter(output_directory)
    rmg.attach(listener)

    Whenever the subject calls the .notify() method, the
    .update() method of the listener will be called.

    To stop listening to the subject, the class can be detached
    from its subject:

    rmg.detach(listener)

    """

    def __init__(self, output_directory=''):
        super(OutputHTMLWriter, self).__init__()
        make_output_subdirectory(output_directory, 'species')

    def update(self, rmg):
        save_output(rmg)
