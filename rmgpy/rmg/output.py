#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains functionality for saving the output of RMG jobs to output
files.
"""

import os.path
import logging
import re

################################################################################

class OutputError(Exception):
    """
    This exception is raised whenever an error occurs while saving output
    information. Pass a string describing the circumstances of the exceptional
    behavior.
    """
    pass

################################################################################

def saveOutputHTML(path, reactionModel, partCoreEdge='core'):
    """
    Save the current set of  species and reactions of `reactionModel` to
    an HTML file `path` on disk. As part of this process, drawings of all 
    species are created in the species folder (if they don't already exist)
    using the :mod:`rmgpy.molecule.draw` module. The :mod:`jinja`
    package is used to generate the HTML; if this package is not found, no
    HTML will be generated (but the program will carry on).
    """

    from model import PDepReaction
    
    from rmgpy.molecule.draw import MoleculeDrawer
    from rmgpy.chemkin import getSpeciesIdentifier

    try:
        import jinja2
    except ImportError:
        logging.warning("jinja2 package not found; HTML output will not be saved.")
        return

    path = os.path.abspath(path)
    dirname = os.path.dirname(path)

    # Prepare parameters to pass to jinja template
    title = 'RMG Output'
    
    if partCoreEdge == 'core':
        species = reactionModel.core.species[:] + reactionModel.outputSpeciesList
        if not os.path.isdir(os.path.join(dirname,'species')):
            os.makedirs(os.path.join(dirname,'species'))
    elif partCoreEdge == 'edge':
        species = reactionModel.edge.species[:] + reactionModel.outputSpeciesList
        if not os.path.isdir(os.path.join(dirname,'species_edge')):
            os.makedirs(os.path.join(dirname,'species_edge'))

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
        if partCoreEdge == 'core':
            fstr = os.path.join(dirname, 'species', '{0}.png'.format(spec))
        elif partCoreEdge == 'edge':
            fstr = os.path.join(dirname, 'species_edge', '{0}.png'.format(spec))
        if not os.path.exists(fstr):
            try:
                MoleculeDrawer().draw(spec.molecule[0], 'png', fstr)
            except IndexError:
                raise OutputError("{0} species could not be drawn because it did not contain a molecular structure. Please recheck your files.".format(getSpeciesIdentifier(spec)))
                
    # We want to keep species sorted in the original order in which they were added to the RMG core.
    # Rather than ordered by index
#    species.sort(key=lambda x: x.index)
    
    if partCoreEdge == 'core': 
        reactions = [rxn for rxn in reactionModel.core.reactions ] + reactionModel.outputReactionList
    elif partCoreEdge == 'edge':
        reactions = [rxn for rxn in reactionModel.edge.reactions ] + reactionModel.outputReactionList

    # We want to keep reactions sorted in original order in which they were added to core
    # rather than ordered by index
    #reactions.sort(key=lambda x: x.index)

    familyCount = {}
    for rxn in reactions:
        
        if isinstance(rxn, PDepReaction):
            family = "PDepNetwork"
        else:
            family = rxn.getSource().label
        if family in familyCount:
            familyCount[family] += 1
        else:
            familyCount[family] = 1
    families = familyCount.keys()
    families.sort()
    
    
    ## jinja2 filters etc.
    to_remove_from_css_names = re.compile('[/.\-+,]')
    def csssafe(input):
        "Replace unsafe CSS class name characters with an underscore."
        return to_remove_from_css_names.sub('_',input)
        
    environment = jinja2.Environment()
    environment.filters['csssafe'] = csssafe
    
    # Make HTML file
    if partCoreEdge == 'core': 
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
            width: 100%;
            border-collapse: collapse;
        }
        table.speciesList th, table.reactionList th {
            text-align: left;
        }
        tr.reaction td {
            border-top: 1px solid #808080;
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
        tr.comment {
            font-size: small;
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
        
        .hide_comment .comment{
            display: none !important;
        }
        .hide_kinetics .kinetics{
           display: none !important;
        }
        .hide_chemkin .chemkin{
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
    $(document).ready(function() {
        checkAllFamilies();
        uncheckAllDetails();
    });

    </script>
</head>

<body>

<h1>{{ title }}</h1>

<h2>Species ({{ species|length }})</h2>

<table class="speciesList">
    <tr><th>Index</th><th>Structure</th><th>Label</th><th>Mol. Wt. (g/mol)</th></tr>
    {% for spec in species %}
    <tr class="species">
        <td class="index">
        {{ spec.index }}.</td>
        
        
        
        <td class="structure"><a href={{ spec.molecule[0].getURL() }}><img src="species/{{ spec|replace('#','%23') }}.png" alt="{{ spec }}" title="{{ spec }}"></a></td>
        <td class="label">{{ spec.label }}</td>
        <td>{{ "%.2f"|format(spec.molecule[0].getMolecularWeight() * 1000) }}</td>
    </tr>
    {% if spec.thermo %}
    
    <tr>
     <td>
            <table align="center">
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
                    {{ "%.2f"|format(spec.thermo.getEnthalpy(298) / 4184) }}
                    {% endif %} </td>
                    <td>{% if spec.thermo.Tmin.value_si <= 298 %}
                    {{ "%.2f"|format(spec.thermo.getEntropy(298) / 4.184) }}
                    {% endif %}</td>
                    <td>{{ "%.2f"|format(spec.thermo.getHeatCapacity(300) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec.thermo.getHeatCapacity(500) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec.thermo.getHeatCapacity(1000) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec.thermo.getHeatCapacity(1500) / 4.184) }}</td>
                </tr>
            </table>
        </td></tr>
        
        {% endif %}
    {% endfor %}
</table>

<h2>Reactions ({{ reactions|length }})</h2>

<form id='familySelector' action="">
    <h4>Reaction families:</h4>
{% for family in families %}    <input type="checkbox" id="{{ family|csssafe }}" name="family" value="{{ family|csssafe }}" checked="checked" onclick="updateFamily(this);"><label for="{{ family|csssafe }}">{{ family }} ({{ familyCount[family] }} rxn{{ 's' if familyCount[family] != 1 }})</label><br>
{% endfor %}
    <a href="javascript:checkAllFamilies();" onclick="checkAllFamilies()">check all</a> &nbsp; &nbsp; <a href="javascript:uncheckAllFamilies();" onclick="uncheckAllFamilies();">uncheck all</a><br>

    <h4>Reaction Details:</h4>
    <input type="checkbox" id="kinetics" name="detail" value="kinetics" onclick="updateDetails(this);"><label for="kinetics">Kinetics</label><br>
    <input type="checkbox" id="comment" name="detail" value="comment" onclick="updateDetails(this);"><label for="comment">Comments</label><br>
    <input type="checkbox" id="chemkin" name="detail" value="chemkin" onclick="updateDetails(this);"><label for="chemkin">Chemkin strings</label><br>
    <a href="javascript:checkAllDetails();" onclick="checkAllDetails()">check all</a> &nbsp; &nbsp; <a href="javascript:uncheckAllDetails();" onclick="uncheckAllDetails();">uncheck all</a>
</form>

<h4>Reaction List:</h4>

<table class="reactionList hide_comment hide_kinetics hide_chemkin">
    <tr><th>Index</th><th colspan="3" style="text-align: center;">Reaction</th><th>Family</th></tr>
    {% for rxn in reactions %}
    <tr class="reaction {{ rxn.getSource().label|csssafe }}">
        <td class="index"><a href="{{ rxn.getURL() }}" title="Search on RMG website" class="searchlink">{{ rxn.index }}.</a></td>
        <td class="reactants">{% for reactant in rxn.reactants %}<a href="{{ reactant.molecule[0].getURL() }}"><img src="species/{{ reactant|replace('#','%23') }}.png" alt="{{ reactant }}" title="{{ reactant }}, MW = {{ "%.2f g/mol"|format(reactant.molecule[0].getMolecularWeight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="reactionArrow">{% if rxn.reversible %}&hArr;{% else %}&rarr;{% endif %}</td>
        <td class="products">{% for product in rxn.products %}<a href="{{ product.molecule[0].getURL() }}"><img src="species/{{ product|replace('#','%23') }}.png" alt="{{ product }}" title="{{ product }}, MW = {{ "%.2f g/mol"|format(product.molecule[0].getMolecularWeight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="family">{{ rxn.getSource().label }}</td>
    </tr>
    <tr class="kinetics {{ rxn.getSource().label|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.kinetics.toHTML() }}</td>
    </tr>
    <tr class="chemkin {{ rxn.getSource().label|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.toChemkin(species) }}</td>
    </tr>
    <tr class="comment {{ rxn.getSource().label|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.kinetics.comment }}</td>
    </tr>
    {% endfor %}

</table>

</body>

</html>
""")
    elif partCoreEdge == 'edge':
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
            width: 100%;
            border-collapse: collapse;
        }
        table.speciesList th, table.reactionList th {
            text-align: left;
        }
        tr.reaction td {
            border-top: 1px solid #808080;
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
        tr.comment {
            font-size: small;
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
        
        .hide_comment .comment{
            display: none !important;
        }
        .hide_kinetics .kinetics{
           display: none !important;
        }
        .hide_chemkin .chemkin{
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
    $(document).ready(function() {
        checkAllFamilies();
        uncheckAllDetails();
    });

    </script>
</head>

<body>

<h1>{{ title }}</h1>

<h2>Species ({{ species|length }})</h2>

<table class="speciesList">
    <tr><th>Index</th><th>Structure</th><th>Label</th><th>Mol. Wt. (g/mol)</th></tr>
    {% for spec in species %}
    <tr class="species">
        <td class="index">
        {{ spec.index }}.</td>
        
        
        
        <td class="structure"><a href={{ spec.molecule[0].getURL() }}><img src="species_edge/{{ spec|replace('#','%23') }}.png" alt="{{ spec }}" title="{{ spec }}"></a></td>
        <td class="label">{{ spec.label }}</td>
        <td>{{ "%.2f"|format(spec.molecule[0].getMolecularWeight() * 1000) }}</td>
    </tr>
    {% if spec.thermo %}
    
    <tr>
     <td>
            <table align="center">
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
                    {{ "%.2f"|format(spec.thermo.getEnthalpy(298) / 4184) }}
                    {% endif %} </td>
                    <td>{% if spec.thermo.Tmin.value_si <= 298 %}
                    {{ "%.2f"|format(spec.thermo.getEntropy(298) / 4.184) }}
                    {% endif %}</td>
                    <td>{{ "%.2f"|format(spec.thermo.getHeatCapacity(300) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec.thermo.getHeatCapacity(500) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec.thermo.getHeatCapacity(1000) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec.thermo.getHeatCapacity(1500) / 4.184) }}</td>
                </tr>
            </table>
        </td></tr>
        
        {% endif %}
    {% endfor %}
</table>

<h2>Reactions ({{ reactions|length }})</h2>

<form id='familySelector' action="">
    <h4>Reaction families:</h4>
{% for family in families %}    <input type="checkbox" id="{{ family|csssafe }}" name="family" value="{{ family|csssafe }}" checked="checked" onclick="updateFamily(this);"><label for="{{ family|csssafe }}">{{ family }} ({{ familyCount[family] }} rxn{{ 's' if familyCount[family] != 1 }})</label><br>
{% endfor %}
    <a href="javascript:checkAllFamilies();" onclick="checkAllFamilies()">check all</a> &nbsp; &nbsp; <a href="javascript:uncheckAllFamilies();" onclick="uncheckAllFamilies();">uncheck all</a><br>

    <h4>Reaction Details:</h4>
    <input type="checkbox" id="kinetics" name="detail" value="kinetics" onclick="updateDetails(this);"><label for="kinetics">Kinetics</label><br>
    <input type="checkbox" id="comment" name="detail" value="comment" onclick="updateDetails(this);"><label for="comment">Comments</label><br>
    <input type="checkbox" id="chemkin" name="detail" value="chemkin" onclick="updateDetails(this);"><label for="chemkin">Chemkin strings</label><br>
    <a href="javascript:checkAllDetails();" onclick="checkAllDetails()">check all</a> &nbsp; &nbsp; <a href="javascript:uncheckAllDetails();" onclick="uncheckAllDetails();">uncheck all</a>
</form>

<h4>Reaction List:</h4>

<table class="reactionList hide_comment hide_kinetics hide_chemkin">
    <tr><th>Index</th><th colspan="3" style="text-align: center;">Reaction</th><th>Family</th></tr>
    {% for rxn in reactions %}
    <tr class="reaction {{ rxn.getSource().label|csssafe }}">
        <td class="index"><a href="{{ rxn.getURL() }}" title="Search on RMG website" class="searchlink">{{ rxn.index }}.</a></td>
        <td class="reactants">{% for reactant in rxn.reactants %}<a href="{{ reactant.molecule[0].getURL() }}"><img src="species_edge/{{ reactant|replace('#','%23') }}.png" alt="{{ reactant }}" title="{{ reactant }}, MW = {{ "%.2f g/mol"|format(reactant.molecule[0].getMolecularWeight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="reactionArrow">{% if rxn.reversible %}&hArr;{% else %}&rarr;{% endif %}</td>
        <td class="products">{% for product in rxn.products %}<a href="{{ product.molecule[0].getURL() }}"><img src="species_edge/{{ product|replace('#','%23') }}.png" alt="{{ product }}" title="{{ product }}, MW = {{ "%.2f g/mol"|format(product.molecule[0].getMolecularWeight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="family">{{ rxn.getSource().label }}</td>
    </tr>
    <tr class="kinetics {{ rxn.getSource().label|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.kinetics.toHTML() }}</td>
    </tr>
    <tr class="chemkin {{ rxn.getSource().label|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.toChemkin(species) }}</td>
    </tr>
    <tr class="comment {{ rxn.getSource().label|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.kinetics.comment }}</td>
    </tr>
    {% endfor %}

</table>

</body>

</html>
""")
        
    f = open(path, 'w')
    f.write(template.render(title=title, species=species, reactions=reactions, families=families, familyCount=familyCount))
    f.close()


def saveDiffHTML(path, commonSpeciesList, speciesList1, speciesList2, commonReactions, uniqueReactions1, uniqueReactions2):
    """
    This function outputs the species and reactions on an HTML page
    for the comparison of two RMG models.
    """
    from model import PDepReaction
    from rmgpy.kinetics import Arrhenius, MultiArrhenius, MultiPDepArrhenius

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

    speciesList = [spec1 for spec1, spec2 in commonSpeciesList] + [spec2 for spec1, spec2 in commonSpeciesList] + speciesList1 + speciesList2
    re_index = re.compile(r'\((\d+)\)$')

        
    if not os.path.isdir(os.path.join(dirname,'species1')):
        os.makedirs(os.path.join(dirname,'species1'))    
    
    if not os.path.isdir(os.path.join(dirname,'species2')):
        os.makedirs(os.path.join(dirname,'species2'))

    for spec1, spec2 in commonSpeciesList:
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
                raise OutputError('{0} species could not be drawn because it did not contain a molecular structure. Please recheck your files.'.format(getSpeciesIdentifier(spec1)))

            
        fstr = os.path.join(dirname, 'species2', '{0}.png'.format(spec2))
        if not os.path.exists(fstr):
            try:
                MoleculeDrawer().draw(spec2.molecule[0], 'png', fstr)
            except IndexError:
                raise OutputError('{0} species could not be drawn because it did not contain a molecular structure. Please recheck your files.'.format(getSpeciesIdentifier(spec2)))
    
                
    for spec in speciesList1:
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
                raise OutputError('{0} species could not be drawn because it did not contain a molecular structure. Please recheck your files.'.format(getSpeciesIdentifier(spec)))
    
    for spec in speciesList2:
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
                raise OutputError('{0} species could not be drawn because it did not contain a molecular structure. Please recheck your files.'.format(getSpeciesIdentifier(spec)))
    


    familyCount1 = {}
    familyCount2 = {}
    for rxn1, rxn2 in commonReactions:
        if isinstance(rxn2.kinetics, (MultiArrhenius,MultiPDepArrhenius)):
            rxn2.duplicate = True   
        if isinstance(rxn1, PDepReaction):
            family = "PDepNetwork"
        else:
            family = rxn1.getSource().label
        if family in familyCount1:
            familyCount1[family] += 1
            familyCount2[family] += 1
        else:
            familyCount1[family] = 1
            familyCount2[family] = 1

    for rxn in uniqueReactions1:
        if isinstance(rxn, PDepReaction):
            family = "PDepNetwork"
        else:
            family = rxn.getSource().label
        if family in familyCount1:
            familyCount1[family] += 1
        else:
            familyCount1[family] = 1

    for rxn in uniqueReactions2:
        if isinstance(rxn, PDepReaction):
            family = "PDepNetwork"
        else:
            family = rxn.getSource().label
        if family in familyCount2:
            familyCount2[family] += 1
        else:
            familyCount2[family] = 1

    families1 = familyCount1.keys()
    families2 = familyCount2.keys()
    families1.sort()
    families2.sort()



    ## jinja2 filters etc.
    to_remove_from_css_names = re.compile('[/.\-+,]')
    def csssafe(input):
        "Replace unsafe CSS class name characters with an underscore."
        return to_remove_from_css_names.sub('_',input)

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
            width: 100%;
            border-collapse: collapse;
        }
        table.speciesList th, table.reactionList th {
            text-align: left;
        }
        tr.reaction td {
            border-top: 1px solid #808080;
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
        tr.comment {
            font-size: small;
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

        .hide_comment .comment{
            display: none !important;
        }
        .hide_kinetics .kinetics{
           display: none !important;
        }
        .hide_chemkin .chemkin{
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
    $(document).ready(function() {
        checkAllFamilies();
        uncheckAllDetails();
    });

    </script>
</head>

<body>

<h1>{{ title }}</h1>

<h2 align="center">Common Species ({{ commonSpecies|length }})</h2>

<table width="100%">
    {% for spec1, spec2 in commonSpecies %}
    <tr>
        <td width="100%" colspan="4">
            <table align="center">
                <tr class="species">
                    <td>{{ spec1.label }}</td>
                    <td class="structure" align="center"><a href="{{spec1.molecule[0].getURL()}}"><img src="species1/{{ spec1|replace('#','%23') }}.png" alt="{{ spec1 }}" title="{{ spec1 }}"></a></td>
                    <td>{{ "%.2f"|format(spec1.molecule[0].getMolecularWeight() * 1000) }}</td>
                </tr>
            </table>
        </td>
    </tr>
    {% if spec1.thermo and spec2.thermo %}
    {% if spec1.thermo.isIdenticalTo(spec2.thermo) %}
    <tr width=100%>
         <td colspan="4" valign="top" width=50%><div align="center"><font color="blue">IDENTICAL THERMO WAS FOUND FOR THIS SPECIES.</font></div>
    </tr>
    {% elif spec1.thermo.isSimilarTo(spec2.thermo) %}
    <tr width=100%>
         <td colspan="4" valign="top" width=50%><div align="center"><font color="green">SIMILAR THERMO WAS FOUND FOR THIS SPECIES.</font></div>
    </tr>
    {% else %}
     <tr width=100%>
         <td colspan="4" valign="top" width=50%><div align="center"><font color="red">DIFFERENT THERMO WAS FOUND FOR THIS SPECIES.</font></div>
    </tr>
    {% endif%}
    <tr>
        <td width="10%">{{ spec1.index }}. </td>
        <td width="40%">
            <table width="100%">
                <tr>
                    <th>H300</th>
                    <th>S300</th>
                    <th>Cp300</th>
                    <th>Cp500</th>
                    <th>Cp1000</th>
                    <th>Cp1500</th>
                </tr>
                <tr>
                    <td>{% if spec1.thermo.Tmin.value_si <= 300 %}
                    {{ "%.2f"|format(spec1.thermo.getEnthalpy(300) / 4184) }}
                    {% endif %}</td>
                    <td>
                    {% if spec1.thermo.Tmin.value_si <= 300 %}
                    {{ "%.2f"|format(spec1.thermo.getEntropy(300) / 4.184) }}
                    {% endif %}</td>
                    <td>{{ "%.2f"|format(spec1.thermo.getHeatCapacity(300) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec1.thermo.getHeatCapacity(500) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec1.thermo.getHeatCapacity(1000) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec1.thermo.getHeatCapacity(1500) / 4.184) }}</td>
                </tr>
            </table>
        </td>
        <td width="10%">{{ spec2.index }}.</td>
        <td width="40%">
            <table width="100%">
                <tr>
                    <th>H300</th>
                    <th>S300</th>
                    <th>Cp300</th>
                    <th>Cp500</th>
                    <th>Cp1000</th>
                    <th>Cp1500</th>
                </tr>
                <tr>
                    <td>{% if spec2.thermo.Tmin.value_si <= 300 %}
                    {{ "%.2f"|format(spec2.thermo.getEnthalpy(300) / 4184) }}
                    {% endif %}</td>
                    <td>{% if spec2.thermo.Tmin.value_si <= 300 %}
                    {{ "%.2f"|format(spec2.thermo.getEntropy(300) / 4.184) }}
                    {% endif %}</td>
                    <td>{{ "%.2f"|format(spec2.thermo.getHeatCapacity(300) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec2.thermo.getHeatCapacity(500) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec2.thermo.getHeatCapacity(1000) / 4.184) }}</td>
                    <td>{{ "%.2f"|format(spec2.thermo.getHeatCapacity(1500) / 4.184) }}</td>
                </tr>
            </table>
        </td>
    </tr>
    {% endif %}
    <tr>
        <td width="100%" colspan="4"><hr/></td>
    </tr>
    {% endfor %}
</table>

<table width=100%>
<tr colspan="2">
<td width=50% valign="top">

<h2>Model 1: Unique Species ({{ speciesList1|length }})</h2>
<table class="speciesList">
    <tr><th>Index</th><th>Structure</th><th>Label</th><th>Mol. Wt. (g/mol)</th></tr>
    {% for spec in speciesList1 %}
    <tr class="species">
        <td class="index">
        {{ spec.index }}.</td>
        <td class="structure"><a href="{{ spec.molecule[0].getURL() }}"><img src="species1/{{ spec|replace('#','%23') }}.png" alt="{{ spec }}" title="{{ spec }}"></a></td>
        <td class="label">{{ spec.label }}</td>
        <td>{{ "%.2f"|format(spec.molecule[0].getMolecularWeight() * 1000) }}</td>
    </tr>
    {% endfor %}
</table>
</td>
<td width=50% valign="top">
<h2>Model 2: Unique Species ({{ speciesList2|length }})</h2>
<table class="speciesList">
    <tr><th>Index</th><th>Structure</th><th>Label</th><th>Mol. Wt. (g/mol)</th></tr>
    {% for spec in speciesList2 %}
    <tr class="species">
        <td class="index">
        {{ spec.index }}.</td>
        <td class="structure"><a href="{{ spec.molecule[0].getURL() }}"><img src="species2/{{ spec|replace('#','%23') }}.png" alt="{{ spec }}" title="{{ spec }}"></a></td>
        <td class="label">{{ spec.label }}</td>
        <td>{{ "%.2f"|format(spec.molecule[0].getMolecularWeight() * 1000) }}</td>
    </tr>
    {% endfor %}
</table>
</td></tr>

<tr colspan="2">
<td width=50% valign="top">
<h2>Model 1 Reactions ({{ commonReactions|length + uniqueReactions1|length}})</h2>

<form id='familySelector' action="">
    <h4>Reaction families:</h4>
{% for family in families1 %}    <input type="checkbox" id="{{ family|csssafe }}" name="family" value="{{ family|csssafe }}" checked="checked" onclick="updateFamily(this);"><label for="{{ family|csssafe }}">{{ family }} ({{ familyCount1[family] }} rxn{{ 's' if familyCount1[family] != 1 }})</label><br>
{% endfor %}
    <a href="javascript:checkAllFamilies();" onclick="checkAllFamilies()">check all</a> &nbsp; &nbsp; <a href="javascript:uncheckAllFamilies();" onclick="uncheckAllFamilies();">uncheck all</a><br>

    <h4>Reaction Details:</h4>
    <input type="checkbox" id="kinetics" name="detail" value="kinetics" onclick="updateDetails(this);"><label for="kinetics">Kinetics</label><br>
    <input type="checkbox" id="comment" name="detail" value="comment" onclick="updateDetails(this);"><label for="comment">Comments</label><br>
    <input type="checkbox" id="chemkin" name="detail" value="chemkin" onclick="updateDetails(this);"><label for="chemkin">Chemkin strings</label><br>
    <a href="javascript:checkAllDetails();" onclick="checkAllDetails()">check all</a> &nbsp; &nbsp; <a href="javascript:uncheckAllDetails();" onclick="uncheckAllDetails();">uncheck all</a>
</form>
</td>

<td width=50% valign="top">
<h2>Model 2 Reactions ({{ commonReactions|length +uniqueReactions2|length}})</h2>

<form id='familySelector' action="">
    <h4>Reaction families:</h4>
{% for family in families2 %}    <input type="checkbox" id="{{ family|csssafe }}" name="family" value="{{ family|csssafe }}" checked="checked" onclick="updateFamily(this);"><label for="{{ family|csssafe }}">{{ family }} ({{ familyCount2[family] }} rxn{{ 's' if familyCount2[family] != 1 }})</label><br>
{% endfor %}
    <a href="javascript:checkAllFamilies();" onclick="checkAllFamilies()">check all</a> &nbsp; &nbsp; <a href="javascript:uncheckAllFamilies();" onclick="uncheckAllFamilies();">uncheck all</a><br>

    <h4>Reaction Details:</h4>
    <input type="checkbox" id="kinetics" name="detail" value="kinetics" onclick="updateDetails(this);"><label for="kinetics">Kinetics</label><br>
    <input type="checkbox" id="comment" name="detail" value="comment" onclick="updateDetails(this);"><label for="comment">Comments</label><br>
    <input type="checkbox" id="chemkin" name="detail" value="chemkin" onclick="updateDetails(this);"><label for="chemkin">Chemkin strings</label><br>
    <a href="javascript:checkAllDetails();" onclick="checkAllDetails()">check all</a> &nbsp; &nbsp; <a href="javascript:uncheckAllDetails();" onclick="uncheckAllDetails();">uncheck all</a>
</form>
</td>
</tr>

</table>


<table width=100%>
<tr><td width=100% align="center">
<h2>Common Reactions ({{ commonReactions|length}})</h2></td></tr>


<tr colspan="1"><td width=100%>

<table class="reactionList hide_comment hide_kinetics hide_chemkin" width=100% cellpadding="10">
    <tr colspan="4" width=100%><th>Index.</th><th>Family</th><th>Index.</th><th>Family</th></tr>

    {% for rxn1, rxn2 in commonReactions %}

<tr>
<td width=100% colspan="4"><hr>
<table align="center">
<tr>
    <td class="reactants" align="right">{% for reactant in rxn1.reactants %}<a href="{{reactant.molecule[0].getURL() }}"><img src="species1/{{ reactant|replace('#','%23') }}.png" alt="{{ reactant }}" title="{{ reactant }}, MW = {{ "%.2f"|format(reactant.molecule[0].getMolecularWeight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
    <td class="reactionArrow" align="center">{% if rxn1.reversible %}&hArr;{% else %}&rarr;{% endif %}</td>
    <td class="products" align="left">{% for product in rxn1.products %}<a href="{{product.molecule[0].getURL()}}"><img src="species1/{{ product|replace('#','%23') }}.png" alt="{{ product }}" title="{{ product }}, MW = {{ "%.2f"|format(product.molecule[0].getMolecularWeight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
</tr>
</table>
</td>
</tr>


{% if rxn1.kinetics.isIdenticalTo(rxn2.kinetics) %}

 <tr width=100%>
     <td colspan="4" valign="top" width=50%><div align="center"><font color="blue">IDENTICAL KINETICS WERE FOUND FOR THIS REACTION.</font></div>

</tr>
{% elif rxn1.kinetics.isSimilarTo(rxn2.kinetics) %}

 <tr width=100%>
     <td colspan="4" valign="top" width=50%><div align="center"><font color="green">SIMILAR KINETICS WERE FOUND FOR THIS REACTION.</font></div>

</tr>

{% else %}
     <tr width=100%>
         <td colspan="4" valign="top" width=50%><div align="center"><font color="red">DIFFERENT KINETICS WERE FOUND FOR THIS REACTION.</font></div>
    </tr>
{% endif%}


<tr width=100%>
     <td class="index" width=10%><a href="{{ rxn1.getURL() }}" title="Search on RMG website" class="searchlink">{{ rxn1.index }}.</a></td>
     <td class="family" width=40%>{{ rxn1.getSource().label }}</td>

     <td class="index" width=10%><a href="{{ rxn2.getURL() }}" title="Search on RMG website" class="searchlink">{{ rxn2.index }}.</a></td>
     <td class="family" width=40%>{{ rxn2.getSource().label }}</td>
 </tr>

<tr width=100%>{% if not rxn1.isIsomorphic(rxn2, eitherDirection=False) %} 
<td colspan="2" width=50%></td>
<td colspan="2" width=50%>* Reaction was found in reverse 


{% if not rxn2.duplicate %}
<P><b>Fitted Reverse Kinetics:</b>
{% if not rxn2.kinetics.isPressureDependent() %}
{{rxn2.generateReverseRateCoefficient().toHTML() }}
{% else %} Pressure dependent
{% endif %}
{% endif %}

<P><b>Original Kinetics:</b>
{% endif %}</td>
</tr>

<tr width=100%>
     <td colspan="2" valign="top" width=50%>
     
     {{ rxn1.kinetics.toHTML() }}</td>
     <td colspan="2" valign="top" width=50%>
     {{ rxn2.kinetics.toHTML() }}</td>
</tr>

<tr width=100%>
    <td colspan="2" valign="top" width=50%><font size="1pt" face="courier">{{ rxn1.toChemkin(speciesList) }}</font></td>
    <td colspan="2" valign="top" width=50%><font size="1pt" face="courier">{{ rxn2.toChemkin(speciesList) }}</font></td>
</tr>

<tr width=100%>
    <td colspan="2" valign="top" width=50%><font size="1pt">{{ rxn1.kinetics.comment }}</font></td>
    <td colspan="2" valign="top" width=50%><font size="1pt">{{ rxn2.kinetics.comment }}</font></td>
</tr>
{% endfor %}

</table>
</td></tr></table>

<table>
<tr>
<td width=50% valign="top">
<h2>Model 1: Unique Reactions ({{ uniqueReactions1|length}})</h2>

<table class="reactionList hide_comment hide_kinetics hide_chemkin">
    <tr><th>Index</th><th colspan="3" style="text-align: center;">Reaction</th><th>Family</th></tr>
    {% for rxn in uniqueReactions1 %}
    <tr class="reaction {{ rxn.getSource().label|csssafe }}">
        <td class="index"><a href="{{ rxn.getURL() }}" title="Search on RMG website" class="searchlink">{{ rxn.index }}.</a></td>
        <td class="reactants">{% for reactant in rxn.reactants %}<a href="{{ reactant.molecule[0].getURL() }}"><img src="species1/{{ reactant|replace('#','%23') }}.png" alt="{{ reactant }}" title="{{ reactant }}, MW = {{ "%.2f"|format(reactant.molecule[0].getMolecularWeight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="reactionArrow">{% if rxn.reversible %}&hArr;{% else %}&rarr;{% endif %}</td>
        <td class="products">{% for product in rxn.products %}<a href="{{ product.molecule[0].getURL() }}"><img src="species1/{{ product|replace('#','%23') }}.png" alt="{{ product }}" title="{{ product }}, MW = {{ "%.2f"|format(product.molecule[0].getMolecularWeight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="family">{{ rxn.getSource().label }}</td>
    </tr>
    <tr class="kinetics {{ rxn.getSource().label|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.kinetics.toHTML() }}</td>
    </tr>
    <tr class="chemkin {{ rxn.getSource().label|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.toChemkin(species) }}</td>
    </tr>
    <tr class="comment {{ rxn.getSource().label|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.kinetics.comment }}</td>
    </tr>
    {% endfor %}
    </table>

</td>

<td width=50% valign="top">
<h2>Model 2: Unique Reactions ({{ uniqueReactions2|length}})</h2>

<table class="reactionList hide_comment hide_kinetics hide_chemkin">
    <tr><th>Index</th><th colspan="3" style="text-align: center;">Reaction</th><th>Family</th></tr>
    {% for rxn in uniqueReactions2 %}
    <tr class="reaction {{ rxn.getSource().label|csssafe }}">
        <td class="index"><a href="{{ rxn.getURL() }}" title="Search on RMG website" class="searchlink">{{ rxn.index }}.</a></td>
        <td class="reactants">{% for reactant in rxn.reactants %}<a href="{{ reactant.molecule[0].getURL() }}"><img src="species2/{{ reactant|replace('#','%23') }}.png" alt="{{ reactant }}" title="{{ reactant }}, MW = {{ "%.2f"|format(reactant.molecule[0].getMolecularWeight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="reactionArrow">{% if rxn.reversible %}&hArr;{% else %}&rarr;{% endif %}</td>
        <td class="products">{% for product in rxn.products %}<a href="{{ product.molecule[0].getURL() }}"><img src="species2/{{ product|replace('#','%23') }}.png" alt="{{ product }}" title="{{ product }}, MW = {{ "%.2f"|format(product.molecule[0].getMolecularWeight() * 1000) }}"></a>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="family">{{ rxn.getSource().label }}</td>
    </tr>
    <tr class="kinetics {{ rxn.getSource().label|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.kinetics.toHTML() }}</td>
    </tr>
    <tr class="chemkin {{ rxn.getSource().label|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.toChemkin(species) }}</td>
    </tr>
    <tr class="comment {{ rxn.getSource().label|csssafe }}">
        <td></td>
        <td colspan="4">{{ rxn.kinetics.comment }}</td>
    </tr>
    {% endfor %}
    </table>

</td>

</table>



</body>

</html>
""")
    f = open(path, 'w')
    f.write(template.render(title=title, commonSpecies=commonSpeciesList, speciesList1=speciesList1, speciesList2 = speciesList2, commonReactions=commonReactions, uniqueReactions1=uniqueReactions1, uniqueReactions2=uniqueReactions2, families1=families1, families2=families2, familyCount1=familyCount1,familyCount2=familyCount2, speciesList=speciesList))
    f.close()
