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

################################################################################

class OutputError(Exception):
    """
    This exception is raised whenever an error occurs while saving output
    information. Pass a string describing the circumstances of the exceptional
    behavior.
    """
    pass

################################################################################

def saveOutputHTML(path, reactionModel):
    """
    Save the current set of core species and reactions of `reactionModel` to
    an HTML file `path` on disk. As part of this process, drawings of all core
    species are created in the species folder (if they don't already exist)
    using the :mod:`rmgpy.chem.ext.molecule_draw` module. The :mod:`jinja`
    package is used to generate the HTML; if this package is not found, no
    HTML will be generated (but the program will carry on).
    """

    from model import PDepReaction
    
    from rmgpy.molecule_draw import drawMolecule
    try:
        import jinja2
    except ImportError:
        logging.warning("jinja package not found; HTML output will not be saved.")
        return

    path = os.path.abspath(path)
    dirname = os.path.dirname(path)

    # Draw molecules if necessary
    for spec in reactionModel.core.species:
        fstr = os.path.join(dirname, 'species', '{0}.png'.format(spec))
        if not os.path.exists(fstr):
            drawMolecule(spec.molecule[0], fstr)

    # Prepare parameters to pass to jinja template
    title = 'RMG Output'

    species = reactionModel.core.species[:]
    species.sort(key=lambda x: x.index)

    reactions = [rxn for rxn in reactionModel.core.reactions if not isinstance(rxn, PDepReaction)]
    reactions.sort(key=lambda x: x.index)

    pdepreactions = [rxn for rxn in reactionModel.core.reactions if isinstance(rxn, PDepReaction)]
    pdepreactions.sort(key=lambda x: x.index)

    familyCount = {}
    for rxn in reactions:
        family = rxn.getSource()
        if family in familyCount:
            familyCount[family] += 1
        else:
            familyCount[family] = 1
    families = familyCount.keys()
    families.sort(key=lambda x: x.label)
    
    # Make HTML file
    template = jinja2.Template(
"""<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN">
<html lang="en">

<head>
    <title>{{ title }}</title>
    <style>
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
            # width: 100%;
            border-collapse: collapse;
        }
        table.speciesList th, table.reactionList th {
            text-align: left;
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
            border-bottom: 1px solid #808080;
        }
        tr.kinetics {
            font-size: small;
        }
        .KineticsData {
            # border: 1px solid grey;
        }
        .KineticsData th {
            width: 15em;
            word-wrap: none;
        }
        .KineticsData td {
            width: 3em;
        }
        .searchlink {
            font-size: large;
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

    function checkAllFamilies() {
        $("#familyForm").find("input").each(function() { this.checked = true; updateFamily(this); });
    }

    function uncheckAllFamilies() {
        $("#familyForm").find("input").each(function() { this.checked = false; updateFamily(this); });
    }

    $(document).ready(function() {
        checkAllFamilies();
    });

    </script>
</head>

<body>

<h1>{{ title }}<h1>

<h2>Species ({{ species|length }})</h2>

<table class="speciesList">
    <tr><th>Index</th><th>Structure</th><th>Label</th></tr>
    {% for spec in species %}
    <tr class="species">
        <td class="index">{{ spec.index }}.</td>
        <td class="structure"><img src="species/{{ spec }}.png" alt="{{ spec }}" title="{{ spec }}"/></td>
        <td class="label">{{ spec.label }}</td>
    </tr>
    {% endfor %}
</table>

<h2>Reactions ({{ reactions|length }})</h2>

<p>
<form id='familySelector'>
    <div><input type="checkbox" id="kinetics" name="family" value="kinetics" checked="checked" onclick="updateFamily(this);"/><label for="kinetics">Kinetics</label></div>
    <div><input type="checkbox" id="comment" name="family" value="comment" checked="checked" onclick="updateFamily(this);"/><label for="comment">Comments</label></div>
    <div>Reaction families:</div>
{% for family in families %}    <div><input type="checkbox" id="{{ family.label }}" name="family" value="{{ family.label }}" checked="checked" onclick="updateFamily(this);"/><label for="{{ family.label }}">{{ family.label }} ({{ familyCount[family] }})</label></div>
{% endfor %}    <div><a href="#" onclick="checkAllFamilies();">check all</a> &nbsp; &nbsp; <a href="#" onclick="uncheckAllFamilies();">uncheck all</a></div>
</form>
</p>

<table class="reactionList">
    <tr><th>Index</th><th colspan="3" style="text-align: center;">Reaction</th><th>Family</th></tr>
    {% for rxn in reactions %}
    <tr class="reaction {{ rxn.getSource().label }}">
        <td class="index">{{ rxn.index }}.</td>
        <td class="reactants">{% for reactant in rxn.reactants %}<img src="species/{{ reactant }}.png" alt="{{ reactant }}" title="{{ reactant }}"/>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="reactionArrow">{% if rxn.reversible %}&hArr;{% else %}&rarr;{% endif %}</td>
        <td class="products">{% for product in rxn.products %}<img src="species/{{ product }}.png" alt="{{ product }}" title="{{ product }}"/>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="family">{{ rxn.getSource().label }}</td>
    </tr>
    <tr class="kinetics">
        <td>
        <a href="{{ rxn.getURL() }}" alt="Search on RMG website" class="searchlink">?</a>
        </td>
        <td colspan="4">{{ rxn.kinetics.toHTML() }}</td>
    </tr>
    <tr class="comment">
        <td></td>
        <td colspan="4">{{ rxn.kinetics.comment }}</td>
    </tr>
    {% endfor %}
    {% for rxn in pdepreactions %}
    <tr class="reaction {{ pdepnetreaction }}">
        <td class="index">{{ rxn.index }}.</td>
        <td class="reactants">{% for reactant in rxn.reactants %}<img src="species/{{ reactant }}.png" alt="{{ reactant }}" title="{{ reactant }}"/>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="reactionArrow">{% if rxn.reversible %}&hArr;{% else %}&rarr;{% endif %}</td>
        <td class="products">{% for product in rxn.products %}<img src="species/{{ product }}.png" alt="{{ product }}" title="{{ product }}"/>{% if not loop.last %} + {% endif %}{% endfor %}</td>
        <td class="family"></td>
    </tr>
    <tr class="kinetics">
        <td>
        <a href="{{ rxn.getURL() }}" alt="Search on RMG website" class="searchlink">?</a>
        </td>
        <td colspan="4">{{ rxn.kinetics.toHTML() }}</td>
    </tr>
    <tr class="comment">
        <td></td>
        <td colspan="4">{{ rxn.kinetics.comment }}</td>
    </tr>
    {% endfor %}
</table>

</body>

</html>
""")
    f = open(path, 'w')
    f.write(template.render(title=title, species=species, reactions=reactions, pdepreactions=pdepreactions, families=families, familyCount=familyCount))
    f.close()
