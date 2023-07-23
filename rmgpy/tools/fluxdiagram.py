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
This module contains functions for generating flux diagrams.
"""

import math
import os.path
import re
import warnings

import numpy as np
import pydot

from rmgpy.kinetics.diffusionLimited import diffusion_limiter
from rmgpy.rmg.settings import SimulatorSettings
from rmgpy.solver.base import TerminationTime, TerminationConversion
from rmgpy.solver.liquid import LiquidReactor
from rmgpy.tools.loader import load_rmg_job

################################################################################

# Here you can set the default values for options that control the generated
# flux diagrams.

# Options controlling the individual flux diagram renderings:
program = "dot"  # The program to use to lay out the nodes and edges
max_node_count = 50  # The maximum number of nodes to show in the diagram
max_edge_count = 50  # The maximum number of edges to show in the diagram
concentration_tol = 1e-6  # The lowest fractional concentration to show (values below this will appear as zero)
species_rate_tol = 1e-6  # The lowest fractional species rate to show (values below this will appear as zero)
max_node_pen_width = 10.0  # The thickness of the border around a node at maximum concentration
max_edge_pen_width = 10.0  # The thickness of the edge at maximum species rate
radius = 1  # The graph radius to plot around a central species
central_reaction_count = None  # The maximum number of reactions to draw from each central species (None draws all)
# If radius > 1, then this is the number of reactions from every species

# Options controlling the ODE simulations:
initial_time = 1e-12  # The time at which to initiate the simulation, in seconds
time_step = 10**0.1  # The multiplicative factor to use between consecutive time points
abs_tol = 1e-16  # The absolute tolerance to use in the ODE simluations
rel_tol = 1e-8  # The relative tolerance to use in the ODE simulations

# Options controlling the generated movie:
video_fps = 6  # The number of frames per second in the generated movie
initial_padding = 5  # The number of seconds to display the initial fluxes at the start of the video
final_padding = 5  # The number of seconds to display the final fluxes at the end of the video


################################################################################


def generate_flux_diagram(
    reaction_model,
    times,
    concentrations,
    reaction_rates,
    output_directory,
    central_species_list=None,
    superimpose=False,
    species_directory=None,
    settings=None,
):
    """
    For a given `reaction_model` and simulation results stored as arrays of
    `times`, species `concentrations`, and `reaction_rates`, generate a series
    of flux diagrams as frames of an animation, then stitch them together into
    a movie. The individual frames and the final movie are saved on disk at
    `output_directory.`
    """
    global max_node_count, max_edge_count, concentration_tol, species_rate_tol, max_node_pen_width, max_edge_pen_width, radius, central_reaction_count
    # Allow user defined settings for flux diagram generation if given
    if settings:
        max_node_count = settings.get("max_node_count", max_node_count)
        max_edge_count = settings.get("max_edge_count", max_edge_count)
        concentration_tol = settings.get("concentration_tol", concentration_tol)
        species_rate_tol = settings.get("species_rate_tol", species_rate_tol)
        max_node_pen_width = settings.get("max_node_pen_width", max_node_pen_width)
        max_edge_pen_width = settings.get("max_edge_pen_width", max_edge_pen_width)
        radius = settings.get("radius", radius)
        central_reaction_count = settings.get("central_reaction_count", central_reaction_count)

    # Get the species and reactions corresponding to the provided concentrations and reaction rates
    species_list = reaction_model.core.species[:]
    num_species = len(species_list)
    reaction_list = reaction_model.core.reactions[:]

    # Search for indices of central species
    central_species_indices = []
    if central_species_list is not None:
        for centralSpecies in central_species_list:
            for i, species in enumerate(species_list):
                if species.index == centralSpecies:
                    central_species_indices.append(i)
                    break
            else:
                raise Exception("Central species '{}' could not be found in species list.".format(centralSpecies))

    # Compute the rates between each pair of species (big matrix warning!)
    species_rates = np.zeros((len(times), num_species, num_species), float)
    for index, reaction in enumerate(reaction_list):
        rate = reaction_rates[:, index]
        if not reaction.pairs:
            reaction.generate_pairs()
        for reactant, product in reaction.pairs:
            reactant_index = species_list.index(reactant)
            product_index = species_list.index(product)
            species_rates[:, reactant_index, product_index] += rate
            species_rates[:, product_index, reactant_index] -= rate

    # Determine the maximum concentration for each species and the maximum overall concentration
    max_concentrations = np.max(np.abs(concentrations), axis=0)
    max_concentration = np.max(max_concentrations)

    # Determine the maximum reaction rates
    max_reaction_rates = np.max(np.abs(reaction_rates), axis=0)

    # Determine the maximum rate for each species-species pair and the maximum overall species-species rate
    max_species_rates = np.max(np.abs(species_rates), axis=0)
    max_species_rate = np.max(max_species_rates)
    species_index = max_species_rates.reshape((num_species * num_species)).argsort()

    # Determine the nodes and edges to keep
    nodes = []
    edges = []
    if not superimpose and central_species_list is not None:
        for central_species_index in central_species_indices:
            nodes.append(central_species_index)
            add_adjacent_nodes(
                central_species_index,
                nodes,
                edges,
                species_list,
                reaction_list,
                max_reaction_rates,
                max_species_rates,
                reactionCount=central_reaction_count,
                rad=radius,
            )
    else:
        for i in range(num_species * num_species):
            product_index, reactant_index = divmod(species_index[-i - 1], num_species)
            if reactant_index > product_index:
                # Both reactant -> product and product -> reactant are in this list,
                # so only keep one of them
                continue
            if max_species_rates[reactant_index, product_index] == 0:
                break
            if reactant_index not in nodes and len(nodes) < max_node_count:
                nodes.append(reactant_index)
            if product_index not in nodes and len(nodes) < max_node_count:
                nodes.append(product_index)
            if [reactant_index, product_index] not in edges and [product_index, reactant_index] not in edges:
                edges.append([reactant_index, product_index])
            if len(nodes) > max_node_count:
                break
            if len(edges) >= max_edge_count:
                break

        if superimpose and central_species_list is not None:
            nodes_copy = nodes[:]
            for central_species_index in central_species_indices:
                if central_species_index not in nodes:  # Only add central species if it doesn't already exist
                    nodes.append(central_species_index)
                    # Recursively add nodes until they connect with main graph
                    add_adjacent_nodes(
                        central_species_index,
                        nodes,
                        edges,
                        species_list,
                        reaction_list,
                        max_reaction_rates,
                        max_species_rates,
                        reactionCount=central_reaction_count,
                        rad=-1,  # "-1" signifies that we add nodes until they connect to the main graph
                        mainNodes=nodes_copy,
                    )

    # Create the master graph
    # First we're going to generate the coordinates for all of the nodes; for
    # this we use the thickest pen widths for all nodes and edges
    graph = pydot.Dot("flux_diagram", graph_type="digraph", overlap="false")
    graph.set_rankdir("LR")
    graph.set_fontname("sans")
    graph.set_fontsize("10")

    # Add a node for each species
    for index in nodes:
        species = species_list[index]
        node = pydot.Node(name=str(species))
        node.set_penwidth(max_node_pen_width)
        graph.add_node(node)
        # Try to use an image instead of the label
        species_index = str(species) + ".png"
        image_path = ""
        if not species_directory or not os.path.exists(species_directory):
            continue
        for root, dirs, files in os.walk(species_directory):
            for f in files:
                if f.endswith(species_index):
                    image_path = os.path.join(root, f)
                    break
        if os.path.exists(image_path):
            node.set_image(image_path)
            node.set_label(" ")
    # Add an edge for each species-species rate
    for reactant_index, product_index in edges:
        if reactant_index in nodes and product_index in nodes:
            reactant = species_list[reactant_index]
            product = species_list[product_index]
            edge = pydot.Edge(str(reactant), str(product))
            edge.set_penwidth(max_edge_pen_width)
            graph.add_edge(edge)

    # Generate the coordinates for all of the nodes using the specified program
    graph = pydot.graph_from_dot_data(graph.create_dot(prog=program).decode("utf-8"))[0]

    # Now iterate over the time points, setting the pen widths appropriately
    # This should preserve the coordinates of the nodes from frame to frame
    frame_number = 1
    for t in range(len(times)):
        # Update the nodes
        slope = -max_node_pen_width / math.log10(concentration_tol)
        for index in nodes:
            species = species_list[index]
            if re.search(r"^[a-zA-Z0-9_]*$", str(species)) is not None:
                species_string = str(species)
            else:
                # species name contains special characters
                species_string = '"{0}"'.format(str(species))

            node = graph.get_node(species_string)[0]
            concentration = concentrations[t, index] / max_concentration
            if concentration < concentration_tol:
                penwidth = 0.0
            else:
                penwidth = round(slope * math.log10(concentration) + max_node_pen_width, 3)
            node.set_penwidth(penwidth)
        # Update the edges
        slope = -max_edge_pen_width / math.log10(species_rate_tol)
        for index in range(len(edges)):
            reactant_index, product_index = edges[index]
            if reactant_index in nodes and product_index in nodes:
                reactant = species_list[reactant_index]
                product = species_list[product_index]

                if re.search(r"^[a-zA-Z0-9_]*$", str(reactant)) is not None:
                    reactant_string = str(reactant)
                else:
                    reactant_string = '"{0}"'.format(str(reactant))

                if re.search(r"^[a-zA-Z0-9_]*$", str(product)) is not None:
                    product_string = str(product)
                else:
                    product_string = '"{0}"'.format(str(product))

                edge = graph.get_edge(reactant_string, product_string)[0]
                # Determine direction of arrow based on sign of rate
                species_rate = species_rates[t, reactant_index, product_index] / max_species_rate
                if species_rate < 0:
                    edge.set_dir("back")
                    species_rate = -species_rate
                else:
                    edge.set_dir("forward")
                # Set the edge pen width
                if species_rate < species_rate_tol:
                    penwidth = 0.0
                    edge.set_dir("none")
                else:
                    penwidth = round(slope * math.log10(species_rate) + max_edge_pen_width, 3)
                edge.set_penwidth(penwidth)
        # Save the graph at this time to a dot file and a PNG image
        if times[t] == 0:
            label = "t = 0 s"
        else:
            label = "t = 10^{0:.1f} s".format(math.log10(times[t]))
        graph.set_label(label)
        if t == 0:
            repeat = video_fps * initial_padding
        elif t == len(times) - 1:
            repeat = video_fps * final_padding
        else:
            repeat = 1
        for r in range(repeat):
            graph.write_dot(os.path.join(output_directory, "flux_diagram_{0:04d}.dot".format(frame_number)))
            graph.write_png(os.path.join(output_directory, "flux_diagram_{0:04d}.png".format(frame_number)))
            frame_number += 1

    # Use ffmpeg to stitch the PNG images together into a movie
    import subprocess

    command = [
        "ffmpeg",
        "-framerate",
        "{0:d}".format(video_fps),  # Duration of each image
        "-i",
        "flux_diagram_%04d.png",  # Input file format
        "-c:v",
        "mpeg4",  # Encoder
        "-r",
        "30",  # Video framerate
        "-pix_fmt",
        "yuv420p",  # Pixel format
        "flux_diagram.avi",
    ]  # Output filename

    subprocess.check_call(command, cwd=output_directory)


################################################################################


def add_adjacent_nodes(
    targetNodeIndex, nodes, edges, species_list, reaction_list, maxReactionRates, maxSpeciesRates, reactionCount=None, rad=0, mainNodes=None
):
    """
    Add adjacent nodes in flux diagram up to a certain radius or
    until they connect with the main graph. Radius should be set to a
    negative value in the latter case.
    """
    if rad == 0:  # Base case if using radius
        return
    elif rad < 0 and targetNodeIndex in mainNodes:  # Base case if connecting to main graph
        return
    else:  # Recurse until all nodes up to desired radius have been added or until they connect to the main graph
        # Select all reactions involving target node
        target_reactions_indices = []
        for index, reaction in enumerate(reaction_list):
            reactant_indices = [species_list.index(reactant) for reactant in reaction.reactants]
            product_indices = [species_list.index(product) for product in reaction.products]
            if targetNodeIndex in reactant_indices or targetNodeIndex in product_indices:
                target_reactions_indices.append(index)

        # Sort by maximum reaction rates and only extract top reactions if desired
        target_reactions_indices.sort(key=lambda index: maxReactionRates[index], reverse=True)
        if reactionCount is None:
            target_reaction_list = [reaction_list[index] for index in target_reactions_indices]
        else:
            target_reaction_list = [reaction_list[index] for i, index in enumerate(target_reactions_indices) if i < reactionCount]

        for reaction in target_reaction_list:
            for reactant, product in reaction.pairs:
                reactant_index = species_list.index(reactant)
                product_index = species_list.index(product)
                if reactant_index == targetNodeIndex:
                    if product_index not in nodes:
                        nodes.append(product_index)
                        add_adjacent_nodes(
                            product_index,
                            nodes,
                            edges,
                            species_list,
                            reaction_list,
                            maxReactionRates,
                            maxSpeciesRates,
                            reactionCount=reactionCount,
                            rad=rad - 1,
                            mainNodes=mainNodes,
                        )
                    if [reactant_index, product_index] not in edges and [product_index, reactant_index] not in edges:
                        edges.append([reactant_index, product_index])
                if product_index == targetNodeIndex:
                    if reactant_index not in nodes:
                        nodes.append(reactant_index)
                        add_adjacent_nodes(
                            reactant_index,
                            nodes,
                            edges,
                            species_list,
                            reaction_list,
                            maxReactionRates,
                            maxSpeciesRates,
                            reactionCount=reactionCount,
                            rad=rad - 1,
                            mainNodes=mainNodes,
                        )
                    if [reactant_index, product_index] not in edges and [product_index, reactant_index] not in edges:
                        edges.append([reactant_index, product_index])


################################################################################


def simulate(reaction_model, reaction_system, settings=None):
    """
    Generate and return a set of core and edge species and reaction fluxes
    by simulating the given `reaction_system` using the given `reaction_model`.
    """
    global time_step
    # Allow user defined settings for flux diagram generation if given
    if settings:
        time_step = settings.get("time_step", time_step)

    core_species = reaction_model.core.species
    core_reactions = reaction_model.core.reactions
    edge_species = reaction_model.edge.species
    edge_reactions = reaction_model.edge.reactions

    species_index = {}
    for index, spec in enumerate(core_species):
        species_index[spec] = index

    simulator_settings = SimulatorSettings(atol=abs_tol, rtol=rel_tol)

    # Enable constant species for LiquidReactor
    if isinstance(reaction_system, LiquidReactor):
        if reaction_system.const_spc_names is not None:
            reaction_system.get_const_spc_indices(core_species)

    reaction_system.initialize_model(
        core_species,
        core_reactions,
        edge_species,
        edge_reactions,
        atol=simulator_settings.atol,
        rtol=simulator_settings.rtol,
        sens_atol=simulator_settings.sens_atol,
        sens_rtol=simulator_settings.sens_rtol,
        conditions=None,
    )

    # Copy the initial conditions to use in evaluating conversions
    y0 = reaction_system.y.copy()

    time = []
    core_species_concentrations = []
    core_reaction_rates = []

    next_time = initial_time
    step_time = initial_time
    terminated = False

    while not terminated:
        # Integrate forward in time to the next time point
        reaction_system.step(step_time)

        if reaction_system.t >= 0.9999 * next_time:
            next_time *= time_step
            time.append(reaction_system.t)
            core_species_concentrations.append(reaction_system.core_species_concentrations)
            core_reaction_rates.append(reaction_system.core_reaction_rates)

        # Finish simulation if any of the termination criteria are satisfied
        for term in reaction_system.termination:
            if isinstance(term, TerminationTime):
                if reaction_system.t > term.time.value_si:
                    terminated = True
                    break
            elif isinstance(term, TerminationConversion):
                index = species_index[term.species]
                if (y0[index] - reaction_system.y[index]) / y0[index] > term.conversion:
                    terminated = True
                    break

        # Increment destination step time if necessary
        if reaction_system.t >= 0.9999 * step_time:
            step_time *= 10.0

    time = np.array(time)
    core_species_concentrations = np.array(core_species_concentrations)
    core_reaction_rates = np.array(core_reaction_rates)

    return time, core_species_concentrations, core_reaction_rates


################################################################################


def load_chemkin_output(output_file, reaction_model):
    """
    Load the species concentrations from a Chemkin Output file in a simulation
    and generate the reaction rates at each time point.
    """
    import rmgpy.constants as constants
    from rmgpy.quantity import Quantity

    core_reactions = reaction_model.core.reactions
    species_list = reaction_model.core.species

    time = []
    core_species_concentrations = []
    core_reaction_rates = []

    with open(output_file, "r") as f:
        line = f.readline()
        while line != "" and "SPECIFIED END" not in line:
            line.strip()
            tokens = line.split()
            if " TIME " in line:
                # Time is in seconds
                time.append(float(tokens[-2]))
            elif " PRESSURE " in line:
                # Pressure from Chemkin is in atm
                P = Quantity(float(tokens[-2]), "atm")
            elif " TEMPERATURE " in line:
                # Temperature from Chemkin in in K
                T = Quantity(float(tokens[-2]), "K")
            elif " MOLE FRACTIONS " in line:
                # Species always come in the same order as listed in chem.inp
                molefractions = []
                line = f.readline()  # This one reads the blank line which follows
                line = f.readline()
                while line.strip() != "":
                    tokens = line.split()
                    for value in tokens[2::3]:
                        # Make all concentrations positive
                        if value.find("-") == 0:
                            value = value.replace("-", "", 1)
                        # Sometimes chemkin removes the `E` in scientific notation due to lack of space,
                        # rendering invalid float values.  If this is the case, add it in.
                        if value.find("-") != -1:
                            if value.find("E") == -1:
                                value = value.replace("-", "E-")

                        molefractions.append(float(value))

                    line = f.readline()

                total_concentration = P.value_si / constants.R / T.value_si
                core_species_concentrations.append([molefrac * total_concentration for molefrac in molefractions])
                core_rates = []
                for reaction in core_reactions:
                    rate = reaction.get_rate_coefficient(T.value_si, P.value_si)
                    for reactant in reaction.reactants:
                        rate *= molefractions[species_list.index(reactant)] * total_concentration
                    core_rates.append(rate)

                if core_rates:
                    core_reaction_rates.append(core_rates)

            line = f.readline()

    time = np.array(time)
    core_species_concentrations = np.array(core_species_concentrations)
    core_reaction_rates = np.array(core_reaction_rates)

    return time, core_species_concentrations, core_reaction_rates


################################################################################


def create_flux_diagram(
    input_file,
    chemkin_file,
    species_dict,
    save_path=None,
    species_path=None,
    java=False,
    settings=None,
    chemkin_output="",
    central_species_list=None,
    superimpose=False,
    save_states=False,
    read_states=False,
    diffusion_limited=True,
    check_duplicates=True,
):
    """
    Generates the flux diagram based on a condition 'input_file', chemkin.inp chemkin_file,
    a species_dict txt file, plus an optional chemkin_output file.
    """
    if species_path is None:
        species_path = os.path.join(os.path.dirname(input_file), "species")
        generate_images = True
    else:
        generate_images = False

    print("Loading RMG job...")
    rmg = load_rmg_job(input_file, chemkin_file, species_dict, generate_images=generate_images, use_java=java, check_duplicates=check_duplicates)

    if save_path is None:
        save_path = os.path.join(rmg.output_directory, "flux")

    # if you have a chemkin output, then you only have one reaction_system
    if chemkin_output:
        out_dir = os.path.join(save_path, "1")
        try:
            os.makedirs(out_dir)
        except OSError:
            pass

        print("Extracting species concentrations and calculating reaction rates from chemkin output...")
        time, core_species_concentrations, core_reaction_rates = load_chemkin_output(chemkin_output, rmg.reaction_model)

        print("Generating flux diagram for chemkin output...")
        generate_flux_diagram(
            rmg.reaction_model,
            time,
            core_species_concentrations,
            core_reaction_rates,
            out_dir,
            central_species_list=central_species_list,
            superimpose=superimpose,
            species_directory=species_path,
            settings=settings,
        )

    else:
        # Generate a flux diagram video for each reaction system
        for index, reaction_system in enumerate(rmg.reaction_systems):
            out_dir = os.path.join(save_path, "{0:d}".format(index + 1))
            try:
                os.makedirs(out_dir)
            except OSError:
                # Fail silently on any OS errors
                pass

            # If there is no termination time, then add one to prevent jobs from
            # running forever
            if not any([isinstance(term, TerminationTime) for term in reaction_system.termination]):
                reaction_system.termination.append(TerminationTime((1e10, "s")))

            states_file = os.path.join(out_dir, "states.npz")
            if read_states:
                print("Reading simulation states from file...")
                states = np.load(states_file)
                time = states["time"]
                core_species_concentrations = states["core_species_concentrations"]
                core_reaction_rates = states["core_reaction_rates"]
            else:
                # Enable diffusion-limited rates
                if diffusion_limited and isinstance(reaction_system, LiquidReactor):
                    rmg.load_database()
                    solvent_data = rmg.database.solvation.get_solvent_data(rmg.solvent)
                    diffusion_limiter.enable(solvent_data, rmg.database.solvation)

                print("Conducting simulation of reaction system {0:d}...".format(index + 1))
                time, core_species_concentrations, core_reaction_rates = simulate(rmg.reaction_model, reaction_system, settings)

                if save_states:
                    np.savez_compressed(
                        states_file, time=time, core_species_concentrations=core_species_concentrations, core_reaction_rates=core_reaction_rates
                    )

            print("Generating flux diagram for reaction system {0:d}...".format(index + 1))
            generate_flux_diagram(
                rmg.reaction_model,
                time,
                core_species_concentrations,
                core_reaction_rates,
                out_dir,
                central_species_list=central_species_list,
                superimpose=superimpose,
                species_directory=species_path,
                settings=settings,
            )
