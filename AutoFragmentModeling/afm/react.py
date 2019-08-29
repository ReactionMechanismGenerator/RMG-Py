

def react_fragments(kinetics_db, 
                    fragment_tuple,
                    products=None, 
                    only_families=None,
                    prod_resonance=True):

    """
    Given a tuple of fragment objects, generates all possible reactions
    from the loaded reaction families and combines degenerate reactions.

    The generated reactions are deflated.
    """

    reactions = kinetics_db.react_molecules(fragment_tuple, 
                                            products=products, 
                                            only_families=only_families, 
                                            prod_resonance=prod_resonance)

    return reactions
