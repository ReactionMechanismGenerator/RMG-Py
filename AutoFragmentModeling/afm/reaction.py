
class FragmentReaction(object):

	def __init__(self,
				index=-1,
				reactants=None,
				products=None,
				kinetics=None,
				reversible=False,
				pairs=None,
				family=None,
				reaction_repr=None
				):

		
		self.index = index
		self.reactants = reactants
		self.products = products
		self.kinetics = kinetics
		self.reversible = reversible
		self.pairs = pairs
		self.family = family
		self.reaction_repr = reaction_repr
	
	def __str__(self):
		"""
		Return a string representation of the reaction, in the form 'A + B <=> C + D'.
		"""
		arrow = ' <=> '
		if not self.reversible: arrow = ' => '
		return arrow.join([' + '.join([str(s) for s in self.reactants]), ' + '.join([str(s) for s in self.products])])
