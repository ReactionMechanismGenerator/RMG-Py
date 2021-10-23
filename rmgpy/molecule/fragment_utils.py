import numpy as np

def shuffle(conc, seed=None):
    """
    Randomly shuffle a list of fragments
    """
	idx_arr = np.arange(len(conc))

	if seed is not None:
		np.random.seed(seed)
	np.random.shuffle(idx_arr)

	return [conc[idx] for idx in idx_arr]

def grind(conc, size):
    """
    Split fragment concentrations into several repeating concentration units with specified size 
    """
	grinded_conc = []
	for label, c in conc:
		times = int(c/size)
		grinded_conc.extend([(label, size)]*times)

		if c-size*times > 0:
			grinded_conc.append((label, c-size*times))

	return grinded_conc

def flatten(combo):
	"""
	Given a combo nested `tuple`, e.g., 
	((('LY', 'XR'), ('LWL', 'RUR'))
	return a list of labels contained in 
	the combo ['LY', 'XR', 'LWL', 'RUR']
	"""
	return_list = []
	for i in combo:
		if isinstance(i, tuple):
			return_list.extend(flatten(i))
		else:
			return_list.append(i)
	return return_list
