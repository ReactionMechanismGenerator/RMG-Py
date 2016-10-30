
from cnn_model import build_model

class Predictor(object):

	def __init__(self):

		self.model = None

	def build_model(self):

		self.model = build_model()

	def train(self):

		pass

	def load_parameters(self):

		pass

	def save_parameters(self):

		pass

	def predict(self):

		pass

def readInputFile(path, predictor):

	pass
