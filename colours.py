# coding=utf-8
# Provides colouring of text
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021
class Colours:

	R = '\033[91m'
	G = '\033[92m'
	B = '\033[34m'
	RST = '\033[0m'

	@staticmethod
	def red(x):
		return Colours.R+x+Colours.RST

	@staticmethod
	def green(x):
		return Colours.G+x+Colours.RST

	@staticmethod
	def blue(x):
		return Colours.B+x+Colours.RST
