import matplotlib.pyplot as plt
from sys import stdin
import numpy as np
from drawing import Drawing

drawing = Drawing()
drawing.read_input_arrangement()
drawing.draw_arrangement()
drawing.read_guards_paths()
drawing.draw_guards_paths()

plt.axis('off')

plt.show()