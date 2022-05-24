import matplotlib.pyplot as plt
from sys import stdin
import numpy as np
from drawing import Drawing

pos = -1

drawing = Drawing()
drawing.read_input_arrangement()
drawing.draw_arrangement()
drawing.read_guards_paths()
# drawing.draw_guards()
# drawing.draw_guards_paths()
# drawing.draw_visibility_regions(pos)
# drawing.draw_guards_dfs(pos)
drawing.draw_guard_visibility_dfs()

# plt.axis('off')

plt.show()