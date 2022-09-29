import matplotlib.pyplot as plt
from sys import stdin
import numpy as np
from drawing import Drawing

# drawing = Drawing()
# drawing.read_input_arrangement()
# drawing.draw_arrangement()
# drawing.read_iteration()
# drawing.draw_guards()
# drawing.draw_guards_paths()
# drawing.draw_visibility_regions(pos)
# drawing.draw_guards_dfs(pos)
# drawing.draw_guard_visibility_dfs()
# drawing.plot_area_time()

for comb_number in ['_2', '_3', '', '5', '6', '7', '8', '9', '10']:
    drawing_comb = Drawing(file = f'../results/comb{comb_number}.out')
    drawing_comb.plot_area_combs()
# plt.axis('off')
plt.show()
