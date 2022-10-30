import matplotlib.pyplot as plt
from sys import stdin
import numpy as np
from drawing import Drawing
import time
import os

PATH = 'results/'
DATE = time.strftime("%Y-%m-%d")
if not os.path.exists(PATH + DATE):
    os.makedirs(PATH + DATE)
drawing = Drawing()
drawing.read_input_arrangement()
drawing.draw_arrangement()
drawing.read_iteration()
drawing.draw_guard_visibility_dfs()
drawing.plot_area_time()

# for comb_number in ['_2', '_3', '', '5', '6', '7', '8', '9', '10', '15', '20']:
# # for comb_number in ['_2', '_3', '', '6']:
# for polygon in ['_broken', '_good']:
#     drawing_comb = Drawing(file = f'results/love{polygon}.out')
#     drawing_comb.read_input_arrangement()
#     drawing_comb.read_iteration()
#     drawing_comb.plot_area_multiple()

# # # plt.legend(['2 teeth', '3 teeth', '4 teeth', '6 teeth'])
# # plt.legend(['2 teeth', '3 teeth', '4 teeth', '5 teeth', '6 teeth', '7 teeth', '8 teeth', '9 teeth', '10 teeth', '15 teeth', '20 teeth'])
# plt.legend(['arbitrary positions', 'irrational positions'])
# plt.title('Total Areas Seen')
# plt.xlabel('# iterations')
# plt.ylabel('total area seen (%)')
# # plt.xlim(0, 90)
# # plt.xticks(np.arange(0, 19))
# plt.grid()
# # plt.axis('off')
# plt.savefig(f'{PATH + DATE}/{time.strftime("%H%M")}_area.png', format = 'png', dpi = 300, bbox_inches = 'tight')
# # plt.show()
