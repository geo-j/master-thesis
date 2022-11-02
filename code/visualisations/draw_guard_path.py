import matplotlib.pyplot as plt
from sys import stdin
import numpy as np
from drawing import Drawing
import time
import os

PATH = 'results/experiments/'
shape = 'random'
DATE = time.strftime("%Y-%m-%d")
if not os.path.exists(PATH + DATE):
    os.makedirs(PATH + DATE)

def compute_average(heuristic: str) -> np.array:
    average_areas = []
    max_iterations = None
    for polygon_number in range(1, 21):
        polygon = Drawing(file = f'{PATH + heuristic}/{shape}_{polygon_number}.out')
        polygon.read_input_arrangement()
        polygon.read_iteration()
        
        average_areas.append([area * 100 / float(polygon.max_area) for area in polygon.areas])

        if max_iterations is None or polygon.n_iterations > max_iterations:
            max_iterations = polygon.n_iterations

    for i in range(len(average_areas)):
        average_areas[i] = average_areas[i] + [100.0] * (max_iterations + 1 - len(average_areas[i]))

    return np.array(average_areas)

def plot_average_seen_area(all_heuristics: np.array, no_pull: np.array) -> None:
    plt.plot(np.average(all_heuristics, axis = 0))
    plt.plot(np.average(no_pull, axis = 0))
    plt.xlabel('# iterations')
    plt.ylabel('% total seen area')
    plt.xlim(0, len(no_pull[-1]))
    plt.grid()
    plt.legend(['all heuristics', 'no pull'])
    plt.title(f'Average Total Seen Area for the {shape} Polygon')
    plt.savefig(f'{PATH + DATE}/{time.strftime("%H%M")}_avg_area_{shape}.png', format = 'png', dpi = 300, bbox_inches = 'tight')
    plt.show()

def plot_average_n_iterations(all_heuristics: np.array, no_pull: np.array) -> None:
    plt.boxplot([[np.count_nonzero(areas != 100) for areas in all_heuristics], [np.count_nonzero(areas != 100) for areas in no_pull]])
    plt.ylabel('# iterations')
    plt.title(f'Average Number of Iterations for the {shape} Polygon')
    plt.grid()
    plt.xticks([1, 2], ['all heuristics', 'no pull'])
    plt.savefig(f'{PATH + DATE}/{time.strftime("%H%M")}_avg_iterations_{shape}.png', format = 'png', dpi = 300, bbox_inches = 'tight')
    plt.show()

all_heuristics = compute_average('all_heuristics')
no_pull = compute_average('no_pull')
# print(np.average(all_heuristics, axis = 0), np.average(no_pull, axis = 0))
# plt.plot(np.average(all_heuristics, axis = 0))
# plt.show()
plot_average_seen_area(all_heuristics, no_pull)
plot_average_n_iterations(all_heuristics, no_pull)


# drawing = Drawing()
# drawing.read_input_arrangement()
# drawing.draw_arrangement()
# drawing.read_iteration()
# drawing.draw_guard_visibility_dfs()
# drawing.plot_area_time()


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
