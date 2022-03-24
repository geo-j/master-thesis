import matplotlib.pyplot as plt
from sys import stdin
import numpy as np
from drawing import Drawing

drawing = Drawing()
drawing.read_input_arrangement()
drawing.draw_arrangement()

xs, ys = [], []

for line in stdin:
    x, y = map(float, line.strip().split())

    xs.append(x)
    ys.append(y)

plt.scatter(xs, ys, color = 'magenta')

u = np.diff(xs)
v = np.diff(ys)
pos_x = xs[:-1] + u / 2
pos_y = ys[:-1] + v / 2
norm = np.sqrt(u ** 2 + v ** 2) 

# TODO: scale still needs to be figured out
plt.quiver(pos_x, pos_y, u / norm, v / norm, angles = 'xy', pivot = 'mid', width = 0.005, scale = 50)

plt.show()