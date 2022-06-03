from skgeom.draw import draw
from skgeom import Segment2, Point2, arrangement, RotationalSweepVisibility, TriangularExpansionVisibility, intersection, Vector2
import matplotlib.pyplot as plt
from numpy import random, diff, sqrt
from sys import stdin
from collections import defaultdict
import time
import os

PATH = 'results/'
DATE = time.strftime("%Y-%m-%d")
if not os.path.exists(PATH + DATE):
    os.makedirs(PATH + DATE)

class Drawing(object):
    def __init__(self) -> None:
        self.arrangement = arrangement.Arrangement()
        self.vs = TriangularExpansionVisibility(self.arrangement)
        self.guards = []
        self.halfedges = []
        # dict for each guard's x and y coord
        self.n_iterations = None
        self.xs = defaultdict(list)
        self.ys = defaultdict(list)
        self.dfs_x = defaultdict(lambda: defaultdict(list))
        self.dfs_y = defaultdict(lambda: defaultdict(list))
        self.local_areas = defaultdict(lambda: defaultdict(list))
        self.areas = []

    
    ###########
    #  Input #
    ###########
    """The format of the input file is:
        E                     * number of edges
        p1.x p1.y p2.x p2.y   * edge with endpoints coordinates separated by spaces p1(x, y)p2(x, y)
        p3.x p3.y p4.x p4.y
        ...
        IV                    * number of isolated vertices (intended to be guards)
        q1.x q1.y             * isolated vertex with coordinates q1(x, y) separated by spaces
        q2.x q2.y
        ...
    """
    def read_input_arrangement(self) -> None:
        E = int(input())

        for _ in range(E):
            segment = list(map(float, input().split()))
            p1 = Point2(*segment[:2])
            p2 = Point2(*segment[2:])

            self.arrangement.insert(Segment2(p1, p2))
            self.halfedges.append(Segment2(p1, p2))

    def read_input_guards(self) -> None:
        IV = int(input())

        for i in range(IV):
            vertex = list(map(float, input().split()))
            q = Point2(*vertex)
            # self.guards[f'g{i}'].append(q)
            self.guards.append(q)
    
    def read_guards_paths(self) -> None:
        iteration = None

        for line in stdin:
            if line.startswith('i='):   # get current iteration
                iteration = int(line[2:].strip())
            elif line.startswith('g'):  # get current guard coords
                i = int(line[1])    # get the guard index
                x, y = map(float, line[3:].strip().split())     # get the coords after removing the guard info

                self.xs[f'g{i}'].append(x)
                self.ys[f'g{i}'].append(y)
            elif line.startswith('D'):  # get current guard gradient info
                i = int(line[2])    # get guard index
                x, y = map(float, line[4:].strip().split())     # get the coords after removing the guard info
                self.dfs_x[f'g{i}'][iteration].append(x)
                self.dfs_y[f'g{i}'][iteration].append(y)
            elif line.startswith('area='):  # get current total seen area
                area = float(line[5:].strip())
                self.areas.append(area)
            elif line.startswith('area'):   # get current guard's area
                i = int(line[4])
                area = float(line[6:].strip())
                self.local_areas[f'g{i}'][iteration].append(area)
        
        self.n_iterations = iteration - 1

    
    def read_all_input(self) -> None:
        self.read_input_arrangement()
        self.read_input_guards()
        self.read_guards_paths()
    
    ###########
    #  Output #
    ###########
    def draw_visibility_regions(self, pos: int = -1) -> None:
        # draw the guards and their visibility regions
        color_list = plt.rcParams['axes.prop_cycle'].by_key()['color']
        i = 0
        for guard in self.xs.keys():
            # print(guard, pos)
            g = Point2(self.xs[guard][pos], self.ys[guard][pos])
            face = self.arrangement.find(g)

            # print(face)
            if type(face) is arrangement.Face and face.is_unbounded():
                for half_edge in self.arrangement.halfedges:
                    segment = Segment2(half_edge.source().point(), half_edge.target().point())
                    if segment.collinear_has_on(g):
                        face = half_edge
            
            if type(face) is arrangement.Halfedge and face.face().is_unbounded():
                face = face.twin()
            
            try:
                vx = self.vs.compute_visibility(g, face)
                # use a random colour for each guard
                # color = random.rand(3, )
                color = color_list[i]
                for v in vx.halfedges:
                    draw(v.curve(), point = g, visible_point = False, fill = True, color = color)
            except:
                pass

            i += 1


    
    def draw_guards(self) -> None:
        for guard in self.xs.keys():
            g = Point2(self.xs[guard][-1], self.ys[guard][-1])
            draw(g, color='magenta')

    def draw_arrangement(self) -> None:
        # draw the polygon's boundaries
        for he in self.arrangement.halfedges:
            draw(he.curve(), visible_point = True, color = 'gray')

    def draw_guards_paths(self) -> None:
        for guard in self.xs.keys():
            plt.scatter(self.xs[guard], self.ys[guard])

            u = diff(self.xs[guard])
            v = diff(self.ys[guard])
            pos_x = self.xs[guard][:-1] + u / 2
            pos_y = self.ys[guard][:-1] + v / 2
            norm = sqrt(u ** 2 + v ** 2) 

            plt.quiver(pos_x, pos_y, u / norm, v / norm, angles = 'xy', pivot = 'mid', width = 0.005, scale = 5, scale_units = 'inches')

    def draw_guards_dfs(self, pos: int = -1) -> None:
        color_list = plt.rcParams['axes.prop_cycle'].by_key()['color']
        i = 0

        for guard in self.xs.keys():
            color = color_list[i]
            i += 1
            plt.scatter(self.xs[guard][pos], self.ys[guard][pos], color = color)

            # print(f'guard: {guard}')
            # print(f'guard {guard}: {self.dfs_x[guard], self.dfs_y[guard]}')
            # print(self.dfs_x[guard][pos])
            # print([self.xs[guard][pos]] * len(self.dfs_x[guard][pos]))
            plt.quiver([self.xs[guard][pos]] * len(self.dfs_x[guard][pos]), [self.ys[guard][pos]] * len(self.dfs_x[guard][pos]), [self.dfs_x[guard][pos]], [self.dfs_y[guard][pos]], scale = 1, scale_units = 'xy', angles = 'xy', width = 0.0055, color = ['g'] * (len(self.dfs_x[guard][pos]) - 1) + ['r'])


    def draw_guard_visibility_dfs(self) -> None:
        self.draw_arrangement()
        for pos in range(self.n_iterations + 1):
            self.draw_visibility_regions(pos)
            self.draw_guards_dfs(pos)
            plt.title(f'Gradient Computation for Iteration #{pos}')
            # plt.axis('off')

            plt.savefig(f'{PATH + DATE}/{time.strftime("%H%M")}_pos{pos}.png', format = 'png')
            plt.clf()
            # plt.show()
            self.draw_arrangement()
        plt.clf()

    def plot_area_time(self) -> None:
        plt.plot(self.areas)
        plt.axhline(y = max(self.areas), color = 'r', linestyle = '--')
        plt.ylim(0, max(self.areas) + 1)
        plt.xlabel('# iterations')
        plt.ylabel('total area seen')
        plt.title('Total Area Seen')
        plt.savefig(f'{PATH + DATE}/{time.strftime("%H%M")}_area.png', format = 'png')
        # plt.show()

    def draw_all(self) -> None:
        self.draw_visibility_regions()
        self.draw_guards()
        self.draw_arrangement()
        self.draw_guards_paths()