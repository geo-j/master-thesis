from skgeom.draw import draw
from skgeom import Segment2, Point2, arrangement, RotationalSweepVisibility, TriangularExpansionVisibility, intersection, Vector2, Line2
import matplotlib.pyplot as plt
from numpy import random, diff, sqrt
from sys import stdin
from collections import defaultdict
import time
import os
from numpy import arange

PATH = 'results/'
DATE = time.strftime("%Y-%m-%d")
if not os.path.exists(PATH + DATE):
    os.makedirs(PATH + DATE)

def distance(p1: Point2, p2: Point2) -> float:
    return (p1.x() - p2.x()) * (p1.x() - p2.x()) + (p1.y() - p2.y()) * (p1.y() - p2.y())

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
        self.hs_x = defaultdict(lambda: defaultdict(list))
        self.hs_y = defaultdict(lambda: defaultdict(list))
        self.local_areas = defaultdict(list)
        self.areas = []
        self.max_area = None

    
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
        p1 = None
        p2 = None

        for _ in range(E):
            line_segment = input().split()
            x = None
            y = None
            for i in range(len(line_segment)):
                if '/' in line_segment[i]:
                    fraction = list(map(float, line_segment[i].split('/')))
                    fraction_result = fraction[0] / fraction[1]
                    if i % 2 == 0:
                        x = fraction_result
                    else:
                        y = fraction_result
                else:
                    if i % 2 == 0:
                        x = float(line_segment[i])
                    else:
                        y = float(line_segment[i])
                if x is not None and y is not None:
                    if i == 1:
                        p1 = Point2(x, y)
                    else:
                        p2 = Point2(x, y)
                



            # segment = list(map(float, input().split()))
            # p1 = Point2(*segment[:2])
            # p2 = Point2(*segment[2:])

            self.arrangement.insert(Segment2(p1, p2))
            self.halfedges.append(Segment2(p1, p2))

    def read_input_guards(self) -> None:
        IV = int(input())

        for i in range(IV):
            vertex = list(map(float, input().split()))
            q = Point2(*vertex)
            # self.guards[f'g{i}'].append(q)
            self.guards.append(q)
    
    def read_iteration(self) -> None:
        iteration = None
        i = None

        for line in stdin:
            # if iteration is not None and iteration < 10:
            #     print(iteration, line)
            if line.startswith('total'):
                self.max_area = float(line[11:].strip())
                # print(self.max_area)
            elif line.startswith('i='):   # get current iteration
                iteration = int(line[2:].strip())
            elif line.startswith('g'):  # get current guard coords
                i = int(line[1])    # get the guard index
                x, y = map(float, line[3:].strip().split())     # get the coords after removing the guard info

                self.xs[f'g{i}'].append(x)
                self.ys[f'g{i}'].append(y)
            elif line.startswith('D'):  # get current guard gradient info
                # i = int(line[2])    # get guard index
                x, y = map(float, line[3:].strip().split())     # get the coords after removing the guard info
                self.dfs_x[f'g{i}'][iteration].append(x)
                self.dfs_y[f'g{i}'][iteration].append(y)
            elif line.startswith('h'):  # get current guard gradient info
                # i = int(line[2])    # get guard index
                x, y = map(float, line[2:].strip().split())     # get the coords after removing the guard info
                self.hs_x[f'g{i}'][iteration].append(x)
                self.hs_y[f'g{i}'][iteration].append(y)
            elif line.startswith('area='):  # get current total seen area
                area = float(line[5:].strip())
                self.areas.append(area)
            elif line.startswith('area'):   # get current guard's area
                # i = int(line[4])
                area = float(line[6:].strip())
                self.local_areas[f'g{i}'].append(area)
        
        self.n_iterations = iteration

    
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
            # print(guard, pos, self.xs[guard])
            if len(self.xs[guard]) <= pos:
                pos = -1
            g = Point2(self.xs[guard][pos], self.ys[guard][pos])
            face = self.arrangement.find(g)

            if (type(face) is arrangement.Face and face.is_unbounded()) or (type(face) is arrangement.Halfedge and face.face().is_unbounded()) or (type(face) is arrangement.Vertex):
                for half_edge in self.arrangement.halfedges:
                    segment = Segment2(half_edge.source().point(), half_edge.target().point())
                    # print(segment)
                    # print(f'trying to place {g} on {segment}')

                    if half_edge.target().point() == guard and not half_edge.face().is_unbounded():
                        face = half_edge
                        break
                else:
                    # if the point is sliiightly outside the boundary due to rounding
                    min_edge = None
                    min_p = None
                    for half_edge in self.arrangement.halfedges:
                        # take inner half-edges
                        if not half_edge.face().is_unbounded():
                            segment = Segment2(half_edge.source().point(), half_edge.target().point())
                            line = Line2(half_edge.source().point(), half_edge.target().point())
                            # project the guard on the boundary segment
                            p = line.projection(g)
                            # print(f'looking at {segment}')

                            # init min with the first half-edge found; if it's not a good one, it will be overwritten anyway
                            if min_edge is None:
                                min_edge = half_edge
                                min_p = p

                            # if the projection of the guard is still inside the polygon
                            if segment.collinear_has_on(p):
                                if distance(p, g) <= distance(min_p, g):
                                    # the projection has to be either on the target edge of the segment, or anywhere in between
                                    if half_edge.target().point() == g or (half_edge.target().point() != g and half_edge.source().point() != g):
                                        min_edge = half_edge
                                        min_p = p
                                        # print(f'\tplaced {g} on {p}')
                            else:
                                # if the projection is outside of the polygon, place it on the closest half-edge target edge
                                if distance(half_edge.target().point(), g) <= distance(min_p, g):
                                    min_edge = half_edge
                                    min_p = half_edge.target().point()
                                    # print(f'\tplaced {g} on {p}')

                    face = min_edge
                    g = min_p
                    self.xs[guard][pos] = float(g.x())
                    self.ys[guard][pos] = float(g.y())

            
            
            # try:
            vx = self.vs.compute_visibility(g, face)
            # use a random colour for each guard
            # color = random.rand(3, )
            color = color_list[i]
            for v in vx.halfedges:
                draw(v.curve(), point = g, visible_point = False, fill = True, color = color)
            # except:
            #     pass

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
            if len(self.xs[guard]) <= pos:
                pos = -1
            plt.scatter(self.xs[guard][pos], self.ys[guard][pos], color = color)

            plt.quiver([self.xs[guard][pos]] * len(self.dfs_x[guard][pos]), [self.ys[guard][pos]] * len(self.dfs_x[guard][pos]), [self.dfs_x[guard][pos]], [self.dfs_y[guard][pos]], scale = 1, scale_units = 'xy', angles = 'xy', width = 0.0055, color = ['g'] * (len(self.dfs_x[guard][pos]) - 1) + ['r'])

            plt.quiver([self.xs[guard][pos]] * len(self.hs_x[guard][pos]), [self.ys[guard][pos]] * len(self.hs_x[guard][pos]), [self.hs_x[guard][pos]], [self.hs_y[guard][pos]], scale = 1, scale_units = 'xy', angles = 'xy', width = 0.0055, color = ['b'] * (len(self.hs_x[guard][pos]) - 1) + ['purple'])



    def draw_guard_visibility_dfs(self) -> None:
        self.draw_arrangement()
        for pos in range(self.n_iterations + 1):
            self.draw_visibility_regions(pos)
            self.draw_guards_dfs(pos)
            plt.title(f'Gradient Computation for Iteration #{pos}')
            plt.axis('off')

            plt.savefig(f'{PATH + DATE}/{time.strftime("%H%M")}_pos{pos}.png', format = 'png', dpi = 300)
            plt.clf()
            # plt.show()
            self.draw_arrangement()
        plt.clf()

    def plot_area_time(self) -> None:
        plt.plot([x * 100 / float(self.max_area) for x in self.areas])
        plt.xlim(0, len(self.areas))
        # plt.margins(x = 0)
        plt.xticks(arange(0, len(self.areas)))
        # plt.axhline(y = self.max_area, color = 'r', linestyle = '--')
        # plt.ylim(min(self.areas) * 100 / float(self.max_area) - 1, 101)

        for guard in self.xs.keys():
            # print([x for x in self.local_areas[guard]])
            plt.plot([x * 100 / float(self.max_area) for x in self.local_areas[guard]])

        plt.xlabel('# iterations')
        plt.ylabel('total area seen (%)')
        plt.title('Total Area Seen')
        plt.grid()
        plt.legend(['total'] + [x for x in self.local_areas.keys()])
        plt.savefig(f'{PATH + DATE}/{time.strftime("%H%M")}_area.png', format = 'png', dpi = 300, bbox_inches = 'tight')
        # plt.show()

    def draw_all(self) -> None:
        self.draw_visibility_regions()
        self.draw_guards()
        self.draw_arrangement()
        self.draw_guards_paths()