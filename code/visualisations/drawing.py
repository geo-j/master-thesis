from skgeom.draw import draw
from skgeom import Segment2, Point2, arrangement, RotationalSweepVisibility, TriangularExpansionVisibility, intersection
import matplotlib.pyplot as plt
from numpy import random, diff, sqrt
from sys import stdin
from collections import defaultdict

class Drawing(object):
    def __init__(self) -> None:
        self.arrangement = arrangement.Arrangement()
        self.vs = TriangularExpansionVisibility(self.arrangement)
        self.guards = []
        self.halfedges = []
        # dict for each guard's x and y coord
        self.xs = defaultdict(list)
        self.ys = defaultdict(list)
        self.dfs_x = defaultdict(lambda: defaultdict(list))
        self.dfs_y = defaultdict(lambda: defaultdict(list))

    
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
        for line in stdin:
            if line.startswith('g'):
                i = int(line[1])    # get the guard index
                x, y = map(float, line[3:].strip().split())     # get the coords after removing the guard info

                self.xs[f'g{i}'].append(x)
                self.ys[f'g{i}'].append(y)
            else:
                i = int(line[2])    # get guard index
                x, y = map(float, line[4:].strip().split())     # get the coords after removing the guard info
                # print(i, self.xs[f'g{i}'])
                self.dfs_x[f'g{i}'][(len(self.xs[f'g{i}']) - 1)].append(x)
                self.dfs_y[f'g{i}'][(len(self.xs[f'g{i}']) - 1)].append(y)

    
    def read_all_input(self) -> None:
        self.read_input_arrangement()
        self.read_input_guards()
        self.read_guards_paths()
    
    ###########
    #  Output #
    ###########
    def draw_visibility_regions(self, pos: int = -1) -> None:
        # draw the guards and their visibility regions
        for guard in self.xs.keys():
            g = Point2(self.xs[guard][pos], self.ys[guard][pos])
            # print(g)
            face = self.arrangement.find(g)
            # print(face)
            if type(face) is arrangement.Face and face.is_unbounded():
                # for half_edge in self.halfedges:
                for half_edge in self.arrangement.halfedges:
                    segment = Segment2(half_edge.source().point(), half_edge.target().point())
                    if segment.collinear_has_on(g):
                        face = half_edge
            
            if type(face) is arrangement.Face:
                        vx = self.vs.compute_visibility(g, face)
            else:
                # print(face.face().is_unbounded())
                if not face.face().is_unbounded():
                    try:
                        vx = self.vs.compute_visibility(g, face)
                    except:
                        g = Point2((face.source().point().x() + face.target().point().x()) / 2, (face.source().point().y() + face.target().point().y()) / 2)
                        vx = self.vs.compute_visibility(g, face)
                else:
                    vx = self.vs.compute_visibility(g, face.twin())

            # use a random colour for each guard
            color = random.rand(3, )
            for v in vx.halfedges:
                draw(v.curve(), point = g, color = color, visible_point = False, fill = True)
    
    def draw_guards(self) -> None:
        for guard in self.xs.keys():
            g = Point2(self.xs[guard][-1], self.ys[guard][-1])
            draw(g, color='magenta')

    def draw_arrangement(self) -> None:
        # draw the polygon's boundaries
        for he in self.arrangement.halfedges:
            draw(he.curve(), visible_point = True)

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
        for guard in self.xs.keys():
            plt.scatter(self.xs[guard][pos], self.ys[guard][pos])
            # print(f'guard: {guard}')
            # print(f'guard {guard}: {self.dfs_x[guard], self.dfs_y[guard]}')
            # print(self.dfs_x[guard][pos])
            # print([self.xs[guard][pos]] * len(self.dfs_x[guard][pos]))
            plt.quiver([self.xs[guard][pos]] * len(self.dfs_x[guard][pos]), [self.ys[guard][pos]] * len(self.dfs_x[guard][pos]), [self.dfs_x[guard][pos]], [self.dfs_y[guard][pos]], scale = 2, scale_units = 'inches', width = 0.005, color = ['g'] * (len(self.dfs_x[guard][pos]) - 1) + ['r'])


    def draw_all(self) -> None:
        self.draw_visibility_regions()
        self.draw_guards()
        self.draw_arrangement()
        self.draw_guards_paths()