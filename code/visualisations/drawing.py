from skgeom.draw import draw
from skgeom import Segment2, Point2, arrangement, RotationalSweepVisibility, TriangularExpansionVisibility
import matplotlib.pyplot as plt
from numpy import random

class Drawing(object):
    def __init__(self) -> None:
        self.arrangement = arrangement.Arrangement()
        self.vs = TriangularExpansionVisibility(self.arrangement)
        self.guards = []
    
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
            segment = list(map(int, input().split()))
            p1 = Point2(*segment[:2])
            p2 = Point2(*segment[2:])

            self.arrangement.insert(Segment2(p1, p2))

    def read_input_guards(self) -> None:
        IV = int(input())

        for _ in range(IV):
            vertex = list(map(float, input().split()))
            q = Point2(*vertex)
            self.guards.append(q)
    
    def read_all_input(self) -> None:
        self.read_input_arrangement()
        self.read_input_guards()
    
    ###########
    #  Output #
    ###########
    def draw_visibility_regions(self) -> None:
        # draw the guards and their visibility regions
        for guard in self.guards:
            face = self.arrangement.find(guard)
            vx = self.vs.compute_visibility(guard, face)

            # use a random colour for each guard
            color = random.rand(3, )
            for v in vx.halfedges:
                draw(v.curve(), point = guard, color = color, visible_point = False, fill = True)
    
    def draw_guards(self) -> None:
        for guard in self.guards:
            draw(guard, color='magenta')

    def draw_arrangement(self) -> None:
        # draw the polygon's boundaries
        for he in self.arrangement.halfedges:
            draw(he.curve(), visible_point = True)

    def draw_all(self) -> None:
        self.draw_visibility_regions()
        self.draw_guards()
        self.draw_arrangement()