from skgeom.draw import draw
from skgeom import Segment2, Point2, arrangement, RotationalSweepVisibility, TriangularExpansionVisibility
import matplotlib.pyplot as plt
from numpy import random


INPUT = "pentagram.out"


arrangement = arrangement.Arrangement()
vs = TriangularExpansionVisibility(arrangement)
guards = []

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
with open(INPUT, 'r') as f:
    E = int(input())

    for e in range(E):
        segment = list(map(int, input().split()))
        p1 = Point2(*segment[:2])
        p2 = Point2(*segment[2:])

        arrangement.insert(Segment2(p1, p2))
    
    IV = int(input())

    for iv in range(IV):
        vertex = list(map(float, input().split()))
        q = Point2(*vertex)
        guards.append(q)

# draw the guards and their visibility regions
for guard in guards:
    face = arrangement.find(guard)
    vx = vs.compute_visibility(guard, face)

    # use a random colour for each guard
    color = random.rand(3,)
    for v in vx.halfedges:
        draw(v.curve(), point = guard, color=color, visible_point=False, fill = True)
    
    draw(guard, color='magenta')

# draw the polygon's boundaries
for he in arrangement.halfedges:
    draw(he.curve(), visible_point=True)

plt.show()

