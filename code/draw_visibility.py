from skgeom.draw import draw
from skgeom import Segment2, Point2, arrangement, RotationalSweepVisibility
import matplotlib.pyplot as plt


INPUT = "arrangement_main_out"


arr = arrangement.Arrangement()

"""The format of the input file is:
    E                     * number of edges
    p1.x p1.y p2.x p2.y   * edge with endpoints coordinates separated by spaces p1(x, y)p2(x, y)
    p3.x p3.y p4.x p4.y
    ...
"""
with open(INPUT, "r") as f:
    E = int(f.readline())

    for e in range(E):
        segment = list(map(int, f.readline().split()))
        p1 = Point2(*segment[:2])
        p2 = Point2(*segment[2:])

        arr.insert(Segment2(p1, p2))


vs = RotationalSweepVisibility(arr)

q = Point2(0.5, 2)
z = Point2(1, 1.5)
face1 = arr.find(q)
face2 = arr.find(z)
vx1 = vs.compute_visibility(q, face1)
vx2 = vs.compute_visibility(z, face2)


for v in vx1.halfedges:
    draw(v.curve(), point = q, color='red', visible_point=False, fill = True)
    # draw(v.curve(), z, color='green', visible_point=False)

for v in vx2.halfedges:
    draw(v.curve(), point = z, color='blue', visible_point=False, fill = True)
    # draw(v.curve(), z, color='green', visible_point=False)
for he in arr.halfedges:
    draw(he.curve(), visible_point=True)

draw(q,  color='magenta')
draw(z, color='magenta')


plt.show()

