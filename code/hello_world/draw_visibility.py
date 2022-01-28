from skgeom.draw import draw
from skgeom import Segment2, Point2, arrangement, RotationalSweepVisibility
import matplotlib.pyplot as plt

arr = arrangement.Arrangement()

with open("arrangement", "r") as f:
    V = int(f.readline())
    # print(V)

    for v in range(V):
        point = list(map(int, f.readline().split()))
        # x = point[0]
        # y = point[1]
        p = Point2(*point)

    E = int(f.readline())

    for e in range(E):
        segment = list(map(int, f.readline().split()))
        p1 = Point2(*segment[:2])
        p2 = Point2(*segment[2:])

        arr.insert(Segment2(p1, p2))
        # arr.insert(Segment2(p2, p1))

arr.insert(Segment2(Point2(1, 3), Point2(2, 2)))
arr.insert(Segment2(Point2(3, 3), Point2(4, 4)))
arr.insert(Segment2(Point2(3, 1), Point2(4, 2)))

# arr.insert(Segment2(Point2(2, 2), Point2(1, 3)))

# draw(arr)
# plt.show()

# M = 10
# outer = [
#     Segment2(Point2(-M, -M), Point2(-M, M)), Segment2(Point2(-M, M), Point2(M, M)),
#     Segment2(Point2(M, M), Point2(M, -M)), Segment2(Point2(M, -M), Point2(-M, -M))
# ]

# segments = [
#     Segment2(Point2(0, 0), Point2(0, 4)), Segment2(Point2(2, 4), Point2(8, 4)),
#     Segment2(Point2(3, 4), Point2(-8, 0)), Segment2(Point2(1, 0), Point2(0, 5)),
# ]

# for s in outer:
#     arr.insert(s)

# for s in segments:
#     arr.insert(s)

vs = RotationalSweepVisibility(arr)

q = Point2(1.5, 3)
face = arr.find(q)
vx = vs.compute_visibility(q, face)

for he in arr.halfedges:
    draw(he.curve(), visible_point=True)
for v in vx.halfedges:
    draw(v.curve(), color='red', visible_point=False)

draw(q, color='magenta')

plt.show()

