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
z = Point2(3, 1.5)
face1 = arr.find(q)
face2 = arr.find(z)
vx1 = vs.compute_visibility(q, face1)
vx2 = vs.compute_visibility(z, face2)


for v in vx1.halfedges:
    draw(v.curve(), point = q, color='red',   visible_point=False, fill = True)
    # draw(v.curve(), z, color='green', visible_point=False)

for v in vx2.halfedges:
    draw(v.curve(), point = z, color='yellow',  visible_point=False, fill = True)
    # draw(v.curve(), z, color='green', visible_point=False)
for he in arr.halfedges:
    draw(he.curve(), visible_point=True)

draw(q,  color='magenta')
draw(z, color='magenta')


plt.show()

