"""This file converts the counterclockwise vertex-only input from the https://www.ic.unicamp.br/~cid/Problem-instances/Art-Gallery/AGPVG/index.html library into my own way of reading input.

Each file consists of one line divided in two parts. The first part is an integer value that represents the number of vertices of the polygon.
The second part is a counterclockwise sequence of the vertices. Each vertex is represented by its x and y coordinates each of which is written as the quotient of two integers int/int.

As an example, here is the representation of a square with coordinates (1, 1) (50, 1) (50, 50) and (1, 50):
        4   1/1   1/1   100/2   1/1   500/10   50/1   1/1   100/2
        

the input is converted into 
E                     number of edges
p1.x p1.y p2.x p2.y   edge with endpoints coordinates separated by spaces p1(x, y)p2(x, y)
p3.x p3.y p4.x p4.y
...

The previous input would be then converted to
4
1/1 100/2 500/10 50/1
500/10 50/1 100/2 1/1
100/2 1/1 1/1 1/1
1/1 1/1 1/1 100/2"""

# read input line
line = input().strip().split()

# get number of vertices
V = int(line[0])
print(V)

# read and remember first vertex to add finishing segment at the end
x1 = line[len(line) - 2]
y1 = line[len(line) - 1]
first_x = x1
first_y = y1

# read vertices clockwise and print the segment between the current vertex and previous one
for v in range(2 * V - 2, 1, -2):
    x2 = line[v - 1]
    y2 = line[v]

    print(x1, y1, x2, y2)

    x1 = x2
    y1 = y2

print(x1, y1, first_x, first_y)