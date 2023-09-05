import numpy as np
from Gridify import Mesh
import matplotlib.pyplot as plt 

M = Mesh(5,5)
#print(M.element_list)
#M.update_element(0.5,0.5,1)
"""print(f"The shape function is evaluate as {M.N}")
print(f"The derivative of shape function is evaluate as {M.dN}")
print(f"The Jacobian is evaluate as {M.J}")
print(M.element_list[1])"""
for key, value in M.element_list.items():
    print(f"{key}: {value}")
    print(f"node {value[0]} = {M.nl[value[0]]}")
    print(f"node {value[1]} = {M.nl[value[1]]}")
    print(f"node {value[2]} = {M.nl[value[2]]}")
    print(f"node {value[3]} = {M.nl[value[3]]}")
    print("\n")

print(M.nl)
x_cord = [points[0] for points in M.nl]
y_cord = [points[1] for points in M.nl]
print(f"The size of the node list is {np.shape(M.nl)}")
#print(M.boundary_nodes)
plt.scatter(x_cord, y_cord, marker='o')
for i, point in enumerate(M.nl):
    plt.annotate(str(i), (point[0], point[1]), textcoords="offset points", xytext=(0,10), ha='center')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.grid(True)
plt.title('Meshed Geometry')
plt.show()
