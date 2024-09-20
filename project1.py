import numpy as np
import json

# 1 TODO: Write a system of equations Ax = b to solve for the node voltages
# DONE (More Details of the generic compositions of A & b in written report)
length = 25 # rows = columns
A = np.zeros((length, length))
# print(np.matrix(A))
b = np.zeros(length)
b.shape = (length, 1) # column vector
# print(np.matrix(b))
# print(A.dot(b))


# 2 TODO: Write a function that reads in a file to obtain the resistances on each link
# DONE
def read_resistances_json(file_name):
    resistances = {}
    with open(file_name, 'r') as file:
        data = json.load(file)
        for value in data:
            node1 = value['node1']
            node2 = value['node2']
            resistance = value['resistance']
            resistances[(node1, node2)] = resistance
    return resistances

resistances = read_resistances_json('node_resistances.json')
# print(resistance)


# 3 TODO: Write a function that reads in a file to obtain the set of voltages at fixed nodes
# DONE
def read_voltages_json(file_name):
    voltages = {}
    with open(file_name, 'r') as file:
        data = json.load(file)
        for value in data:
            node = value['node']
            voltage = value['voltage']
            voltages[node] = voltage
    return voltages

voltages = read_voltages_json('node_voltages.json')
# print(voltage)


# 4 TODO: Compute the matrix A using the data read in previously
# IN PROGRESS (ALMOST DONE)
for i in range(length):
    if i + 1 in voltages:
        A[i, i] = 1
        b[i] = voltages[i + 1]
# print(np.matrix(A))
# print(np.matrix(b))

def get_neighbors(node_self):
    neighbors = []
    for key in resistances:
        if key[0] == node_self:
            neighbors.append(key[1])
        elif key[1] == node_self:
            neighbors.append(key[0])
    return neighbors

for node_self in range(1, length + 1):
    neighbors = get_neighbors(node_self)
    # print('Node ' + str(node) + ' is connected to node(s) ' + str(neighbor))
    sum_resistances = 0
    for node_neighbor in neighbors:
        if (node_self, node_neighbor) in resistances:
            resistance = resistances[(node_self, node_neighbor)]
        else:
            resistance = resistances[(node_neighbor, node_self)]
    # print('The resistance between ' + str(node) + ' and ' + str(neighbor) + ' is ' + str(resistance))
    if resistance > 0:
        A[node_self - 1, node_neighbor - 1] = -1 / resistance
        sum_resistances += 1 / resistance
    A[node_self - 1, node_self - 1] = sum_resistances
    b[node_self - 1] = 0
# print(np.matrix(A))
# print(np.matrix(b))


# 5 TODO: Compute the LU factorization of A


# 6 TODO: Compute & output the node voltages & currents through each link


# 7 TODO: Write the output of the previous three steps to a file


# 8 TODO: Repeat but with a tree/graph network