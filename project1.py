import numpy as np
import json

# 1 TODO: Write a system of equations Ax = b to solve for the node voltages
# DONE
length = 25 # rows = columns = # of nodes
A = np.zeros((length, length))
b = np.zeros(length)
b.shape = (length, 1) # column vector


# 2 TODO: Write a function that reads in a file to obtain the resistances on each link
# DONE
def read_resistances_json(file_name):
    resistances = {}
    with open(file_name, 'r') as file:
        data = json.load(file)
        for value in data:
            node1 = value['node1'] # 1st element in key tuple
            node2 = value['node2'] # 2nd element in key tuple
            resistance = value['resistance']
            resistances[(node1, node2)] = resistance # assign resistance value to key i.e. pair of nodes
    return resistances

resistances = read_resistances_json('node_resistances.json')
# print(resistances)


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
# print(voltages)


# 4 TODO: Compute the matrix A using the data read in previously
# DONE
def get_neighbors(node_self):
    neighbors = []
    for key in resistances:
        if key[0] == node_self:
            neighbors.append(key[1])
        elif key[1] == node_self:
            neighbors.append(key[0])
    return neighbors

def calculate_A(resistances, voltages, length):
    for node_self in range(1, length + 1):
        if node_self in voltages: # node voltage already given by node_voltages.json
            A[node_self - 1, node_self - 1] = 1
            b[node_self - 1] = voltages[node_self]
        else: # node voltage unknown; must apply KCL
            neighbors = get_neighbors(node_self)
            sum_resistances = 0
            for node_neighbor in neighbors:
                if (node_self, node_neighbor) in resistances:
                    resistance = resistances[(node_self, node_neighbor)]
                else:
                    resistance = resistances[(node_neighbor, node_self)]
                # print('The resistance between ' + str(node_self) + ' and ' + str(node_neighbor) + ' is ' + str(resistance))
                if resistance > 0:
                    A[node_self - 1, node_neighbor - 1] = -1 / resistance
                    sum_resistances += 1 / resistance
                A[node_self - 1, node_self - 1] = sum_resistances
                b[node_self - 1] = 0
    return A, b

A, b = calculate_A(resistances, voltages, length)
# print(np.matrix(np.round(A, 3)))
# print(np.matrix(np.round(b, 3)))


# 5 TODO: Compute the LU factorization of A
# DONE
L = np.identity(length) # 25x25 identity matrix
U = A.copy() # start with A

def LU_decomposition(L, U, length):
    for i in range(length - 1):
        for j in range(i + 1, length):
            scalar = U[j, i] / U[i, i] # multiplier used for elimination below the pivot
            L[j, i] = scalar # store in lower triangle of L i.e. inverse of all E's
            for k in range(i, length):
                U[j, k] -= scalar * U[i, k] # zero out the lower triangle column by column, moving down & to the right
    return L, U

L, U = LU_decomposition(L, U, length)
# print(np.matrix(np.round(U, 3)))
# print(np.matrix(np.round(L, 3)))


# 6 TODO: Compute & output the node voltages & currents through each link
# DONE
x = np.linalg.solve(A, b) # voltage at each node
node_voltages = [0] * length # 25 nodes in 5x5 grid
for i in range(len(node_voltages)):
    node_voltages[i] = float(x[i][0])
    # print('The voltage at Node ' + str(i + 1) + ' is ' + str(node_voltages[i]))

link_currents = [] # 40 resistors connecting all nodes (will append one by one)
for i in resistances:
    v1 = node_voltages[i[0] - 1]
    v2 = node_voltages[i[1] - 1]
    current = (v1 - v2) / resistances.get(i) # define current from 1st node to 2nd node in resistances dictionary key
    link_currents.append(current)
    # print('The current from Node ' + str(i[0]) + ' to ' + str(i[1]) + ' is ' + str(current))


# 7 TODO: Write the output of the previous three steps to a file
# DONE
grid_output = open('grid_output.txt', 'w')
grid_output.write('Matrix A (from Step 4):\n')
grid_output.write(str(np.matrix(np.round(A, 3)))) # rounded for readability
grid_output.write('\n\nMatrix L (from Step 5):\n')
grid_output.write(str(np.matrix(np.round(L, 3))))
grid_output.write('\n\nMatrix U (from Step 5):\n')
grid_output.write(str(np.matrix(np.round(U, 3))))
grid_output.write('\n')
for i in range(len(node_voltages)):
    grid_output.write('\nThe voltage at Node ' + str(i + 1) + ' is ' + str(node_voltages[i]))
grid_output.write('\n')
counter = 0
for i in resistances:
    grid_output.write('\nThe current from Node ' + str(i[0]) + ' to ' + str(i[1]) + ' is ' + str(link_currents[counter]))
    counter += 1


# 8 TODO: Repeat but with a tree/graph network
# IN PROGRESS
length1 = 7
A1 = np.zeros((length1, length1))
b1 = np.zeros(length1)
b1.shape = (length1, 1)
resistances1 = {(1, 2): 10.0, (1, 3): 12.0, (2, 4): 2.0, (2, 5): 4.0, (3, 6): 6.0, (3, 7): 8.0}
voltages1 = {1: 0.0, 4: 2.0} # hard-coded for our tree network instead of taking input from json files again

def get_neighbors1(node_self, length): # modified get_neighbors to account for tree structure (avoids IndexError)
    neighbors = []
    for key in resistances:
        if key[0] == node_self and key[1] <= length:
            neighbors.append(key[1])
        elif key[1] == node_self and key[0] <= length:
            neighbors.append(key[0])
    return neighbors

def calculate_A1(resistance, voltages, length): # modified calculate_A to account for tree structure (avoids KeyError)
    for node_self in range(1, length + 1):
        if node_self in voltages: # node voltage already given by node_voltages.json
            A1[node_self - 1, node_self - 1] = 1
            b1[node_self - 1] = voltages[node_self]
        else: # node voltage unknown; must apply KCL
            neighbors = get_neighbors1(node_self, length)
            sum_resistances = 0
            for node_neighbor in neighbors:
                if (node_self, node_neighbor) in resistances:
                    resistance = resistances[(node_self, node_neighbor)]
                elif (node_neighbor, node_self) in resistances:
                    resistance = resistances[(node_neighbor, node_self)]
                else:
                    resistance = None
                # print('The resistance between ' + str(node_self) + ' and ' + str(node_neighbor) + ' is ' + str(resistance))
                if resistance is not None and resistance > 0:
                    A1[node_self - 1, node_neighbor - 1] = -1 / resistance
                    sum_resistances += 1 / resistance
                A1[node_self - 1, node_self - 1] = sum_resistances
                b1[node_self - 1] = 0
    return A1, b1

A1, b1 = calculate_A1(resistances1, voltages1, length1)

L1 = np.identity(length1)
U1 = A1.copy()
L1, U1 = LU_decomposition(L1, U1, length1)

x1 = np.linalg.solve(A1, b1)
node_voltages1 = [0] * length1
for i in range(len(node_voltages1)):
    node_voltages1[i] = float(x1[i][0])
link_currents1 = []
for i in resistances1:
    v1_1 = node_voltages1[i[0] - 1]
    v2_1 = node_voltages1[i[1] - 1]
    current1 = (v1_1 - v2_1) / resistances1.get(i)
    link_currents1.append(current1)

tree_output = open('tree_output.txt', 'w')
tree_output.write('Matrix A (from Step 4):\n')
tree_output.write(str(np.matrix(np.round(A1, 3))))
tree_output.write('\n\nMatrix L (from Step 5):\n')
tree_output.write(str(np.matrix(np.round(L1, 3))))
tree_output.write('\n\nMatrix U (from Step 5):\n')
tree_output.write(str(np.matrix(np.round(U1, 3))))
tree_output.write('\n')
for i in range(len(node_voltages1)):
    tree_output.write('\nThe voltage at Node ' + str(i + 1) + ' is ' + str(node_voltages1[i]))
tree_output.write('\n')
counter = 0
for i in resistances1:
    tree_output.write('\nThe current from Node ' + str(i[0]) + ' to ' + str(i[1]) + ' is ' + str(link_currents1[counter]))
    counter += 1