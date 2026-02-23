# Truss Analysis Solver - User Guide

## Overview
This package includes two Python scripts for Finite Element Method (FEM) analysis on 2D truss structures:

1. **truss_solver.py** - Edit the code directly to change inputs
2. **truss_solver_interactive.py** - Interactive version that asks for inputs when you run it

## Requirements
- Python 3.x
- NumPy library


## How to Use

### Option 1: Interactive Version
Run the interactive script:
```
python truss_solver_interactive.py
```

The program will:
1. Ask you to choose between example problem or custom input
2. Guide you through entering all required data:
   - Material properties (E, A)
   - Node coordinates
   - Element connectivity
   - Boundary conditions
   - External forces
3. Solve and display results
4. Optionally save results to JSON file


### Option 2: Direct Code Editing
Simply run the script:
```
python truss_solver.py
```

### 2. Editing Input Parameters

#### Material Properties
```python
E = 200*10**3   # Young's modulus (MPa)
A = 500         # Cross-sectional area (mm^2)
```

#### Node Coordinates
Define your node positions (0-indexed):
```python
nodes = np.array([
    [0.0, 0.0],   # Node 1 (index 0)
    [0.0, 2.0],   # Node 2 (index 1)
    [2.0, 0.0],   # Node 3 (index 2)
    # ... add more nodes
])
```

#### Element Connectivity
Define which nodes are connected (0-indexed):
```python
elements = np.array([
    [0, 1],  # Element connecting node 0 to node 1
    [0, 2],  # Element connecting node 0 to node 2
    # ... add more elements
])
```

#### Boundary Conditions
Define fixed supports (0-indexed for nodes):
```python
bc = np.array([
    [1, 1, 1],  # Node at index 1: fixed in X (1) and Y (1)
    [4, 0, 1]   # Node at index 4: free in X (0), fixed in Y (1)
])
```
Format: `[node_index, fix_x, fix_y]` where `1 = fixed`, `0 = free`

#### External Forces
Define forces at each node (N):
```python
forces = np.array([
    [0.0, -100000.0],   # Node 1: Fx=0, Fy=-100kN
    [0.0, 0.0],         # Node 2: no force
    # ... one force pair per node
])
```

### 3. Understanding the Output

The script prints:
- **System Information**: DOFs, nodes, elements
- **Nodal Displacements**: Movement of each node in X and Y directions (mm)
- **Reaction Forces**: Forces at supports (N)
- **Internal Forces**: Axial forces in each element with TENSION/COMPRESSION indication (N)


## Modifying the Code

### Adding New Calculations
You can modify the `solve_truss()` function to add:
- Stress calculations: `stress = internal_force / A`
- Strain calculations: `strain = stress / E`
- Safety factors
- Additional output parameters

### Changing the Algorithm
Edit the calculation sections:
- Line 34-43: Stiffness matrix assembly
- Line 60-64: Boundary condition application
- Line 75-77: Displacement solver
- Line 83-102: Internal force calculation

### Author
Ahmet Hakan KARAAHMET
