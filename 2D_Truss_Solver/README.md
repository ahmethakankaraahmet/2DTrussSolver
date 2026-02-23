# Truss Analysis Solver - User Guide

## Overview
This package includes two Python scripts for Finite Element Method (FEM) analysis on 2D truss structures:

1. **truss_solver.py** - Edit the code directly to change inputs
2. **truss_solver_interactive.py** - Interactive version that asks for inputs when you run it

## Requirements
- Python 3.x
- NumPy library

### Installation
```bash
pip install numpy
```

## How to Use

### Option 1: Interactive Version (Recommended for Beginners)
Run the interactive script:
```bash
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

**Example Session:**
```
Choose an option:
1. Use example problem (5 nodes, 7 elements)
2. Enter custom problem
3. Exit

Your choice (1/2/3): 2

--- MATERIAL PROPERTIES ---
Enter Young's Modulus (MPa) [default: 200000]: 200000
Enter Cross-sectional Area (mm²) [default: 500]: 500

--- NODE COORDINATES ---
Node 0 (or 'done'): 0,0
Node 1 (or 'done'): 0,2
Node 2 (or 'done'): 2,0
Node 3 (or 'done'): done

--- ELEMENT CONNECTIVITY ---
Element 0 (or 'done'): 0,1
Element 1 (or 'done'): 1,2
Element 2 (or 'done'): done

... (continues with boundary conditions and forces)
```

### Option 2: Direct Code Editing
Simply run the script:
```bash
python truss_solver.py
```

### 2. Editing Input Parameters

Open `truss_solver.py` in any text editor and modify the values in the `if __name__ == "__main__":` section:

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

### 4. Advanced Usage

#### Using as a Module
Import the function in your own script:
```python
from truss_solver import solve_truss
import numpy as np

# Define your inputs
E = 200e3
A = 500
nodes = np.array([[0,0], [1,0], [1,1]])
elements = np.array([[0,1], [1,2]])
bc = np.array([[0,1,1]])
forces = np.array([[0,0], [0,-1000], [0,0]])

# Solve
results = solve_truss(E, A, nodes, elements, bc, forces)

# Access results
displacements = results['displacements']
internal_forces = results['internal_forces']
```

#### Saving Results to JSON
Uncomment the code at the end of the script to save results to a JSON file:
```python
with open('truss_results.json', 'w') as f:
    json.dump(results_serializable, f, indent=2)
```

## Example Problem

The default example in the script is a 5-node, 7-element truss with:
- Pinned support at node 2
- Roller support at node 5
- 100 kN downward forces at nodes 1 and 3

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

## Troubleshooting

**Error: "module 'numpy' has no attribute..."**
- Make sure NumPy is installed: `pip install numpy`

**Error: "Singular matrix"**
- Check that your structure has sufficient supports
- Ensure boundary conditions are properly defined
- Verify element connectivity

**Unexpected results:**
- Verify node indexing (0-indexed in code)
- Check force directions (positive = right/up, negative = left/down)
- Confirm boundary conditions (1 = fixed, 0 = free)
- Ensure consistent units (mm, N, MPa)

## Contact & Support
For questions or issues, refer to your course materials or consult with your instructor.
