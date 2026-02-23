import numpy as np
import json

def solve_truss(E, A, nodes, elements, bc, forces):
    """
    Finite Element Method solver for 2D truss structures
    
    Parameters:
    -----------
    E : float
        Young's modulus (MPa)
    A : float
        Cross-sectional area (mm²)
    nodes : ndarray
        Node coordinates [[x1,y1], [x2,y2], ...]
    elements : ndarray
        Element connectivity [[n1,n2], ...] (0-indexed)
    bc : ndarray
        Boundary conditions [[node, fixX, fixY], ...] (0-indexed for node)
    forces : ndarray
        External forces [[Fx1,Fy1], [Fx2,Fy2], ...] per node
    
    Returns:
    --------
    dict : Dictionary containing all results
    """
    
    # Setup global stiffness matrix (DOFS)
    ndof = 2 * len(nodes)
    K = np.zeros((ndof, ndof))
    F = np.array([forces.flatten()]).T  # Column vector format
    
    # Assemble global stiffness matrix
    for el in elements:
        n1, n2 = el 
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]
        L = np.sqrt((x2-x1)**2 + (y2-y1)**2)
        c = (x2-x1)/L
        s = (y2-y1)/L

        # k**e in global coordinates
        k_local = (E*A/L) * np.array([
            [ c*c,  c*s, -c*c, -c*s],
            [ c*s,  s*s, -c*s, -s*s],
            [-c*c, -c*s,  c*c,  c*s],
            [-c*s, -s*s,  c*s,  s*s]
        ])

        dof_map = [2*n1, 2*n1+1, 2*n2, 2*n2+1]
        for i in range(4):
            for j in range(4):
                K[dof_map[i], dof_map[j]] += k_local[i, j]
    
    # Apply boundary conditions
    fixed_dofs = []
    for i, bx, by in bc:
        if bx: fixed_dofs.append(2*i)
        if by: fixed_dofs.append(2*i+1)

    free_dofs = np.setdiff1d(np.arange(ndof), fixed_dofs)
    
    # Formation of the submatrices
    K_kk = K[np.ix_(fixed_dofs, fixed_dofs)]
    K_ku = K[np.ix_(fixed_dofs, free_dofs)]
    K_uk = K_ku.T 
    K_uu = K[np.ix_(free_dofs, free_dofs)]
    
    F_k = F[free_dofs]
    
    # Solve displacements
    D = np.array([np.zeros(ndof).flatten()]).T
    D_k = D[fixed_dofs]
    D_u = np.linalg.solve(K_uu, F_k) - K_uk @ D_k
    D[free_dofs] = D_u
    
    # Solve the reaction forces
    F_u = np.matmul(K_ku, D_u)
    F[fixed_dofs] = F_u
    
    # Calculate internal forces
    internal_forces = []
    for el in elements:
        n1, n2 = el 
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]
        L = np.sqrt((x2-x1)**2 + (y2-y1)**2)
        c = (x2-x1)/L
        s = (y2-y1)/L

        dN1_dx = -1/L
        dN2_dx = 1/L
        dN_dx = np.array([np.array([dN1_dx, 0, dN2_dx, 0])])

        T = np.array([
            [ c, s, 0, 0],
            [-s, c, 0, 0],
            [ 0, 0, c, s],
            [ 0, 0,-s, c]
        ])

        dof_map = [2*n1, 2*n1+1, 2*n2, 2*n2+1]
        d = T @ D[dof_map]
        internal_force = E*(dN_dx @ d)*A
        internal_forces.append(float(internal_force[0][0]))
    
    # Package results
    results = {
        "displacements": D.flatten(),
        "reactions": F.flatten(),
        "internal_forces": np.array(internal_forces),
        "K_uu": K_uu,
        "K_uk": K_uk,
        "K_ku": K_ku,
        "K_kk": K_kk,
        "K_global": K,
        "fixed_dofs": fixed_dofs,
        "free_dofs": free_dofs,
        "num_nodes": len(nodes),
        "num_elements": len(elements)
    }
    
    return results


def print_results(results):
    """Print formatted results"""
    print("\n" + "="*70)
    print("TRUSS ANALYSIS RESULTS")
    print("="*70)
    
    print("\n--- SYSTEM INFORMATION ---")
    print(f"Total DOFs: {results['num_nodes'] * 2}")
    print(f"Nodes: {results['num_nodes']}")
    print(f"Elements: {results['num_elements']}")
    print(f"Fixed DOFs: {results['fixed_dofs']}")
    print(f"Free DOFs: {results['free_dofs']}")
    
    print("\n--- NODAL DISPLACEMENTS ---")
    for i in range(results['num_nodes']):
        dx = results['displacements'][2*i]
        dy = results['displacements'][2*i+1]
        print(f"Node {i+1}: dx = {dx:12.6e} mm, dy = {dy:12.6e} mm")
    
    print("\n--- REACTION FORCES ---")
    for i in range(results['num_nodes']):
        fx = results['reactions'][2*i]
        fy = results['reactions'][2*i+1]
        if abs(fx) > 1e-6 or abs(fy) > 1e-6:
            print(f"Node {i+1}: Fx = {fx:12.2f} N, Fy = {fy:12.2f} N")
    
    print("\n--- INTERNAL FORCES ---")
    for i, force in enumerate(results['internal_forces']):
        status = "TENSION" if force > 0 else "COMPRESSION"
        print(f"Element {i+1}: {force:12.2f} N ({status})")
    
    print("\n" + "="*70)


if __name__ == "__main__":
    # Material and section properties 
    E = 200*10**3   # Young's modulus (MPa) -> (N/mm^2)
    A = 500         # Cross-sectional area (mm^2)
    
    # Define geometry, connectivity and boundary conditions
    # Node coordinates: [x, y]
    nodes = np.array([
        [0.0, 0.0],   # Node 1
        [0.0, 2.0],   # Node 2
        [2.0, 0.0],
        [2.0, 2.0],
        [4.0, 0.0],
        [4.0, 2.0],
        [6.0, 0.0]
    ]) 

    
    # Element connectivity: [start_node, end_node] (ZERO INDEXING)
    elements = np.array([
        [0, 1],
        [0, 2],
        [1, 2],
        [1, 3],
        [2, 3],
        [2, 4],
        [3, 4],
        [3, 5],
        [4, 5],
        [4, 6],
        [5, 6]
    ])
    
    # Boundary Conditions (ZERO INDEXING FOR THE NODES)
    # Format: [node_index, fix_x, fix_y] where 1 = fixed, 0 = free
    bc = np.array([
        [0, 1, 0],
        [1, 1, 1]
    ])
    
    # Nodal (EXTERNAL) forces [Fx, Fy] at each node (N)
    forces = np.array([
        [0.0, -100.0*10**3],      # Node 1 
        [0.0, 0.0],      # Node 2 
        [0.0, -100.0*10**3],      # Node 3
        [0.0, 0.0],      # Node 4
        [0.0, -100.0*10**3],    # Node 5
        [0.0, 0.0],
        [0.0, -100.0*10**3]
    ])
    
    # Solve the truss
    results = solve_truss(E, A, nodes, elements, bc, forces)
    
    # Print results
    print_results(results)
    
    # Optional: Save results to JSON file
    # Uncomment the following lines to save results
    """
    results_serializable = {
        "displacements": results['displacements'].tolist(),
        "reactions": results['reactions'].tolist(),
        "internal_forces": results['internal_forces'].tolist(),
        "fixed_dofs": results['fixed_dofs'],
        "free_dofs": results['free_dofs'].tolist(),
        "num_nodes": results['num_nodes'],
        "num_elements": results['num_elements']
    }
    
    with open('truss_results.json', 'w') as f:
        json.dump(results_serializable, f, indent=2)
    print("\nResults saved to truss_results.json")
    """
