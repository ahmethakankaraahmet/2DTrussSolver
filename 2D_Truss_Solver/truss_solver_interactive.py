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


def get_user_input():
    """Interactive input collection from user"""
    
    print("="*70)
    print("         TRUSS ANALYSIS - FINITE ELEMENT METHOD")
    print("="*70)
    print("\nThis program will guide you through defining your truss structure.")
    print("Note: All node indices are 0-based (first node is index 0)")
    print()
    
    # Material Properties
    print("\n--- MATERIAL PROPERTIES ---")
    E = float(input("Enter Young's Modulus (MPa) [default: 200000]: ") or "200000")
    E = E * 10**3  # Convert to N/mm²
    A = float(input("Enter Cross-sectional Area (mm²) [default: 500]: ") or "500")
    
    # Node Coordinates
    print("\n--- NODE COORDINATES ---")
    print("Enter node coordinates. Type 'done' when finished.")
    print("Format for each node: x,y")
    print("Example: 0,0")
    
    nodes_list = []
    node_count = 0
    while True:
        node_input = input(f"Node {node_count} (or 'done'): ").strip()
        if node_input.lower() == 'done':
            if len(nodes_list) < 2:
                print("Error: You need at least 2 nodes. Please continue.")
                continue
            break
        try:
            x, y = map(float, node_input.split(','))
            nodes_list.append([x, y])
            node_count += 1
        except:
            print("Invalid format. Please use: x,y (e.g., 0,0)")
    
    nodes = np.array(nodes_list)
    print(f"\nTotal nodes defined: {len(nodes)}")
    
    # Element Connectivity
    print("\n--- ELEMENT CONNECTIVITY ---")
    print("Define elements by connecting nodes (use 0-based indexing).")
    print("Type 'done' when finished.")
    print("Format: startNode,endNode")
    print("Example: 0,1 (connects node 0 to node 1)")
    
    elements_list = []
    element_count = 0
    while True:
        elem_input = input(f"Element {element_count} (or 'done'): ").strip()
        if elem_input.lower() == 'done':
            if len(elements_list) < 1:
                print("Error: You need at least 1 element. Please continue.")
                continue
            break
        try:
            n1, n2 = map(int, elem_input.split(','))
            if n1 < 0 or n1 >= len(nodes) or n2 < 0 or n2 >= len(nodes):
                print(f"Error: Node indices must be between 0 and {len(nodes)-1}")
                continue
            elements_list.append([n1, n2])
            element_count += 1
        except:
            print("Invalid format. Please use: n1,n2 (e.g., 0,1)")
    
    elements = np.array(elements_list)
    print(f"\nTotal elements defined: {len(elements)}")
    
    # Boundary Conditions
    print("\n--- BOUNDARY CONDITIONS ---")
    print("Define fixed supports (1=fixed, 0=free).")
    print("Type 'done' when finished.")
    print("Format: nodeIndex,fixX,fixY")
    print("Example: 0,1,1 (node 0 fixed in both X and Y)")
    print("Example: 3,0,1 (node 3 free in X, fixed in Y - roller)")
    
    bc_list = []
    bc_count = 0
    while True:
        bc_input = input(f"Support {bc_count} (or 'done'): ").strip()
        if bc_input.lower() == 'done':
            if len(bc_list) < 1:
                print("Error: You need at least 1 boundary condition. Please continue.")
                continue
            break
        try:
            node_idx, fix_x, fix_y = map(int, bc_input.split(','))
            if node_idx < 0 or node_idx >= len(nodes):
                print(f"Error: Node index must be between 0 and {len(nodes)-1}")
                continue
            if fix_x not in [0, 1] or fix_y not in [0, 1]:
                print("Error: Fix values must be 0 or 1")
                continue
            bc_list.append([node_idx, fix_x, fix_y])
            bc_count += 1
        except:
            print("Invalid format. Please use: nodeIndex,fixX,fixY (e.g., 0,1,1)")
    
    bc = np.array(bc_list)
    print(f"\nTotal supports defined: {len(bc)}")
    
    # External Forces
    print("\n--- EXTERNAL FORCES ---")
    print("Enter forces for EACH node (N).")
    print("Format: Fx,Fy")
    print("If no force at a node, enter: 0,0")
    
    forces_list = []
    for i in range(len(nodes)):
        while True:
            force_input = input(f"Force at Node {i}: ").strip()
            try:
                fx, fy = map(float, force_input.split(','))
                forces_list.append([fx, fy])
                break
            except:
                print("Invalid format. Please use: Fx,Fy (e.g., 0,-10000)")
    
    forces = np.array(forces_list)
    
    return E, A, nodes, elements, bc, forces


def main():
    """Main program with user input"""
    
    print("\n" + "="*70)
    print("Welcome to the Truss Analysis Solver!")
    print("="*70)
    
    while True:
        print("\nChoose an option:")
        print("1. Use example problem (5 nodes, 7 elements)")
        print("2. Enter custom problem")
        print("3. Exit")
        
        choice = input("\nYour choice (1/2/3): ").strip()
        
        if choice == '1':
            # Use example problem
            E = 200*10**3
            A = 500
            nodes = np.array([
                [0.0, 0.0],
                [0.0, 2.0],
                [2.0, 0.0],
                [2.0, 2.0],
                [4.0, 0.0]
            ])
            elements = np.array([
                [0, 1],
                [0, 2],
                [1, 2],
                [1, 3],
                [2, 3],
                [2, 4],
                [3, 4]
            ])
            bc = np.array([
                [1, 1, 1],
                [4, 0, 1]
            ])
            forces = np.array([
                [0.0, -100.0*10**3],
                [0.0, 0.0],
                [0.0, -100.0*10**3],
                [0.0, 0.0],
                [0.0, 0.0]
            ])
            print("\nUsing example problem...")
            
        elif choice == '2':
            # Get custom input
            E, A, nodes, elements, bc, forces = get_user_input()
            
        elif choice == '3':
            print("\nGoodbye!")
            break
            
        else:
            print("Invalid choice. Please enter 1, 2, or 3.")
            continue
        
        # Solve the truss
        print("\n" + "="*70)
        print("SOLVING...")
        print("="*70)
        
        try:
            results = solve_truss(E, A, nodes, elements, bc, forces)
            print_results(results)
            
            # Ask if user wants to save results
            save_choice = input("\nDo you want to save results to JSON file? (y/n): ").strip().lower()
            if save_choice == 'y':
                filename = input("Enter filename (without .json): ").strip() or "truss_results"
                results_serializable = {
                    "displacements": results['displacements'].tolist(),
                    "reactions": results['reactions'].tolist(),
                    "internal_forces": results['internal_forces'].tolist(),
                    "fixed_dofs": results['fixed_dofs'],
                    "free_dofs": results['free_dofs'].tolist(),
                    "num_nodes": results['num_nodes'],
                    "num_elements": results['num_elements']
                }
                
                with open(f'{filename}.json', 'w') as f:
                    json.dump(results_serializable, f, indent=2)
                print(f"\nResults saved to {filename}.json")
        
        except Exception as e:
            print(f"\nError during analysis: {str(e)}")
            print("Please check your inputs and try again.")
        
        print("\n")


if __name__ == "__main__":
    main()
