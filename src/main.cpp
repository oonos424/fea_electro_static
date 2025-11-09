#include "electrostatics.h"
#include <iostream>

int main() {
    std::cout << "=== 2D Electrostatics FEM Solver ===" << std::endl;
    std::cout << std::endl;
    
    // Create a rectangular mesh
    double width = 1.0;   // 1 meter
    double height = 0.5;  // 0.5 meter
    int nx = 21;          // Number of nodes in x-direction
    int ny = 11;          // Number of nodes in y-direction
    
    Mesh2D mesh;
    std::cout << "Creating mesh..." << std::endl;
    mesh.createRectangularMesh(width, height, nx, ny);
    mesh.printInfo();
    std::cout << std::endl;
    
    // Create solver
    ElectrostaticsSolver solver(mesh);
    
    // Set uniform charge density (C/m^3)
    // solver.setChargeDensity(1.0e-6);  // Uncomment to add volume charge
    
    // Apply boundary conditions - simulate parallel plate capacitor
    // Bottom plate: V = 0 (ground)
    // Top plate: V = 100 V
    // Left and right sides: V = 0
    
    std::cout << "Setting boundary conditions..." << std::endl;
    
    // Bottom boundary (y = 0): V = 0
    for (int i = 0; i < nx; i++) {
        int node_id = i;  // Bottom row nodes
        solver.addBoundaryCondition(BoundaryCondition(node_id, BCType::DIRICHLET, 0.0));
    }
    
    // Top boundary (y = height): V = 100
    for (int i = 0; i < nx; i++) {
        int node_id = (ny - 1) * nx + i;  // Top row nodes
        solver.addBoundaryCondition(BoundaryCondition(node_id, BCType::DIRICHLET, 100.0));
    }
    
    // Left boundary (x = 0): V = 0
    for (int j = 0; j < ny; j++) {
        int node_id = j * nx;  // Left column nodes
        solver.addBoundaryCondition(BoundaryCondition(node_id, BCType::DIRICHLET, 0.0));
    }
    
    // Right boundary (x = width): V = 0
    for (int j = 0; j < ny; j++) {
        int node_id = j * nx + (nx - 1);  // Right column nodes
        solver.addBoundaryCondition(BoundaryCondition(node_id, BCType::DIRICHLET, 0.0));
    }
    
    std::cout << std::endl;
    
    // Solve the system
    solver.solve();
    std::cout << std::endl;
    
    // Export solution to VTK
    solver.exportToVTK("electrostatics_solution.vtk");
    std::cout << std::endl;
    
    // Print some solution values
    const auto& solution = solver.getSolution();
    std::cout << "Sample solution values:" << std::endl;
    std::cout << "  Node 0 (corner): " << solution[0] << " V" << std::endl;
    std::cout << "  Node " << nx/2 << " (bottom center): " << solution[nx/2] << " V" << std::endl;
    std::cout << "  Node " << (ny/2)*nx + nx/2 << " (center): " << solution[(ny/2)*nx + nx/2] << " V" << std::endl;
    std::cout << "  Node " << (ny-1)*nx + nx/2 << " (top center): " << solution[(ny-1)*nx + nx/2] << " V" << std::endl;
    
    std::cout << std::endl;
    std::cout << "You can visualize the solution in ParaView by opening electrostatics_solution.vtk" << std::endl;
    
    return 0;
}