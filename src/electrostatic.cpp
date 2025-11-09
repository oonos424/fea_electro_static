#include "electrostatics.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

const double EPSILON_0 = 8.854187817e-12;  // Permittivity of free space

// ============================================================================
// Mesh2D Implementation
// ============================================================================

void Mesh2D::createRectangularMesh(double width, double height, int nx, int ny) {
    nodes.clear();
    elements.clear();
    
    // Generate nodes
    double dx = width / (nx - 1);
    double dy = height / (ny - 1);
    
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            nodes.emplace_back(i * dx, j * dy);
        }
    }
    
    // Generate triangular elements (2 triangles per quad)
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
            int n0 = j * nx + i;
            int n1 = j * nx + (i + 1);
            int n2 = (j + 1) * nx + i;
            int n3 = (j + 1) * nx + (i + 1);
            
            // Lower triangle
            elements.emplace_back(n0, n1, n2);
            // Upper triangle
            elements.emplace_back(n1, n3, n2);
        }
    }
}

void Mesh2D::printInfo() const {
    std::cout << "Mesh Information:" << std::endl;
    std::cout << "  Number of nodes: " << nodes.size() << std::endl;
    std::cout << "  Number of elements: " << elements.size() << std::endl;
}

// ============================================================================
// ElectrostaticsSolver Implementation
// ============================================================================

ElectrostaticsSolver::ElectrostaticsSolver(const Mesh2D& mesh_)
    : mesh(mesh_), charge_density(0.0) {
    solution.resize(mesh.getNumNodes(), 0.0);
    rhs.resize(mesh.getNumNodes(), 0.0);
}

void ElectrostaticsSolver::setChargeDensity(double rho) {
    charge_density = rho;
}

void ElectrostaticsSolver::addBoundaryCondition(const BoundaryCondition& bc) {
    boundary_conditions.push_back(bc);
}

void ElectrostaticsSolver::solve() {
    std::cout << "Assembling system..." << std::endl;
    assembleSystem();
    
    std::cout << "Applying boundary conditions..." << std::endl;
    applyBoundaryConditions();
    
    std::cout << "Solving linear system..." << std::endl;
    solveLinearSystem();
    
    std::cout << "Solution completed!" << std::endl;
}

void ElectrostaticsSolver::computeElementMatrix(int elem_id, 
                                                std::array<std::array<double, 3>, 3>& Ke) const {
    const Element& elem = mesh.elements[elem_id];
    
    // Get node coordinates
    double x1 = mesh.nodes[elem.nodes[0]].x;
    double y1 = mesh.nodes[elem.nodes[0]].y;
    double x2 = mesh.nodes[elem.nodes[1]].x;
    double y2 = mesh.nodes[elem.nodes[1]].y;
    double x3 = mesh.nodes[elem.nodes[2]].x;
    double y3 = mesh.nodes[elem.nodes[2]].y;
    
    // Compute element area (using cross product)
    double area = 0.5 * std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
    
    if (area < 1e-12) {
        std::cerr << "Warning: Degenerate element " << elem_id << std::endl;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ke[i][j] = 0.0;
        return;
    }
    
    // Shape function derivatives (constant for linear triangles)
    double b1 = y2 - y3;
    double b2 = y3 - y1;
    double b3 = y1 - y2;
    
    double c1 = x3 - x2;
    double c2 = x1 - x3;
    double c3 = x2 - x1;
    
    // Element stiffness matrix: Ke_ij = (1/(4*Area)) * (b_i*b_j + c_i*c_j)
    double coeff = 1.0 / (4.0 * area);
    
    Ke[0][0] = coeff * (b1 * b1 + c1 * c1);
    Ke[0][1] = coeff * (b1 * b2 + c1 * c2);
    Ke[0][2] = coeff * (b1 * b3 + c1 * c3);
    
    Ke[1][0] = Ke[0][1];  // Symmetric
    Ke[1][1] = coeff * (b2 * b2 + c2 * c2);
    Ke[1][2] = coeff * (b2 * b3 + c2 * c3);
    
    Ke[2][0] = Ke[0][2];  // Symmetric
    Ke[2][1] = Ke[1][2];
    Ke[2][2] = coeff * (b3 * b3 + c3 * c3);
}

void ElectrostaticsSolver::assembleSystem() {
    int n = mesh.getNumNodes();
    
    // Initialize global matrix as dense (for simplicity)
    // In production code, use sparse matrix format
    std::vector<std::vector<double>> K_dense(n, std::vector<double>(n, 0.0));
    std::fill(rhs.begin(), rhs.end(), 0.0);
    
    // Loop over all elements
    for (int e = 0; e < mesh.getNumElements(); e++) {
        std::array<std::array<double, 3>, 3> Ke;
        computeElementMatrix(e, Ke);
        
        const Element& elem = mesh.elements[e];
        
        // Get element area for load vector
        double x1 = mesh.nodes[elem.nodes[0]].x;
        double y1 = mesh.nodes[elem.nodes[0]].y;
        double x2 = mesh.nodes[elem.nodes[1]].x;
        double y2 = mesh.nodes[elem.nodes[1]].y;
        double x3 = mesh.nodes[elem.nodes[2]].x;
        double y3 = mesh.nodes[elem.nodes[2]].y;
        double area = 0.5 * std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
        
        // Assemble into global system
        for (int i = 0; i < 3; i++) {
            int gi = elem.nodes[i];  // Global node index
            
            // Load vector (charge density): fe_i = (rho/epsilon_0) * Area/3
            rhs[gi] += (-charge_density / EPSILON_0) * area / 3.0;
            
            for (int j = 0; j < 3; j++) {
                int gj = elem.nodes[j];
                K_dense[gi][gj] += Ke[i][j];
            }
        }
    }
    
    // Store as dense for now (convert to sparse in production code)
    K_values.clear();
    K_col_indices.clear();
    K_row_ptr.clear();
    K_row_ptr.push_back(0);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (std::abs(K_dense[i][j]) > 1e-14) {
                K_values.push_back(K_dense[i][j]);
                K_col_indices.push_back(j);
            }
        }
        K_row_ptr.push_back(K_values.size());
    }
}

void ElectrostaticsSolver::applyBoundaryConditions() {
    // Apply Dirichlet boundary conditions
    // Method: Set row to identity and RHS to boundary value
    
    for (const auto& bc : boundary_conditions) {
        if (bc.type == BCType::DIRICHLET) {
            int node = bc.node_id;
            
            // Find all entries in this row
            int row_start = K_row_ptr[node];
            int row_end = K_row_ptr[node + 1];
            
            // Zero out the row except diagonal
            for (int idx = row_start; idx < row_end; idx++) {
                if (K_col_indices[idx] == node) {
                    K_values[idx] = 1.0;  // Diagonal
                } else {
                    K_values[idx] = 0.0;
                }
            }
            
            // Set RHS
            rhs[node] = bc.value;
        }
    }
}

void ElectrostaticsSolver::solveLinearSystem() {
    // Simple Gauss-Seidel iterative solver
    // For production, use better solvers (CG, GMRES, direct solvers)
    
    int n = solution.size();
    int max_iter = 10000;
    double tolerance = 1e-8;
    
    std::vector<double> x_old = solution;
    
    for (int iter = 0; iter < max_iter; iter++) {
        double max_diff = 0.0;
        
        for (int i = 0; i < n; i++) {
            double sigma = 0.0;
            double diag = 1.0;
            
            int row_start = K_row_ptr[i];
            int row_end = K_row_ptr[i + 1];
            
            for (int idx = row_start; idx < row_end; idx++) {
                int j = K_col_indices[idx];
                if (j != i) {
                    sigma += K_values[idx] * solution[j];
                } else {
                    diag = K_values[idx];
                }
            }
            
            if (std::abs(diag) > 1e-14) {
                solution[i] = (rhs[i] - sigma) / diag;
            }
            
            max_diff = std::max(max_diff, std::abs(solution[i] - x_old[i]));
        }
        
        x_old = solution;
        
        if (max_diff < tolerance) {
            std::cout << "  Converged in " << iter + 1 << " iterations" << std::endl;
            std::cout << "  Residual: " << max_diff << std::endl;
            return;
        }
    }
    
    std::cout << "  Warning: Did not converge in " << max_iter << " iterations" << std::endl;
}

void ElectrostaticsSolver::exportToVTK(const std::string& filename) const {
    std::ofstream file(filename);
    
    file << "# vtk DataFile Version 3.0\n";
    file << "FEM Electrostatics Solution\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    
    // Write points
    file << "POINTS " << mesh.getNumNodes() << " double\n";
    for (const auto& node : mesh.nodes) {
        file << node.x << " " << node.y << " 0.0\n";
    }
    
    // Write cells
    file << "\nCELLS " << mesh.getNumElements() << " " << mesh.getNumElements() * 4 << "\n";
    for (const auto& elem : mesh.elements) {
        file << "3 " << elem.nodes[0] << " " << elem.nodes[1] << " " << elem.nodes[2] << "\n";
    }
    
    // Write cell types (5 = triangle)
    file << "\nCELL_TYPES " << mesh.getNumElements() << "\n";
    for (int i = 0; i < mesh.getNumElements(); i++) {
        file << "5\n";
    }
    
    // Write solution data
    file << "\nPOINT_DATA " << mesh.getNumNodes() << "\n";
    file << "SCALARS potential double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& val : solution) {
        file << val << "\n";
    }
    
    file.close();
    std::cout << "Solution exported to " << filename << std::endl;
}