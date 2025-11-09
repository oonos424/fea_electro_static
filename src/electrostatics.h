#ifndef FEM_ELECTROSTATICS_H
#define FEM_ELECTROSTATICS_H

#include <vector>
#include <array>
#include <string>
#include <cmath>

// Simple 2D point structure
struct Point2D {
    double x, y;
    Point2D(double x_ = 0.0, double y_ = 0.0) : x(x_), y(y_) {}
};

// Triangle element (3 nodes)
struct Element {
    std::array<int, 3> nodes;  // Node indices
    Element(int n0, int n1, int n2) : nodes{n0, n1, n2} {}
};

// Boundary condition type
enum class BCType {
    NONE,
    DIRICHLET,  // Fixed potential
    NEUMANN     // Fixed flux (not implemented yet)
};

// Boundary condition
struct BoundaryCondition {
    int node_id;
    BCType type;
    double value;
    BoundaryCondition(int id, BCType t, double v) 
        : node_id(id), type(t), value(v) {}
};

// 2D FEM Mesh class
class Mesh2D {
public:
    std::vector<Point2D> nodes;
    std::vector<Element> elements;
    
    void createRectangularMesh(double width, double height, int nx, int ny);
    void printInfo() const;
    int getNumNodes() const { return nodes.size(); }
    int getNumElements() const { return elements.size(); }
};

// FEM Solver for electrostatics
class ElectrostaticsSolver {
public:
    ElectrostaticsSolver(const Mesh2D& mesh);
    
    // Set charge density (source term)
    void setChargeDensity(double rho);
    
    // Add boundary condition
    void addBoundaryCondition(const BoundaryCondition& bc);
    
    // Assemble and solve the system
    void solve();
    
    // Get solution
    const std::vector<double>& getSolution() const { return solution; }
    
    // Export solution to VTK format
    void exportToVTK(const std::string& filename) const;
    
private:
    const Mesh2D& mesh;
    std::vector<double> solution;
    std::vector<BoundaryCondition> boundary_conditions;
    double charge_density;
    
    // Sparse matrix storage (CSR format)
    std::vector<double> K_values;
    std::vector<int> K_col_indices;
    std::vector<int> K_row_ptr;
    std::vector<double> rhs;
    
    void assembleSystem();
    void applyBoundaryConditions();
    void solveLinearSystem();
    
    // Element stiffness matrix computation
    void computeElementMatrix(int elem_id, 
                             std::array<std::array<double, 3>, 3>& Ke) const;
    
    // Gauss elimination for sparse system (simple solver)
    void solveGaussElimination();
};

#endif // FEM_ELECTROSTATICS_H