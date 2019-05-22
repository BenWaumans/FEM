#include "Config.h"
#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <chrono>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
    // 1. Parse command-line options.
    double x = 1;
    double y = 1;
    double z = 1;
    double fx = 0;
    double fy = 0;
    double fz = 0;
    bool static_cond = false;
    bool visualization = true;
    int order = 1;
    double E = 1.0;
    double v = 1.0;

    fprintf(stdout, "%s Version %d.%d.%d\n", argv[0], FEM_SOLVER_VERSION_MAJOR,
            FEM_SOLVER_VERSION_MINOR, FEM_SOLVER_VERSION_PATCH);

    OptionsParser args(argc, argv);
    args.AddOption(&x, "-x", "--x-axis", "Size in the x axis direction.");
    args.AddOption(&y, "-y", "--y-axis", "Size in the y axis direction.");
    args.AddOption(&z, "-z", "--z-axis", "Size in the z axis direction.");
    args.AddOption(&fx, "-fx", "--force-x-axis", "Force in the x axis direction.");
    args.AddOption(&fy, "-fy", "--force-y-axis", "Force in the y axis direction.");
    args.AddOption(&fz, "-fz", "--force-z-axis", "Force in the z axis direction.");
    args.AddOption(&E, "-E", "--young-modulus", "Young's modulus (E) of the material.", true);
    args.AddOption(&v, "-v", "--poisson-ratio", "Poisson's ratio (v) of the material.", true);
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                   "--no-visualization",
                   "Enable or disable GLVis visualization.");
    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(cout);
        return 1;
    }
    args.PrintOptions(cout);

    auto start = std::chrono::high_resolution_clock::now();

    // 2. Read the mesh from the given mesh file. We can handle triangular,
    //    quadrilateral, tetrahedral or hexahedral elements with the same code.
    Mesh *mesh = nullptr;
    {
        double max_value = std::max(std::max(x, y), z);
        mesh = new Mesh(
            ceil(x / max_value * 10.) * 10, ceil(y / max_value * 10.) * 10,
            ceil(z / max_value * 10.) * 10, Element::Type::HEXAHEDRON, 1, x, y, z);
    }
    int dim = mesh->Dimension();

    if (mesh->bdr_attributes.Max() < 2)
    {
        cerr << "\nInput mesh should have at least "
             << "two boundary attributes! (See schematic in test1.cpp)\n"
             << endl;
        return 3;
    }

    // 3. Select the order of the finite element discretization space. For NURBS
    //    meshes, we increase the order by degree elevation.
    if (mesh->NURBSext)
    {
        mesh->DegreeElevate(order, order);
    }

    // 4. Refine the mesh to increase the resolution. In this example we do
    //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
    //    largest number that gives a final mesh with no more than 50,000
    //    elements.
    {
        int ref_levels = (int)floor(log(50000. / mesh->GetNE()) / log(2.) / dim);
        for (int l = 0; l < ref_levels; l++)
        {
            mesh->UniformRefinement();
            cout << "Number of vertices after run " << l << ": "
                 << mesh->GetNV() << std::endl;
        }
    }

    // 5. Define a finite element space on the mesh. Here we use vector finite
    //    elements, i.e. dim copies of a scalar finite element space. The vector
    //    dimension is specified by the last argument of the FiniteElementSpace
    //    constructor. For NURBS meshes, we use the (degree elevated) NURBS
    //    space associated with the mesh nodes.
    FiniteElementCollection *fec;
    FiniteElementSpace *fespace;
    if (mesh->NURBSext)
    {
        fec = NULL;
        fespace = mesh->GetNodes()->FESpace();
    }
    else
    {
        fec = new H1_FECollection(order, dim);
        fespace = new FiniteElementSpace(mesh, fec, dim);
    }
    cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
         << endl
         << "Assembling: " << flush;

    // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
    //    In this example, the boundary conditions are defined by marking only
    //    boundary attribute 1 from the mesh as essential and converting it to a
    //    list of true dofs.
    Array<int> ess_tdof_list, ess_bdr(mesh->bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[0] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    // 7. Set up the linear form b(.) which corresponds to the right-hand side
    // of
    //    the FEM linear system. In this case, b_i equals the boundary integral
    //    of f*phi_i where f represents a "pull down" force on the Neumann part
    //    of the boundary and phi_i are the basis functions in the finite
    //    element fespace. The force is defined by the VectorArrayCoefficient
    //    object f, which is a vector of Coefficient objects. The fact that f is
    //    non-zero on boundary attribute 2 is indicated by the use of piece-wise
    //    constants coefficient for its last component.
    VectorArrayCoefficient f(dim);
    for (int i = 0; i < dim - 1; i++)
    {
        f.Set(i, new ConstantCoefficient(0.0));
    }
    {
        Vector pull_force(mesh->bdr_attributes.Max());
        pull_force = 0.0;
        pull_force(0) = fx;
        pull_force(1) = fy;
        pull_force(2) = fz;
        f.Set(dim - 1, new PWConstCoefficient(pull_force));
    }

    LinearForm *b = new LinearForm(fespace);
    b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));
    cout << "r.h.s. ... " << flush;
    b->Assemble();

    // 7.b Set device config parameters from the command line options and switch
    //    to working on the device.
#ifdef MFEM_USE_CUDA
    Device::Configure("cuda");
    Device::Print();
    Device::Enable();
#endif

    // 8. Define the solution vector x_grid as a finite element grid function
    //    corresponding to fespace. Initialize x_grid with initial guess of zero,
    //    which satisfies the boundary conditions.
    GridFunction x_grid(fespace);
    x_grid = 0.0;

    // 9. Set up the bilinear form a(.,.) on the finite element space
    //    corresponding to the linear elasticity integrator with piece-wise
    //    constants coefficient lambda and mu.
    Vector lambda(mesh->attributes.Max());
    lambda = (E * v) / ((1 + v) * (1 - 2 * v));
    PWConstCoefficient lambda_func(lambda);
    Vector mu(mesh->attributes.Max());
    mu = E / (2 * (1 + v));
    PWConstCoefficient mu_func(mu);

    BilinearForm *a = new BilinearForm(fespace);
    // a->SetAssemblyLevel(AssemblyLevel::PARTIAL);
    a->AddDomainIntegrator(new ElasticityIntegrator(lambda_func, mu_func));

    // 10. Assemble the bilinear form and the corresponding linear system,
    //     applying any necessary transformations such as: eliminating boundary
    //     conditions, applying conforming constraints for non-conforming AMR,
    //     static condensation, etc.
    cout << "matrix ... " << flush;
    if (static_cond)
    {
        a->EnableStaticCondensation();
    }
    a->Assemble();

    SparseMatrix A;
    Vector B, X;
    a->FormLinearSystem(ess_tdof_list, x_grid, *b, A, X, B);
    cout << "done." << endl;

    cout << "Size of linear system: " << A.Height() << endl;

#ifndef MFEM_USE_SUITESPARSE
    // 11. Define a simple symmetric Gauss-Seidel preconditioner and use it to
    //     solve the system Ax=b with PCG.
    GSSmoother M(A);
    PCG(A, M, B, X, 1, 500, 1e-8, 0.0);
#else
    // 11. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the
    // system.
    UMFPackSolver umf_solver;
    umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    umf_solver.SetOperator(A);
    umf_solver.Mult(B, X);
#endif

    // 12. Recover the solution as a finite element grid function.
    a->RecoverFEMSolution(X, *b, x_grid);

    // 12.b Switch back to the host.
#ifdef MFEM_USE_CUDA
    Device::Disable();
#endif

    // 13. For non-NURBS meshes, make the mesh curved based on the finite
    // element
    //     space. This means that we define the mesh elements through a fespace
    //     based transformation of the reference element. This allows us to save
    //     the displaced mesh as a curved mesh when using high-order finite
    //     element displacement field. We assume that the initial mesh (read
    //     from the file) is not higher order curved mesh compared to the chosen
    //     FE space.
    if (!mesh->NURBSext)
    {
        mesh->SetNodalFESpace(fespace);
    }

    std::chrono::duration<float> duration = std::chrono::high_resolution_clock::now() - start;
    std::cout << "It ran for " << duration.count() << "s\n";
    // 15. Send the above data by socket to a GLVis server. Use the "n" and "b"
    //     keys in GLVis to visualize the displacements.
    if (visualization)
    {
        char vishost[] = "localhost";
        int visport = 19916;
        socketstream sol_sock(vishost, visport);
        sol_sock.precision(8);
        sol_sock << "solution\n"
                 << *mesh << x_grid << flush;
    }

    // 16. Free the used memory.
    delete a;
    delete b;
    if (fec)
    {
        delete fespace;
        delete fec;
    }
    delete mesh;

    return 0;
}
