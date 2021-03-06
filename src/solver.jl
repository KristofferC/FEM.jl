abstract Solver

# A Newton Raphsoin solver
immutable NRSolver <: Solver
    rel_tol::Float64
    abs_tol::Float64
    alt_tol::Float64
    max_iters::Int
    err_on_nonconv::Bool
end

NRSolver(;rel_tol=0.0::Float64, abs_tol=0.0::Float64,
          max_iters=10::Int, err_on_nonconv=true::Bool, alt_tol = 0.0) =
    NRSolver(rel_tol, abs_tol, alt_tol, max_iters, err_on_nonconv)

# Solves a nonlinear problem by iterating until convergence
function solve(solver::NRSolver, fp::FEProblem, exporter::AbstractDataExporter)
    println("Starting Newton-Raphson solver..")
    load = extload(fp)
    force_imbalance = similar(load)
    iteration = 0
    tstep = 0
    n_print = 0

    K = create_sparse_structure(fp::FEProblem)
    colptrs = get_colptrs(K, fp::FEProblem)

    for t in 0.0 #TODO Implement proper time stepping
        println("Current time $t")
        iteration = 0
        tstep += 1
        updatebcs!(fp, t)
        while true
            iteration += 1
            if iteration > solver.max_iters
                warn("Max iterations hit.")
                if solver.err_on_nonconv
                    error("Did not converge")
                end
                break
            end

            int_f = assemble_intf(fp)
            @devec force_imbalance[:] = load .- int_f


            abs_res = norm(force_imbalance)
            rel_res = norm(force_imbalance) / norm(load)

            @printf("\tIteration %d, relative residual %1.4e, absolute residual %1.4e", iteration, rel_res, abs_res)


            print("\n\tDofType residuals:\n \t\t")
            for (dof_type, eqs) in fp.doftype_eqs
                @devec r = sum(sqr(int_f[eqs]))
                @printf("%s: %1.2e   \t", dof_type, sqrt(r))
            end

            if  abs_res < solver.abs_tol
                println("\nConverged with absolute tolerance!")
                break
            end

            if rel_res < solver.rel_tol
                println("\nConverged with relative tolerance!")
                break
            end

            assembleK!(K, fp, colptrs, fp.dof_vals)


            #du = cholfact(Symmetric(K, :L)) \ force_imbalance
            du = K \ force_imbalance

            print("\n\tDofType updates:\n \t\t")
            for (dof_type, eqs) in fp.doftype_eqs
                @devec r = sum(sqr(du[eqs]))
                @printf("%s: %1.2e   \t", dof_type, sqrt(r))
            end
            updatedofs!(fp, du)

            if norm(du) < solver.alt_tol
                     println("\n Converged with alternative tolerance!")
                break
            end

            print("\n\n")
        end
    update_feproblem(fp)

    n_print += 1
    write_data(fp, exporter, n_print)
    end

end

