abstract Solver

immutable NRSolver <: Solver
    tol::Float64
    max_iters::Int
end

function solve(solver::NRSolver, fp::FEProblem)
    println("Starting Newton-Raphson solver..")

    load = extload(fp)
    iteration = 0
    while true
        iteration += 1
        if iteration > solver.max_iters
            warn("Max iterations hit.")
            break
        end

        int_f = assemble_intf(fp)
        force_unbalance = load - int_f
        residual = norm(force_unbalance) / norm(load)

        @printf("\t\tIteration %d, relative residual %e\n", iteration, residual)

        if residual < solver.tol
            println("Converged!")
            break
        end

        K = assembleK(fp)
        du = K \ force_unbalance

        updatedofs!(fp, du)
    end
end
