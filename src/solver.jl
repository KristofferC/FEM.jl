abstract Solver

immutable NRSolver <: Solver
    tol::Float64
    max_iters::Int
end

function solve(solver::NRSolver, fp::FEProblem, exporter::AbstractDataExporter)
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
        force_imbalance = load - int_f
        residual = norm(force_imbalance) / norm(load)
        if norm(load) < 1.00^(-15)
            residual = norm(int_f)
        else
            residual = norm(force_imbalance) / norm(load)
        end


        println("\t\tIteration $iteration, relative residual $residual")

        if residual < solver.tol
            println("Converged!")
            break
        end


        K = assembleK(fp)

        #du = cholfact(Symmetric(K, :L)) \ force_imbalance

        du = K \ force_imbalance

        println(du)

        updatedofs!(fp, du)

    end


    update_feproblem(fp)
     write_data(fp, exporter)
end

