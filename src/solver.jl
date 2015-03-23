abstract Solver

immutable NRSolver <: Solver
    tol::Float64
    max_iters::Int
end

function solve(solver::NRSolver, fp::FEProblem)
    println("Starting Newton-Raphson solver..")
    createdofs(fp)
    load = extload(fp)
    iteration = 0



    while true
        iteration += 1
        if iteration > solver.max_iters
            warn("Max iterations hit.")
            break
        end

        # TODO: Check signs and stuff here
        #int_f = assemble_intf(fp)
        force_unbalance = load  #- int_f
        residual = norm(load ) / norm(load)


        println("\t\tIteration $iteration, relative residual $residual")

        #println(load)
        #println(int_f)
        if residual < solver.tol
            println("Converged!")
            break
        end

        K = assembleK(fp)

        println(norm(K - K', 1))

        #print(K)
        #println(force_unbalance)
        du = cholfact(Symmetric(-K, :L)) \ force_unbalance
        #println(du)
        println(force_unbalance' * du)

        updatedofs!(fp, du)
        break
    end
end
