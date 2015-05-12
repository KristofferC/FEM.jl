abstract Solver

include("linsolver.jl")

immutable NRSolver <: Solver
    rel_tol::Float64
    abs_tol::Float64
    max_iters::Int
    err_on_nonconv::Bool
end

NRSolver(;rel_tol=0.0::Float64, abs_tol=0.0::Float64,
          max_iters=10::Int, err_on_nonconv=true::Bool) =
    NRSolver(rel_tol, abs_tol, max_iters, err_on_nonconv)

function solve(solver::NRSolver, fp::FEProblem, exporter::AbstractDataExporter)
    println("Starting Newton-Raphson solver..")
    load = extload(fp)
    force_imbalance = similar(load)
    iteration = 0
    tstep = 0
    n_print = 0


    @time K = create_sparse_structure(fp::FEProblem)
    @time colptrs = get_colptrs(K, fp::FEProblem)

    for t in 1:5
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
            force_imbalance = - int_f
            @devec force_imbalance[:] = load .- int_f

            #
            #residual = norm(force_imbalance) / norm(load)
           # if norm(load) < 1.00^(-15)
            #residual = norm(int_f)
           # else
               residual = norm(force_imbalance) / norm(load)
           # end


            println("\n\t\tIteration $iteration, relative residual $residual")


            println("\n Residual:")
            for (dof_type, eqs) in fp.doftype_eqs
                @devec r = sum(sqr(int_f[eqs]))
                @printf("%s: %1.2e   \n", dof_type, sqrt(r))
            end


            if residual < solver.abs_tol
                println("Converged!")
                break
            end

           if iteration > 2 && norm(du) < 1e-6
                println("Converged with alt conv!")
                break
            end

            println("assemble")

           # if (tstep - 2) % 5 == 0 || iteration % 5 == 0
                assembleK!(K, fp, colptrs)
           # end


            #du = cholfact(Symmetric(K, :L)) \ force_imbalance

            du = K \ force_imbalance

            println("\nUpdated:")
            for (dof_type, eqs) in fp.doftype_eqs
                @devec r = sum(sqr(du[eqs]))
                @printf("%s: %1.2e   \n", dof_type, sqrt(r))
            end
            updatedofs!(fp, du)



        end
        update_feproblem(fp)

       # if (tstep % 5 == 0)
       #     n_print += 1
       #     write_data(fp, exporter, n_print)
       # end
    end

end

