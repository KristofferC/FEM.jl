
immutable LinSolver <: Solver
    tol::Float64
    max_iters::Int
end

function solve(solver::LinSolver, fp::FEProblem, exporter::AbstractDataExporter)
    println("Starting Lin solver..")
    load = extload(fp)
    force_imbalance = similar(load)
    iteration = 0
    tstep = 0
    n_print = 0



    @time K = create_sparse_structure(fp::FEProblem)
    @time colptrs = get_colptrs(K, fp::FEProblem)
    t = 0.0

    updatebcs!(fp, t)

    int_f = assemble_intf(fp)
    force_imbalance = - int_f
    @devec force_imbalance[:] = load .- int_f


            residual = norm(force_imbalance) / norm(load)


        println("\n\t\tIteration $iteration, relative residual $residual")

        assembleK!(K, fp, colptrs)


        du = cholfact(Symmetric(K, :L)) \ force_imbalance

        updatedofs!(fp, du)
         int_f = assemble_intf(fp)
    update_feproblem(fp)

        write_data(fp, exporter, n_print)


end

