include("./initialize.jl")

function driver(DT::ContGal,       #Space discretization type
                inputs::Dict,      #input parameters from src/user_input.jl
                OUTPUT_DIR::String,
                TFloat) 

    sem = sem_setup(inputs)
    
    @time qp = initialize(sem.mesh.SD, sem.PT, sem.mesh, inputs, OUTPUT_DIR, TFloat)
    
    solution = time_loop!(sem.QT, sem.PT, sem.mesh, sem.metrics, sem.basis, sem.ω, qp,
                          sem.matrix.M, sem.matrix.De, sem.matrix.Le,
                          inputs[:Δt],
                          inputs,
                          OUTPUT_DIR,
                          TFloat)

    if (inputs[:lwrite_output])
        write_output(solution, sem.mesh.SD, sem.mesh, OUTPUT_DIR, inputs, inputs[:outformat];
                 nvar=qp.neqs, qexact=qp.qe, case="rtb")
    end
    
end
