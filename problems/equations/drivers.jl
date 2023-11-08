function driver(DT::ContGal,       #Space discretization type
                inputs::Dict,      #input parameters from src/user_input.jl
                OUTPUT_DIR::String,
                TFloat) 

    sem = sem_setup(inputs)
    
    qp = initialize(sem.mesh.SD, sem.PT, sem.mesh, inputs, OUTPUT_DIR, TFloat)
    
    check_length(qp.qn[1,:], qp.neqs+1, "drivers --> initialize.jl")
        
    solution = time_loop!(sem.QT, sem.PT, inputs[:SOL_VARS_TYPE], inputs[:CL], sem.mesh, sem.metrics, sem.basis, sem.ω, qp,
                          sem.matrix.M, sem.matrix.Minv,
                          inputs[:Δt],
                          inputs,
                          OUTPUT_DIR,
                          TFloat;fx=sem.fx,fy=sem.fy)
    
    if (inputs[:ndiagnostics_outputs] > 0)        
        varsout = qp.qvars
        if (isempty(inputs[:outvars]) == false)
            varsout = inputs[:outvars]
        end
        write_output(sem.mesh.SD, solution,  sem.mesh, OUTPUT_DIR, inputs, varsout, inputs[:outformat]; nvar=qp.neqs, qexact=qp.qe, case="rtb")
    end
    
end
