function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_1D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, ip=1)
    γ = 1.4
    γm1 = 0.4

    A  = 1.0 + 2.2*(mesh.x[ip] - 1.5)^2

    ρ  = q[1]/A
    u  = q[2]/q[1]
    T  = γm1*(q[3]/q[1] - 0.5*γ*u*u)
    p  = ρ*T
    
    #F[1] = q[2]
    #F[2] = q[2]*u + p*A/γ
    #F[3] = q[3]*u + p*A*u
    
    # U22oU1 = q[2]*q[2]/q[1]    
    # F[1] = q[2]
    # F[2] = U22oU1 + (γm1/γ)*(q[3] - 0.5*γ*U22oU1)
    # F[3] = γ*q[2]*q[3]/q[1] - 0.5*γ*γm1*q[2]*q[2]*q[2]/(q[1]*q[1])

    #CORRECT###
    F[1] = q[2]
        
    # introduce some useful definitions
    fac1 = (3.0 - γ)/2
    fac2 = (γ - 1.0)/γ
    fac3 = (γ - 1.0)/2
    ratio1 = 1.0/q[1]
    ratio2 = (γ*q[2])/(q[1]*q[1])
    
    term1 = q[1]*q[3]; 
    term2 = q[2]*q[2]
    
    F[2] = ratio1*( fac2 * term1 + fac1*term2 )
    F[3] = ratio2*( term1 - fac3*term2 )
    # @info ip
    # @info q[1]
    # @info q[2]
    # @info q[3]
    # @info F[2]
    # readline();
    
    #@info ip F[2]
end
