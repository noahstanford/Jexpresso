function user_bc_dirichlet!(u,
                            params, 
                            q::SubArray{Float64},
                            x::AbstractFloat,
                            t::AbstractFloat,
                            tag::String,
                            qbdy::AbstractArray,                            
                            qe::SubArray{Float64},
                            ::TOTAL)
    
    U1in = 0.0
    U2in = 0.0
    U3in = 0.0
    
    U1out = 0.0
    U2out = 0.0
    U3out = 0.0
    
    ip2 = 2 #this is the 2nd point of the linear grid
    ip3 = 3 #this is the 3rd point of the linear grid
    #ipN = length(q[:,1]) #last geometric point of the 1D mesh. The H-O node count starts from this one in the first element.
    ipN = size(q)[1] #last geometric point of the 1D mesh. The H-O node count starts from this one in the first element.
        
    
    xin = 0.0
    Ain = 1.0 + 2.2*(xin - 1.5)^2
    xout = 3.0
    Aout = 1.0 + 2.2*(xout - 1.5)^2
    
    Tin = 1.0
    ρin = 1.0


    pin = ρin*Tin
    mass_flow = 0.579
    uin = mass_flow/(ρin*Ain)
    
    γ = 1.4
    γm1 = 0.4
    
    lshock = false #Notice, only try shock if you have some artificial diffusion implemented
    
    if (tag == "left")
        #### DEBUGGING NOTES
        # values of q[1, 1] and q[1, 2] should  match U[2, 1] and U[1, 1], both are incorrect
        # for qbdy[2], params.RHS should be used, not q
        #build_rhs!(params.RHS, u, params, 0.0)
        #inviscid_rhs_el!(u, params, true, NSD_1D(), FD())
        ####### NEW PROBLEM: 
        U1in = Ain#*ρin
        @info params.RHS[ip2, 2]
        U2in = 2*params.RHS[ip2, 2] - params.RHS[ip3, 2] #corrected version
        #U2in = 2*q[ip2,2] - q[ip3,2] #should not use q, should be using RHS
        #fac3 = params.RHS[1, 2]/params.RHS[1, 1]
        fac3 = q[1,2]/q[1,1] # this is the corrected version
        # @info q[ip3, 2]
        # @info params.RHS[ip3, 2]
        # readline();
        #U3in = U1in*(Tin/γm1 + 0.5*γ*uin*uin)
        U3in = U1in*(1/γm1 + 0.5*γ*(U2in/U1in)^2)
        qbdy[1] = 0#?U1in
        qbdy[2] = U2in
        qbdy[3] = γ*fac3*qbdy[2]#U3in
        # @info "q" q[ip2, 2]
        # readline();
        #@info qbdy[1]
        #@info qbdy[2]
        readline()

       # @info "U2 inner points: " U1in U2in U3in
       #@info qbdy[2]
       #@info qbdy[3]
        
    end

    if (tag == "right")
        pout = 0.6784    
        U1out = 2*q[ipN-1,1] - q[ipN-2,1]
        U2out = 2*q[ipN-1,2] - q[ipN-2,2]
        U3out = 2*q[ipN-1,3] - q[ipN-2,3]

        
        #=if lshock

            U1out = q[ipN,1]
            U2out = q[ipN,2]
            uout = U2out/U1out
            
            U3out = pout*Aout/γm1 + 0.5*γ*U2out*uout
            
        end=#
        
        qbdy[1] = U1out
        qbdy[2] = U2out
        qbdy[3] = U3out
    end
    
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end
