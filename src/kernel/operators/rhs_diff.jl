#---------------------------------------------------------------------------
# Fetch equations name to access the user_rhs functions
#---------------------------------------------------------------------------
if (length(ARGS) === 1) #equations
    if isfile(string(@__DIR__, "/../../equations/", ARGS[1], "/user_source.jl"))
        user_source_dir = string(@__DIR__, "/../../equations/", ARGS[1], "/user_source.jl")
    else
        user_source_dir = "../../fallbacks/source.jl"
    end
elseif (length(ARGS) === 2)  #equations/equations_case_name
    user_flux_dir   = string("../../equations/", ARGS[1], "/", ARGS[2], "/user_flux.jl")
    if isfile(string(@__DIR__, "/../../equations/", ARGS[1], "/", ARGS[2], "/user_source.jl"))
        user_source_dir = string(@__DIR__, "/../../equations/", ARGS[1], "/", ARGS[2], "/user_source.jl")
    else
        @info " user_source.jl not defined. The fallback ../../fallbacks/source.jl will be used."
        user_source_dir =  "../../fallbacks/source.jl"
    end
end
include(user_flux_dir)
include(user_source_dir)
include("../ArtificialViscosity/DynSGS.jl")
#---------------------------------------------------------------------------

#
# CompEuler
#
function build_rhs_diff_work_array()
    
    ρel = zeros(mesh.ngl, mesh.ngl)
    uel = zeros(mesh.ngl, mesh.ngl)
    Tel = zeros(mesh.ngl, mesh.ngl)
    Eel = zeros(mesh.ngl, mesh.ngl)
    
end


#function build_rhs_diff!(rhs_diff::SubArray{Float64}, SD::NSD_2D, QT, PT::CompEuler, qp, neqs, basis, ω, inputs, mesh::St_mesh, metrics::St_metrics, μ, T; qoutauxi=zeros(1,1))   
function viscous_rhs_el!(rhs_diff_el, uaux, u, SD::NSD_2D, mesh, metrics, basis, ω, neqs, inputs)#; νx=0.0, νy=0.0)
    
    #rhsdiffξ_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem, neqs)
    #rhsdiffη_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem, neqs)

    ρel = zeros(Float64, mesh.ngl, mesh.ngl)
    uel = zeros(Float64, mesh.ngl, mesh.ngl)
    vel = zeros(Float64, mesh.ngl, mesh.ngl)
    Tel = zeros(Float64, mesh.ngl, mesh.ngl)
        
    #
    # qp[1:npoin]         <-- qq[1:npoin, "ρ"]
    # qp[npoin+1:2npoin]  <-- qq[1:npoin, "ρu"]
    # qp[2npoin+1:3npoin] <-- qq[1:npoin, "ρE"]
    #
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        uaux[:,i] = view(u, idx+1:i*mesh.npoin)
    end
    
    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    if (inputs[:case] === "rtb")
        δenergy = 0.0
    else
        δenergy = 1.0
    end

  #=  
    μ = inputs[:νx]
    ν = 0.0
    κ = μ
    γ = 1.4
    Pr = 0.1
    for iel=1:mesh.nelem
        for j=1:mesh.ngl, i=1:mesh.ngl
            
            m = mesh.connijk[iel,i,j]

            ρel[i,j] = uaux[m,1]          
            uel[i,j] = uaux[m,2]/ρel[i,j]
            vel[i,j] = uaux[m,3]/ρel[i,j]
            Tel[i,j] = uaux[m,4]/ρel[i,j] - δenergy*0.5*(uel[i,j]^2 + vel[i,j]^2)
            
        end

    end
        #ν = Pr*μ[iel]/maximum(ρel[:,:,iel])
        #κ = Pr*μ[iel]/(γ - 1.0)
        #ν = μ[iel]#10.0
        #κ = μ[iel]#10.0
       =#
    #=    for l = 1:params.mesh.ngl
            for k = 1:params.mesh.ngl
                
                ωJkl = params.ω[k]*params.ω[l]*params.metrics.Je[k, l, iel]
                
                dρdξ = 0.0
                dudξ = 0.0
                dvdξ = 0.0
                dTdξ = 0.0

                dρdη = 0.0
                dudη = 0.0
                dvdη = 0.0
                dTdη = 0.0
                for i = 1:params.mesh.ngl
                    dρdξ += params.basis.dψ[i,k]*params.uaux_el[1,i,l,iel]
                    dudξ += params.basis.dψ[i,k]*params.uaux_el[2,i,l,iel]
                    dvdξ += params.basis.dψ[i,k]*params.uaux_el[3,i,l,iel]
                    dTdξ += params.basis.dψ[i,k]*params.uaux_el[4,i,l,iel]

                    dρdη += params.basis.dψ[i,l]*params.uaux_el[1,k,i,iel]
                    dudη += params.basis.dψ[i,l]*params.uaux_el[2,k,i,iel]
                    dvdη += params.basis.dψ[i,l]*params.uaux_el[3,k,i,iel]
                    dTdη += params.basis.dψ[i,l]*params.uaux_el[4,k,i,iel]
                end
                                
                #
                dξdx_kl = params.metrics.dξdx[k,l,iel]
                dξdy_kl = params.metrics.dξdy[k,l,iel]
                dηdx_kl = params.metrics.dηdx[k,l,iel]
                dηdy_kl = params.metrics.dηdy[k,l,iel]
                
                #
                dρdx =  ν*(dρdξ*dξdx_kl + dρdη*dηdx_kl)
                dudx =  μ*(dudξ*dξdx_kl + dudη*dηdx_kl)
                dvdx =  μ*(dvdξ*dξdx_kl + dvdη*dηdx_kl)
                dTdx =  κ*(dTdξ*dξdx_kl + dTdη*dηdx_kl) #+μ∇u⋅u
                
                dρdy =  ν*(dρdξ*dξdy_kl + dρdη*dηdy_kl)
                dudy =  μ*(dudξ*dξdy_kl + dudη*dηdy_kl)
                dvdy =  μ*(dvdξ*dξdy_kl + dvdη*dηdy_kl)
                dTdy =  κ*(dTdξ*dξdy_kl + dTdη*dηdy_kl) #+μ∇u⋅u
                
                ∇ξ∇ρ_kl = dξdx_kl*dρdx + dξdy_kl*dρdy
                ∇η∇ρ_kl = dηdx_kl*dρdx + dηdy_kl*dρdy
                
                ∇ξ∇u_kl = dξdx_kl*dudx + dξdy_kl*dudy
                ∇η∇u_kl = dηdx_kl*dudx + dηdy_kl*dudy            
                ∇ξ∇v_kl = dξdx_kl*dvdx + dξdy_kl*dvdy
                ∇η∇v_kl = dηdx_kl*dvdx + dηdy_kl*dvdy

                ∇ξ∇T_kl = dξdx_kl*dTdx + dξdy_kl*dTdy
                ∇η∇T_kl = dηdx_kl*dTdx + dηdy_kl*dTdy
                
                for i = 1:params.mesh.ngl
                    
                    dhdξ_ik, dhdη_il = params.basis.dψ[i,k], params.basis.dψ[i,l]
                    
                    params.rhs_diffξ_el[i,l,iel,1] -= ωJkl*dhdξ_ik*∇ξ∇ρ_kl
                    params.rhs_diffη_el[k,i,iel,1] -= ωJkl*dhdη_il*∇η∇ρ_kl
                    
                    params.rhs_diffξ_el[i,l,iel,2] -= ωJkl*dhdξ_ik*∇ξ∇u_kl
                    params.rhs_diffη_el[k,i,iel,2] -= ωJkl*dhdη_il*∇η∇u_kl
                    
                    params.rhs_diffξ_el[i,l,iel,3] -= ωJkl*dhdξ_ik*∇ξ∇v_kl
                    params.rhs_diffη_el[k,i,iel,3] -= ωJkl*dhdη_il*∇η∇v_kl
                    
                    params.rhs_diffξ_el[i,l,iel,4] -= ωJkl*dhdξ_ik*∇ξ∇T_kl
                    params.rhs_diffη_el[k,i,iel,4] -= ωJkl*dhdη_il*∇η∇T_kl
                    
                end
            end
        end
    end

    params.rhs_diff_el .= @views (params.rhs_diffξ_el[:,:,:,:] + params.rhs_diffη_el[:,:,:,:])
    =#
end
