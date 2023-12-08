include("jVKI/IPrtn.jl")
include("jVKI/InitParms.jl")
include("jVKI/Set_XGrid.jl")
include("jVKI/IntegrateODE.jl")

function solveODE(mesh::St_mesh, inputs::Dict, q)
    Gamma = 1.4
    Tot_Int_Pts = mesh.npoin - 2
    Tot_X_Pts = mesh.npoin
    Shock_Flag = 0
    Exit_Pressure = 1
    ithroat = mesh.npoin - 2
    a = 0
    rho = 1

    N, delta1, n, k = InitParms(16, 0.005, 0.005, mesh.npoin) # all parameters based on user inputs
    
    t, hbar = IPrtn(inputs[:tinit], inputs[:tend], n, N)

    x, Del_x = Set_XGrid(mesh.xmin, mesh.xmax, mesh.npoin)

    d = 3#mesh.neqs
    r = 2 #tentative

    InitVal = zeros(Float64, d, Tot_X_Pts)

    # InitVal[1,:] = q[1:mesh.npoin]
    # InitVal[2,:] = q[mesh.npoin+1:2*mesh.npoin]
    # InitVal[3,:] .= 0
    InitVal[1,:] .= q[1:mesh.npoin].+1
    InitVal[2,:] .= q[mesh.npoin+1:2*mesh.npoin].+1
    InitVal[3,:] .= q[2*mesh.npoin+1:3*mesh.npoin].+1

    U2 = zeros(Tot_X_Pts, Int(n+1))

    for gridpt = 1:Tot_X_Pts
        U2[gridpt, 1] = InitVal[2, gridpt]
    end
    
    n = Int(n)

    A = ones(Float64, Tot_X_Pts)
    U2_in = U2
    ff0_throat_in = zeros(Float64, d, n)
    ff1_throat_in = zeros(Float64, d, n)
    ff2_throat_in = InitVal#zeros(Float64, d, n)
    Mach_E = zeros(Float64, Tot_X_Pts)
    Mrho_E = zeros(Float64, Tot_X_Pts)
    Press_E = zeros(Float64, Tot_X_Pts)
    Temp_E = zeros(Float64, Tot_X_Pts)
    Vel_E = zeros(Float64, Tot_X_Pts)
    In_Mass_Flow = 0

    FinalVal, U2, Mach_D, Mrho_D, Press_D, Temp_D, Vel_D, Rel_MachErr, Rel_MrhoErr, Rel_PressErr, Rel_TempErr, Rel_VelErr, AvRelTempErr, AvPlusSDevRelTempErr, 
    AvMinusSDevRelTempErr, AvRelMachErr, AvPlusSDevRelMachErr, AvMinusSDevRelMachErr, AvRelMrhoErr, AvPlusSDevRelMrhoErr, AvMinusSDevRelMrhoErr, AvRelPressErr, 
    AvPlusSDevRelPressErr, AvMinusSDevRelPressErr, AvU2, ff0_throat, ff1_throat, ff2_throat = IntegrateODE(d, 
        n, N, hbar, r, Del_x, Gamma, Tot_Int_Pts, k, Tot_X_Pts, Shock_Flag, 
        Exit_Pressure, ithroat, a, delta1, rho, InitVal, A, t, U2_in, ff0_throat_in, 
        ff1_throat_in, ff2_throat_in, Mach_E, Mrho_E, Press_E, Temp_E, Vel_E, In_Mass_Flow)

    print(FinalVal)

    return FinalVal
end