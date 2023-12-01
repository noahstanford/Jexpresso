include("../InitCalcParms.jl");
include("../IntegrateODE.jl");
include("./WriteResults.jl");

function QNavStokes_Solvr( a, Tot_TSteps, in_n, d, r,err1, delta, rho, x_min, x_max, Tot_X_Pts, Gamma, ICMFlowErrScale, ICrhoErrScale, ICtempErrScale, In_Mass_Flow, Exit_Pressure, Shock_Flag )
   
    #QNAVSTOKES_SOLVR solves ODEs for discretized 1D Navier-Stokes flow.
    #
    #   QNavStokes_Solvr solves the coupled set of non-linear ODEs 
    #   resulting from discretizing the Navier_Stokes dynamics for a quasi-1D 
    #   nozzle flow. It returns an approximate solution z[t] which satisfies  
    #   an error upper bound with probability 1 - delta. All physical 
    #   quantities scaled appropriately & so are dimensionless.
    #
    #   This algorithm is based on quantum algorithm for solving ODEs 
    #   due to B. Kacewicz, J. Complexity vol. 22 [2006] 676-690. 
    #
    #   INPUTS: 
    #       a = initial time [final time b = Tot_TSteps*Delta_t]
    #       Tot_TSteps = total number of time-steps taken by solver
    #       in_n = initial number of time subintervals in primary partition 
    #               of [a,b]
    #       d = number of components in ODE solution z[t] at an x-location 
    #       r = number of derivatives allowed for ODE driver function f(z) 
    #       err1 = upper bound on error of approximate value found for an 
    #               integral using quantum integration algorithm 
    #       delta = probability that approximate solution z[t] violates its 
    #                   error upper bound bound [Eq. (13) in above paper] 
    #       rho = basic Holder class parameter; must be estimated for each
    #               application - tentatively set to 1 in nozzle flow problem
    #       x_min = x-location of nozzle entrance
    #       x_max = x-location of nozzle exit()
    #       Tot_X_Pts = total number of x-grid points in nozzle
    #       Gamma = ratio of specific heats c_p/c_v
    #       ICMFlowErrScale = scale of random shift to initial mass flow rate 
    #                           away from exact result
    #       ICrhoErrScale = scale of random shift to initial mass density 
    #                           away from exact result
    #       ICtempErrScale = scale of random shift to initial temperature 
    #                           away from exact result
    #       In_Mass_Flow = initial value for mass flow through nozzle
    #       Exit_Pressure = pressure an nozzle exit()
    #       Shock_Flag = 0 [1] if shock wave absent [present]
    #
    #   Support functions: InitCalcParms; IntegrateODE; WriteResults
    # call InitCalcParms to initialize input parameters for quantum
    #   Navier-Stokes algorithm

    #addpath("./")

    Del_x, A, Mach_E, Mrho_E, Press_E, Temp_E, Vel_E, InitVal, n, N, t, hbar, Tot_Int_Pts, U2_in, delta1, k, ithroat, ff0_throat_in, ff1_throat_in, ff2_throat_in = InitCalcParms(x_min, x_max, Tot_X_Pts, Gamma, Exit_Pressure, Shock_Flag, In_Mass_Flow, d, ICMFlowErrScale, ICrhoErrScale, ICtempErrScale, in_n, delta, err1, Tot_TSteps, a)

    t = @elapsed begin
    
    U2, Mach_D, Mrho_D, Press_D, Temp_D, Vel_D, Rel_MachErr, Rel_MrhoErr, Rel_PressErr, Rel_TempErr, Rel_VelErr, AvRelTempErr, AvPlusSDevRelTempErr, AvMinusSDevRelTempErr, AvRelMachErr, AvPlusSDevRelMachErr, AvMinusSDevRelMachErr, AvRelMrhoErr, AvPlusSDevRelMrhoErr, AvMinusSDevRelMrhoErr, AvRelPressErr, AvPlusSDevRelPressErr, AvMinusSDevRelPressErr, AvU2, ff0_throat, ff1_throat, ff2_throat = IntegrateODE(d, n, N, hbar, r, Del_x, Gamma, Tot_Int_Pts, k, Tot_X_Pts, Shock_Flag, Exit_Pressure, ithroat, a, delta1, rho, InitVal, A, t, U2_in, ff0_throat_in, ff1_throat_in, ff2_throat_in, Mach_E, Mrho_E, Press_E, Temp_E, Vel_E, In_Mass_Flow)

    end

    runtime = t/60
    timepersubint = runtime/n

    # write computational results to files for eventual plotting:

    WriteResults(n, Tot_X_Pts, d, U2, Mach_D, Mach_E, Mrho_D, Mrho_E, Press_D, Press_E, Temp_D, Temp_E, Rel_MachErr, Rel_MrhoErr, Rel_PressErr, Rel_TempErr, AvRelTempErr, AvPlusSDevRelTempErr, AvMinusSDevRelTempErr, AvRelMachErr, AvPlusSDevRelMachErr, AvMinusSDevRelMachErr, AvRelMrhoErr, AvPlusSDevRelMrhoErr, AvMinusSDevRelMrhoErr, AvRelPressErr, AvPlusSDevRelPressErr, AvMinusSDevRelPressErr, AvU2, ff0_throat, ff1_throat, ff2_throat);

    message = "QNavStokes_solvr has finished; results written to files.";

    print(message)

    print("Program runtime (minutes) = ", string(runtime))

    print("Program runtime per subinterval(minutes) = ", string(timepersubint))
end