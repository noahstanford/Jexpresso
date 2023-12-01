using Printf

function WriteResults(n, Tot_X_Pts, d, U2, Mach_D, Mach_E, Mrho_D, Mrho_E, Press_D, Press_E, Temp_D, Temp_E, Rel_MachErr, Rel_MrhoErr, Rel_PressErr, Rel_TempErr, AvRelTempErr, AvPlusSDevRelTempErr, AvMinusSDevRelTempErr, AvRelMachErr, AvPlusSDevRelMachErr, AvMinusSDevRelMachErr, AvRelMrhoErr, AvPlusSDevRelMrhoErr, AvMinusSDevRelMrhoErr, AvRelPressErr, AvPlusSDevRelPressErr, AvMinusSDevRelPressErr, AvU2, ff0_throat, ff1_throat, ff2_throat)
    # @printfRESULTS @printfs results of Navier-Stokes solution to files.
    #   @printfResults @printfs primary flow variables and related results to files
    #       for eventual plotting.
    #
    #   INPUTS:
    #       n = number of subintervals used in time partition
    #       Tot_X_Pts = total number of spatial grid-points used
    #       d = number of components of ODE solution z(t)
    #       U2 = Tot_X_Pts x (n+1) array that stores calculated mass flow rate
    #               at beginning of each time subinterval, and the mass flow
    #               rate at the final time of the ODE integration.
    #       Mach_D = 1 x Tot_X_Pts array storing calculated Mach number at each
    #                   grid-point at end of ODE integration
    #       Mach_E = 1 x Tot_X_Pts array storing exact steady-state Mach 
    #                   number at each grid-point at end of ODE integration
    #       Mrho_D = 1 x Tot_X_Pts array storing calculated (dimensionless)
    #                   mass density at each grid-point at end of ODE 
    #                   integration
    #       Mrho_E = 1 x Tot_X_Pts array storing exact (dimensionless) mass 
    #                   density at each grid-point at end of ODE integration
    #       Press_D = 1 x Tot_X_Pts array storing calculated (dimensionless)
    #                   pressure at each grid-point at end of ODE integration
    #       Press_E = 1 x Tot_X_Pts array storing exact (dimensionless)
    #                   pressure at each grid-point at end of ODE integration
    #       Temp_D = 1 x Tot_X_Pts array storing calculated (dimensionless)
    #                   temperature at each grid-point at end of ODE 
    #                   integration
    #       Temp_E = 1 x Tot_X_Pts array storing exact (dimensionless)
    #                   temperature at each grid-point at end of ODE 
    #                   integration
    #       Rel_MachErr = 1 x Tot_X_Pts array storing relative error in Mach
    #                       number at end of ODE integration
    #       Rel_MrhoErr = 1 x Tot_X_Pts array storing relative error in 
    #                       dimensionless mass density at end of ODE 
    #                       integration
    #       Rel_PressErr = 1 x Tot_X_Pts array storing relative error in 
    #                       dimensionless pressure at end of ODE 
    #                       integration
    #       Rel_TempErr = 1 x Tot_X_Pts array storing relative error in 
    #                       dimensionless temperature at end of ODE 
    #                       integration
    #       AvRelTempErr = 1 x Tot_X_Pts array storing average relative
    #                       temperature error, same value to each grid-point
    #       AvPlusSDevRelTempErr = 1 x Tot_X_Pts array storing (average plus
    #                       standard deviation) relative temperature error, 
    #                       same value to each grid-point
    #       AvMinusSDevRelTempErr = 1 x Tot_X_Pts array storing (average minus
    #                       standard deviation) relative temperature error, 
    #                       same value to each grid-point
    #       AvRelMachErr = 1 x Tot_X_Pts array storing average relative Mach
    #                       number error, same value to each grid-point
    #       AvPlusSDevRelMachErr = 1 x Tot_X_Pts array storing (average plus
    #                       standard deviation) relative Mach number error, 
    #                       same value to each grid-point
    #       AvMinusSDevRelMachErr = 1 x Tot_X_Pts array storing (average minus
    #                       standard deviation) relative Mach number error, 
    #                       same value to each grid-point
    #       AvRelMrhoErr = 1 x Tot_X_Pts array storing average relative mass
    #                       density error, same value to each grid-point
    #       AvPlusSDevRelMrhoErr = 1 x Tot_X_Pts array storing (average plus
    #                       standard deviation) relative mass density error, 
    #                       same value to each grid-point
    #       AvMinusSDevRelMrhoErr = 1 x Tot_X_Pts array storing (average minus
    #                       standard deviation) relative mass density error, 
    #                       same value to each grid-point
    #       AvRelPressErr = 1 x Tot_X_Pts array storing average relative 
    #                       pressure error, same value to each grid-point
    #       AvPlusSDevRelPressErr = 1 x Tot_X_Pts array storing (average plus
    #                       standard deviation) relative pressure error, 
    #                       same value to each grid-point
    #       AvMinusSDevRelPressErr = 1 x Tot_X_Pts array storing (average minus
    #                       standard deviation) relative pressure error, 
    #                       same value to each grid-point
    #       AvU2 = 1 x Tot_X_Pts array storing average mass flow rate (over 
    #                 nozzle) at final time of ODE integration at each grid-pt
    #       ff0_throat = d x n array storing driver function value at throat at
    #                       end of each subinterval
    #       ff1_throat = d x n array storing first derivative of driver 
    #                       function at throat at end of each subinterval
    #       ff2_throat = d x n array storing second derivative of driver 
    #                       function at throat at end of each subinterval

    # @printf results for various flow variables to file:

    filenameU2 = open("./Output/U2vals.txt", "w+");
    filenameAvU2 = open("./Output/AvU2vals.txt", "w+");
    filenameMachD = open("./Output/MachDvals.txt", "w+");
    filenameMachE = open("./Output/MachEvals.txt", "w+");
    filenameMrhoD = open("./Output/MrhoDvals.txt", "w+");
    filenameMrhoE = open("./Output/MrhoEvals.txt", "w+");
    filenamePressD = open("./Output/PressDvals.txt", "w+");
    filenamePressE = open("./Output/PressEvals.txt", "w+");
    filenameTempD = open("./Output/TempDvals.txt", "w+");
    filenameTempE = open("./Output/TempEvals.txt", "w+");

    column = n+1;
    row = Tot_X_Pts;

    for gridpt = 1:row
        @printf(filenameAvU2, "%8.3f", AvU2[Int(gridpt)]);
        
        if gridpt == 1
            for tyme = 1:column
                @printf(filenameU2, "%8.3f", U2[Int(gridpt), Int(tyme)]);
            end
        elseif gridpt != 1
            @printf(filenameU2, "\n");
            
            for tyme = 1:column
                @printf(filenameU2, "%8.3f", U2[Int(gridpt), Int(tyme)]);
            end
        end
    end

    for gridpt = 1:Tot_X_Pts
        gridpt = Int(gridpt);
        @printf(filenameMachD, "%6.3f", Mach_D[gridpt]);
        @printf(filenameMachE, "%6.3f", Mach_E[gridpt]);
        @printf(filenameMrhoD, "%6.3f", Mrho_D[gridpt]);
        @printf(filenameMrhoE, "%6.3f", Mrho_E[gridpt]);
        @printf(filenamePressD, "%6.3f", Press_D[gridpt]);
        @printf(filenamePressE, "%6.3f", Press_E[gridpt]);
        @printf(filenameTempD, "%6.3f", Temp_D[gridpt]);
        @printf(filenameTempE, "%6.3f", Temp_E[gridpt]);
    end

    close(filenameU2);
    close(filenameAvU2);
    close(filenameMachD);
    close(filenameMachE);
    close(filenameMrhoD);
    close(filenameMrhoE);
    close(filenamePressD);
    close(filenamePressE);
    close(filenameTempD);
    close(filenameTempE);


    # open files for writing flow variable relative errors
    filenameRelMachErr = open("./Output/RelMachErrvals.txt", "w+");
    filenameRelMrhoErr = open("./Output/RelMrhoErrvals.txt", "w+");
    filenameRelPressErr = open("./Output/RelPressErrvals.txt", "w+");
    filenameRelTempErr = open("./Output/RelTempErrvals.txt", "w+");

    # @printf data to files
    for gridpt = 1:Tot_X_Pts
        gridpt = Int(gridpt);
        @printf(filenameRelMachErr, "%6.3f", Rel_MachErr[gridpt]);
        @printf(filenameRelMrhoErr, "%6.3f", Rel_MrhoErr[gridpt]);
        @printf(filenameRelPressErr, "%6.3f", Rel_PressErr[gridpt]);
        @printf(filenameRelTempErr, "%6.3f", Rel_TempErr[gridpt]);
    end

    # close files 
    close(filenameRelMachErr);
    close(filenameRelMrhoErr);
    close(filenameRelPressErr);
    close(filenameRelTempErr);

    # open files to @printf mean and mean +/- standard deviation of relative
    #   errors
    filenameAvRelTempErr = open("./Output/AvRelTempErr.txt", "w+");
    filenameAvRelMachErr = open("./Output/AvRelMachErr.txt", "w+");
    filenameAvRelMrhoErr = open("./Output/AvRelMachErrAvRelMrhoErr.txt", "w+");
    filenameAvRelPressErr = open("./Output/AvRelMachErrAvRelPressErr.txt", "w+");
    filenameAvPlusSDevRelTempErr = open("./Output/AvRelMachErrAvPlusSDevRelTempErr.txt", "w+");
    filenameAvPlusSDevRelMachErr = open("./Output/AvRelMachErrAvPlusSDevRelMachErr.txt", "w+");
    filenameAvPlusSDevRelMrhoErr = open("./Output/AvRelMachErrAvPlusSDevRelMrhoErr.txt", "w+");
    filenameAvPlusSDevRelPressErr = open("./Output/AvRelMachErrAvPlusSDevRelPressErr.txt", "w+");
    filenameAvMinusSDevRelTempErr = open("./Output/AvRelMachErrAvMinusSDevRelTempErr.txt", "w+");
    filenameAvMinusSDevRelMachErr = open("./Output/AvRelMachErrAvMinusSDevRelMachErr.txt", "w+");
    filenameAvMinusSDevRelMrhoErr = open("./Output/AvRelMachErrAvMinusSDevRelMrhoErr.txt", "w+");
    filenameAvMinusSDevRelPressErr = open("./Output/AvRelMachErrAvMinusSDevRelPressErr.txt", "w+");

    for gridpt = 1:Tot_X_Pts
        gridpt = Int(gridpt);
        @printf(filenameAvRelTempErr, "%8.6f", AvRelTempErr[gridpt]);
        @printf(filenameAvRelMachErr, "%8.6f", AvRelMachErr[gridpt]);
        @printf(filenameAvRelMrhoErr, "%8.6f", AvRelMrhoErr[gridpt]);
        @printf(filenameAvRelPressErr, "%8.6f", AvRelPressErr[gridpt]);
        @printf(filenameAvPlusSDevRelTempErr, "%8.6f", AvPlusSDevRelTempErr[gridpt]);
        @printf(filenameAvPlusSDevRelMachErr, "%8.6f", AvPlusSDevRelMachErr[gridpt]);
        @printf(filenameAvPlusSDevRelMrhoErr, "%8.6f", AvPlusSDevRelMrhoErr[gridpt]);
        @printf(filenameAvPlusSDevRelPressErr, "%8.6f", AvPlusSDevRelPressErr[gridpt]);
        @printf(filenameAvMinusSDevRelTempErr, "%8.6f", AvMinusSDevRelTempErr[gridpt]);
        @printf(filenameAvMinusSDevRelMachErr, "%8.6f", AvMinusSDevRelMachErr[gridpt]);
        @printf(filenameAvMinusSDevRelMrhoErr, "%8.6f", AvMinusSDevRelMrhoErr[gridpt]);
        @printf(filenameAvMinusSDevRelPressErr, "%8.6f", AvMinusSDevRelPressErr[gridpt]);
    end

    # close files
    close(filenameAvRelTempErr);
    close(filenameAvRelMachErr);
    close(filenameAvRelMrhoErr);
    close(filenameAvRelPressErr);
    close(filenameAvPlusSDevRelTempErr);
    close(filenameAvPlusSDevRelMachErr);
    close(filenameAvPlusSDevRelMrhoErr);
    close(filenameAvPlusSDevRelPressErr);
    close(filenameAvMinusSDevRelTempErr);
    close(filenameAvMinusSDevRelMachErr);
    close(filenameAvMinusSDevRelMrhoErr);
    close(filenameAvMinusSDevRelPressErr);

    # open files to @printf residual and its first r time derivatives at throat
    filenameff0 = open("./Output/ff0Throatvals.txt","w+");
    filenameff1 = open("./Output/ff1Throatvals.txt","w+");
    filenameff2 = open("./Output/ff2Throatvals.txt","w+");

    # @printf values to files
    for varble = 1:d
        if varble == 1
            for subint = 1:n
                @printf(filenameff0, "%6.3f", ff0_throat[Int(varble), Int(subint)]);
            end
        elseif varble != 1
            @printf(filenameff0, "\n");
            
            for subint = 1:n
                @printf(filenameff0, "%6.3f",  ff0_throat[Int(varble), Int(subint)]);
            end
        end
    end

    for varble = 1:d
        if varble == 1
            for subint = 1:n
                @printf(filenameff1, "%6.3f", ff1_throat[Int(varble), Int(subint)]);
            end
        elseif varble != 1
            @printf(filenameff1, "\n");
            
            for subint = 1:n
                @printf(filenameff1, "%6.3f",  ff1_throat[Int(varble), Int(subint)]);
            end
        end
    end

    for varble = 1:d
        if varble == 1
            for subint = 1:n
                @printf(filenameff2, "%6.3f", ff2_throat[Int(varble), Int(subint)]);
            end
        elseif varble != 1
            @printf(filenameff2, "\n");
            
            for subint = 1:n
                @printf(filenameff2, "%6.3f",  ff2_throat[Int(varble), Int(subint)]);
            end
        end
    end

    # close files and exit
    close(filenameff0);
    close(filenameff1);
    close(filenameff2);
end