function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # User define your inputs below: the order doesn't matter
        # IMPORTANT NOTICE: DO NOT FORGET the "," at the end of each entry!!!
        #---------------------------------------------------------------------------
        :tend                 => 3.0,
        :ode_solver          => "Tsit5", #"SSPRK53",
        :Δt                   => 1e-3,
        :ndiagnostics_outputs => 20, #these are steps, not seconds
        #:output_dir          => "/Users/simone/runs/",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl", # Choice: "lgl", "cg", "cgl"
        :nop                 => 4,     # Polynomial order
        :lexact_integration  => false,
        :luser_bc            => true,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :case                 => "kopriva.5.3.5",
        :lvisc                => true, #true is defualt
        :νx                   => 0.01, #kinematic viscosity constant
        :νy                   => 0.01, #kinematic viscosity constant
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/2d-grid.msh", 
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_25x25.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/circle_TFI.msh",
        #:gmsh_filename        => "./meshes/gmsh_grids/hexa_TFI_10x10_periodic.msh",
        #---------------------------------------------------------------------------
        # Output formats: "png" -> plots to png file. "ascii" -> data to npoin file
        #---------------------------------------------------------------------------
        :outformat     => "png", #choice: "png", "ascii" (default is ascii)
        #---------------------------------------------------------------------------
        # 1D (lread_gmsh => faluse): the grid is built by jexpresso
        #---------------------------------------------------------------------------
        :xmin          =>   0.0,
        :xmax          =>   1.0,
        :nelx          =>   25,        
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
