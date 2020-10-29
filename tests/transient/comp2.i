# rho  = {rho = 100}
# k    = {k   = 1}
# cp   = {cp  = 1}
# TR   = {TR  = 500}
# init = {initialize = 0}
# {outputfilename = "comp2_output.e"}

Begin sierra mySquare
    Title Simple Conduction Problem

    BEGIN TPETRA EQUATION SOLVER linearsolver
        BEGIN CG SOLVER
            BEGIN JACOBI PRECONDITIONER
            END

            MAXIMUM ITERATIONS = 1000
            RESIDUAL SCALING = NONE
            CONVERGENCE TOLERANCE = 1.000000e-12
        END
    END TPETRA EQUATION SOLVER

    Begin ARIA MATERIAL mat3
        Density = constant rho = {rho}
        Thermal Conductivity = constant k = {k}
        Specific Heat = constant cp = {cp}
        Heat Conduction = Fouriers_Law
    End ARIA MATERIAL mat3

    Begin ARIA MATERIAL mat2
        Density = constant rho = 100
        Thermal Conductivity = constant k = 1
        Specific Heat = constant cp = 1
        Heat Conduction = Fouriers_Law
    End ARIA MATERIAL mat2

    Begin finite element model FEModel
        database name = mesh.g
        database type = exodusII

        use material mat2 for block_2
        use material mat3 for block_3
        omit block block_1
    End finite element model FEModel

    {if(initialize==0)}
    Begin finite element model comp1_FEModel
        database name = comp1_output.e
        database type = exodusII
    End finite element model comp1_FEModel
    {endif}

    Begin average value postprocessor myave
        use function solution->temperature
        volumes block_3
    End average value postprocessor myave
    Begin min max postprocessor mymin
        use function solution->temperature
        compute min
        volumes block_3
    End min max postprocessor mymin
    Begin min max postprocessor mymax
        use function solution->temperature
        compute max
        volumes block_3
    End min max postprocessor mymax

    Begin average value postprocessor overlapave
        use function solution->temperature
        volumes block_2
    End average value postprocessor overlapave
    Begin min max postprocessor overlapmin
        use function solution->temperature
        compute min
        volumes block_2
    End min max postprocessor overlapmin
    Begin min max postprocessor overlapmax
        use function solution->temperature
        compute max
        volumes block_2
    End min max postprocessor overlapmax


    Begin PROCEDURE Aria_procedure
        Begin Solution Control Description
            Use System main
            Begin System main
                Simulation Start Time = 0.0
                Begin Transient TimeBlock
                    {if(initialize==0)}
                    Advance Input_Region
                    transfer comp1_to_comp2
                    {endif}
                    Advance myRegion
                End Transient TimeBlock
            End
            Begin Parameters For Transient TimeBlock
                Start Time = 0.0
                Termination Time = 4
                Begin Parameters For Aria Region myRegion
                    Time Step Variation = Adaptive
                    Time Integration Method = First_Order
                    Initial time step size  = 1
                    Maximum Time Step Size Ratio  = 2
                    minimum resolved time step size = 1.e-6
                    minimum time step size = 1.e-6
                    maximum time step size = 20.0
                    Predictor Order = 1
                    Predictor-corrector tolerance = 1e-3
                    Predictor-Corrector Begin After Step = 4
                End
            End Parameters For Transient TimeBlock
        End

        {if(initialize==0)}
        Begin Input_Output region Input_Region
            Use Finite Element Model comp1_FEModel
        END INPUT_OUTPUT REGION Input_Region

        begin transfer comp1_to_comp2
            interpolate surface NODES FROM Input_Region TO myRegion
            SEND BLOCK surface_1 TO surface_1
            SEND field solution->temperature state none TO temperature_from_comp1 state none
            SEARCH GEOMETRIC TOLERANCE IS 0.01
            SEARCH SURFACE GAP TOLERANCE IS 0.1
        end transfer comp1_to_comp2
        {endif}

        Begin Aria region myRegion
            Use Linear Solver linearsolver
            Use Finite Element Model FEModel

            nonlinear solution strategy = newton
            use dof averaged nonlinear residual
            maximum nonlinear iterations = 8
            nonlinear residual tolerance = 1.0e-8
            nonlinear relaxation factor       = 1.0


            EQ Energy for Temperature on all_blocks using Q1 with Diff lumped_mass src
            IC for Temperature ON all_blocks = constant value = 300.            

            Begin Temperature Boundary Condition surface_4
                add surface surface_4
                temperature = {TR}
            End Temperature Boundary Condition surface_4

            {if(initialize==0)}
            Begin User Variable temperature_from_comp1
                Type = Node real length=1
                Add Part surface_1
            End

            Begin Temperature Boundary Condition interface
                add surface surface_1
                temperature node variable = temperature_from_comp1
            End
            {else}
            Begin Temperature Boundary Condition interface
                add surface surface_1
                temperature = 300
            End
            {endif}

            evaluate postprocessor myave
            evaluate postprocessor mymin
            evaluate postprocessor mymax
            evaluate postprocessor overlapave
            evaluate postprocessor overlapmin
            evaluate postprocessor overlapmax
            Begin Data Probe myprobe
                Nodal_Location solution->temperature Location is 0. 0. 0. Label centerpoint
            End


            Begin Results Output data
                Database Name = {outputfilename}
                #Include = block_2
                Property compose_results = yes
                #at time 0, increment = 5
                #at time 50, increment = 10
                #timestep adjustment interval = 2
                at step 0, increment = 1

                Title test results
                nodal Variables = solution->temperature

                global variables = myave
                global variables = mymin
                global variables = mymax
                global variables = overlapave
                global variables = overlapmin
                global variables = overlapmax
                global variables = centerpoint
            End Results Output data
        End Aria region myRegion
    End PROCEDURE Aria_procedure
End sierra mySquare

