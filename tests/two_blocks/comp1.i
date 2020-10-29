# rho  = {rho = 1}
# k    = {k   = 1}
# cp   = {cp  = 1}
# TL   = {TL  = 1000}
# Tmid = {Tmid = 1000}
# init = {initialize = 0}
# {outputfilename = "comp1_output.e"}

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

    Begin ARIA MATERIAL mat1
        Density = constant rho = {rho}
        Thermal Conductivity = constant k = {k}
        Specific Heat = constant cp = {cp}
        Heat Conduction = Fouriers_Law
    End ARIA MATERIAL mat1

    Begin ARIA MATERIAL mat2
        Density = constant rho = 1
        Thermal Conductivity = constant k = 1
        Specific Heat = constant cp = 1
        Heat Conduction = Fouriers_Law
    End ARIA MATERIAL mat2

    Begin finite element model FEModel
        database name = mesh.g
        database type = exodusII

        use material mat1 for block_1
        use material mat2 for block_2
        omit block block_3
    End finite element model FEModel

    {if(initialize==0)}
    Begin finite element model comp2_FEModel
        database name = comp2_output.e
        database type = exodusII
    End finite element model comp2_FEModel
    {endif}

    Begin average value postprocessor myave
        use function solution->temperature
        volumes block_1
    End average value postprocessor myave
    Begin min max postprocessor mymin
        use function solution->temperature
        compute min
        volumes block_1
    End min max postprocessor mymin
    Begin min max postprocessor mymax
        use function solution->temperature
        compute max
        volumes block_1
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
                Begin Sequential solution_block
                    {if(initialize==0)}
                    Advance Input_Region
                    transfer comp2_to_comp1
                    {endif}
                    Advance myRegion
                End Sequential solution_block
            End
        End

        {if(initialize==0)}
        Begin Input_Output region Input_Region
            Use Finite Element Model comp2_FEModel
        END INPUT_OUTPUT REGION Input_Region

        begin transfer comp2_to_comp1
            interpolate surface NODES FROM Input_Region TO myRegion
            SEND BLOCK surface_2 TO surface_2
            SEND field solution->temperature state none TO temperature_from_comp2 state none
            SEARCH GEOMETRIC TOLERANCE IS 0.01
            SEARCH SURFACE GAP TOLERANCE IS 0.1
        end transfer comp2_to_comp1
        {endif}

        Begin Aria region myRegion
            Use Linear Solver linearsolver
            Use Finite Element Model FEModel

            nonlinear solution strategy = newton
            use dof averaged nonlinear residual
            maximum nonlinear iterations = 8
            nonlinear residual tolerance = 1.0e-8
            nonlinear relaxation factor       = 1.0


            EQ Energy for Temperature on all_blocks using Q1 with Diff
            IC for Temperature ON all_blocks = constant value = 300.            

            Begin Temperature Boundary Condition surface_3
                add surface surface_3
                temperature = {TL}
            End Temperature Boundary Condition surface_3

            {if(initialize==0)}
            Begin User Variable temperature_from_comp2
                Type = Node real length=1
                Add Part surface_2
            End

            Begin Temperature Boundary Condition interface
                add surface surface_2
                temperature node variable = temperature_from_comp2
            End
            {else}
            Begin Temperature Boundary Condition surface_2
                add surface surface_2
                temperature = {Tmid}
            End Temperature Boundary Condition surface_2
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
                Include = block_2 #surface_1
                Property compose_results = yes
                at step 1, increment= 1

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

