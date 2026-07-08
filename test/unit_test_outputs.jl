@testset "outputs" verbose=true begin

    ac = load_default_model()
    size_aircraft!(ac; printiter=false);

    file_weights = joinpath(TASOPT.__TASOPTroot__, "../test/weights.txt")
    file_aero = joinpath(TASOPT.__TASOPTroot__, "../test/aero.txt")
    file_geom = joinpath(TASOPT.__TASOPTroot__, "../test/geom.txt")
    
    @testset "output summaries" begin
        f = open(io->TASOPT.weight_buildup(ac,io=io), "temp.txt", "w")
        content = read(file_weights, String)
        content = replace(content, "\r\n"=> "\n") # Convert \r\n to \n
        @test content == read("temp.txt", String)
        rm("temp.txt")

        f = open(io->TASOPT.aero(ac,io=io), "temp.txt", "w")
        content = read(file_aero, String)
        content = replace(content, "\r\n"=> "\n") # Convert \r\n to \n
        @test content == read("temp.txt", String)
        rm("temp.txt")

        f = open(io->TASOPT.geometry(ac,io=io), "temp.txt", "w")
        content = read(file_geom, String)
        content = replace(content, "\r\n"=> "\n") # Convert \r\n to \n
        @test content == read("temp.txt", String)
        rm("temp.txt")
    end

    @testset "output plots" begin
        
        p1 = TASOPT.stickfig(ac)
        @test p1 isa Plots.Plot
        p2 = TASOPT.plot_details(ac)
        @test p2 isa Plots.Plot
        p3 = TASOPT.plot737compare(ac)
        @test p3 isa Plots.Plot
        p4 = TASOPT.MomentShear(ac)
        @test p4 isa Plots.Plot
        pr_state = (
            parg = copy(ac.parg),
            parm = copy(ac.parm),
            para = copy(ac.para),
            pare = copy(ac.pare),
        )
        pr = TASOPT.payload_range_data(ac; Rpts=5, Ppts=5, printTO=false)

        @test !isempty(pr.range)
        @test issorted(pr.range)
        @test all(diff(pr.range) .> 0)
        @test all(diff(pr.payload) .<= 1e-8)
        @test ac.parm[imRange, 1] in pr.range
        @test all(pr.payload .<= ac.parg[igWpaymax] .+ 1.0)
        @test all(pr.fuel .<= pr.Fuelmax .+ 1.0)
        @test all(pr.WTO .<= pr.Wmax .+ 1.0)
        @test all(isfinite.(pr.pfei[pr.payload .> 0.0]))
        @test isequal(ac.parg, pr_state.parg)
        @test isequal(ac.parm, pr_state.parm)
        @test isequal(ac.para, pr_state.para)
        @test isequal(ac.pare, pr_state.pare)

        p5 = TASOPT.PayloadRange(ac; Rpts=5, Ppts=5, printTO=false)
        @test p5 isa Plots.Plot
        p6 = DragPolar(ac)
        @test p6 isa Plots.Plot
        p7 = plot_airf(ac)
        @test p7 isa Plots.Plot
        p8 = plot_drag_breakdown(ac, show_values=true)
        p9 = plot_drag_breakdown(ac, ip = :, show_fractions=false)
        @test p8 isa Plots.Plot
        @test p9 isa Plots.Plot
    end

end
