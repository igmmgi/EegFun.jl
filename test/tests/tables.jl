using Test
using EegFun
using DataFrames
import Tables

@testset "Tables.jl Interface" begin
    # Create test data
    layout_file = joinpath(pkgdir(EegFun), "resources", "layouts", "biosemi", "biosemi72.csv")
    data_file = joinpath(pkgdir(EegFun), "resources", "data", "bdf", "example1.bdf")

    dat = EegFun.create_eegfun_data(EegFun.read_raw_data(data_file), EegFun.read_layout(layout_file))

    @testset "ContinuousData" begin
        @test Tables.istable(dat)
        @test Tables.istable(typeof(dat))
        @test Tables.columnaccess(typeof(dat))

        df_cont = DataFrame(dat)
        @test "file" in names(df_cont)
        @test "time" in names(df_cont)
        @test df_cont.file[1] == dat.file
        @test nrow(df_cont) == nrow(dat.data)
    end

    @testset "EpochData & ErpData" begin
        epoch_cfg = [EegFun.EpochCondition(name = "Test", trigger_sequences = [[1]])]
        epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.1, 0.2))

        @test Tables.istable(epochs[1])
        df_epoch = DataFrame(epochs[1])

        @test "file" in names(df_epoch)
        @test "condition" in names(df_epoch)
        @test "condition_name" in names(df_epoch)
        @test "epoch" in names(df_epoch) # verified that epoch is present in df internally
        @test df_epoch.condition[1] == epochs[1].condition
        @test df_epoch.condition_name[1] == epochs[1].condition_name
        @test nrow(df_epoch) == sum(nrow.(epochs[1].data))

        erp = EegFun.average_epochs(epochs)

        @test Tables.istable(erp[1])
        df_erp = DataFrame(erp[1])

        @test "file" in names(df_erp)
        @test "condition" in names(df_erp)
        @test "condition_name" in names(df_erp)
        @test "n_epochs" in names(df_erp)
        @test df_erp.n_epochs[1] == erp[1].n_epochs
        @test nrow(df_erp) == nrow(erp[1].data)

        # Test Vector{<:EegData}
        @test Tables.istable(epochs)
        df_epochs_vec = DataFrame(epochs)
        @test "condition" in names(df_epochs_vec)
        # Should be equal to the sum of rows of all conditions
        @test nrow(df_epochs_vec) == sum(nrow.(DataFrame.(epochs)))

        @test Tables.istable(erp)
        df_erp_vec = DataFrame(erp)
        @test nrow(df_erp_vec) == sum(nrow.(DataFrame.(erp)))
    end

    @testset "Tables.rowaccess and DataAPI forwarding" begin
        @test Tables.rowaccess(dat)
        @test Tables.rowaccess(typeof(dat))
        rows = collect(Tables.rows(dat))
        @test length(rows) == nrow(dat.data)
        @test Tables.schema(dat) !== nothing

        # Test first and last
        @test nrow(first(dat, 5)) == 5
        @test nrow(last(dat, 3)) == 3
        @test first(dat) isa DataFrameRow
        @test last(dat) isa DataFrameRow

        # Test eachcol and eachrow
        @test length(eachcol(dat)) == ncol(dat.data)
        @test length(eachrow(dat)) == nrow(dat.data)

        # Test convert and parent
        @test convert(DataFrame, dat) isa DataFrame
        @test nrow(convert(DataFrame, dat)) == nrow(dat.data)
        @test parent(dat) === dat.data

        # Test names with cols selector
        @test names(dat, r"^time") == ["time"]
    end
end
