using Test
using EegFun
using DataFrames

@testset "DimensionalData Extension" begin
    ext_path = joinpath(dirname(dirname(@__DIR__)), "ext", "EegFunDimensionalDataExt.jl")
    @test isfile(ext_path)

    dim_array = if Base.find_package("DimensionalData") !== nothing
        @eval using DimensionalData
        x -> Base.invokelatest(DimensionalData.DimArray, x)
    else
        # Evaluate the extension against a lightweight sandbox DimensionalData shim when
        # DimensionalData is not installed in the active environment (avoids polluting Base.loaded_modules)
        sandbox = Module(:_DimDataExtSandbox)
        Core.eval(sandbox, Meta.parseall("""
                           module DimensionalData
                           export DimArray, Dim, Ti
                           struct Ti{T}
                               val::T
                           end
                           struct Dim{S,T}
                               val::T
                               Dim{S}(v::T) where {S,T} = new{S,T}(v)
                           end
                           struct DimArray{T,N,A<:AbstractArray{T,N},D<:Tuple}
                               data::A
                               dims::D
                           end
                           dims(da::DimArray) = da.dims
                           end
                           """))
        ext_src = replace(read(ext_path, String), "using DimensionalData" => "using ..DimensionalData")
        Core.eval(sandbox, Meta.parseall(ext_src))
        dd_mod = Base.invokelatest(getproperty, sandbox, :DimensionalData)
        dd_ctor = Base.invokelatest(getproperty, dd_mod, :DimArray)
        x -> Base.invokelatest(dd_ctor, x)
    end

    # 1. ContinuousData & ErpData (SingleDataFrameEeg: Time × Channel)
    cont = EegFun.create_test_continuous_data(n = 10, n_channels = 3)
    da_cont = dim_array(cont)
    @test size(da_cont.data) == (10, 3)

    erp = EegFun.create_test_erp_data(n_channels = 2)
    da_erp = dim_array(erp)
    @test size(da_erp.data) == (EegFun.n_samples(erp), 2)

    # 2. EpochData (MultiDataFrameEeg: Epoch × Time × Channel) + 0-epoch guard
    ep = EegFun.create_test_epoch_data(n_epochs = 4, n_channels = 2)
    da_ep = dim_array(ep)
    @test size(da_ep.data) == (4, EegFun.n_samples(ep), 2)

    empty_ep = EegFun.EpochData("empty", 1, "cond_1", DataFrame[], copy(ep.layout), 250, EegFun.AnalysisInfo())
    da_empty_ep = dim_array(empty_ep)
    @test size(da_empty_ep.data) == (0, 0, 2)

    # 3. SpectrumData (Freq × Channel)
    spec_lay = EegFun.create_test_layout(n_channels = 2)
    spec_df = DataFrame(freq = [2.0, 4.0, 8.0], Ch1 = [1.0, 2.0, 3.0], Ch2 = [4.0, 5.0, 6.0])
    spec = EegFun.SpectrumData("spec", 1, "cond_1", spec_df, spec_lay, 250, :welch, EegFun.AnalysisInfo())
    da_spec = dim_array(spec)
    @test size(da_spec.data) == (3, 2)
    @test da_spec.data[:, 1] == [1.0, 2.0, 3.0]

    # 4. TimeFreqData (Freq × Time × Channel)
    tf = EegFun.create_test_tf_data(n_channels = 2, frequencies = [4.0, 8.0, 12.0], time_points = [0.0, 0.1, 0.2, 0.3])
    da_tf = dim_array(tf)
    @test size(da_tf.data) == (3, 4, 2)
    # First timepoint across all 3 frequencies matches first 3 rows of data_power
    @test da_tf.data[:, 1, 1] == tf.data_power.Ch1[1:3]

    # 5. TimeFreqEpochData (Epoch × Freq × Time × Channel) + 0-epoch guard
    tf_ep =
        EegFun.create_test_tf_epoch_data(n_epochs = 3, n_channels = 2, frequencies = [4.0, 8.0, 12.0], time_points = [0.0, 0.1, 0.2, 0.3])
    da_tf_ep = dim_array(tf_ep)
    @test size(da_tf_ep.data) == (3, 3, 4, 2)
    @test da_tf_ep.data[2, :, 1, 1] == tf_ep.data_power[2].Ch1[1:3]

    empty_tf_ep = EegFun.TimeFreqEpochData(
        "empty_tf",
        1,
        "cond_1",
        DataFrame[],
        DataFrame[],
        copy(tf_ep.layout),
        256,
        :wavelet,
        nothing,
        EegFun.AnalysisInfo(),
    )
    da_empty_tf_ep = dim_array(empty_tf_ep)
    @test size(da_empty_tf_ep.data) == (0, 0, 0, 2)
end
