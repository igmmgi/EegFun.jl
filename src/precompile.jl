import PrecompileTools

PrecompileTools.@compile_workload begin

    # Get the package directory to construct absolute paths
    pkg_dir = pkgdir(EegFun)
    data_file = joinpath(pkg_dir, "resources", "data", "bdf", "example1.bdf")
    layout_file_path = joinpath(pkg_dir, "resources", "layouts", "biosemi", "biosemi72.csv")

    if isfile(data_file) && isfile(layout_file_path)
        println("Precompiling EegFun...")

        # TODO: how much does it really help? Initial tests seem to show it helps a bit
        # TODO: what should we put in here? Most common use cases? Everything?
        # TODO: GitHub actions does not see to like the plot_xxx calls?

        # Reading *.bdf files and opening databrowser
        dat = EegFun.read_raw_data(data_file)
        layout_file = EegFun.read_layout(layout_file_path)
        dat = EegFun.create_eegfun_data(dat, layout_file)

        # Some minimal preprocessing
        EegFun.rereference!(dat, :avg)
        EegFun.highpass_filter!(dat, 0.5)
        EegFun.lowpass_filter!(dat, 50)

        # Artifact detection
        EegFun.is_extreme_value!(dat, 100)

        # Epoch extraction & ERP workflow
        epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
        local epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.5, 1.0))
        local erp = EegFun.average_epochs(epochs)

        # # Key plot types
        # EegFun.plot_databrowser(dat)
        # EegFun.plot_epochs(epochs[1])
        # EegFun.plot_erp(erp)
        # EegFun.plot_topography(erp)

        println("Precompilation complete!")
    else
        println("Skipping precompilation workload (example data files not found)")
        println("  Expected data file: $data_file")
        println("  Expected layout file: $layout_file_path")
    end

end
