import PrecompileTools

PrecompileTools.@compile_workload begin

    # Get the package directory to construct absolute paths
    pkg_dir = pkgdir(EegFun)
    data_file = joinpath(pkg_dir, "resources", "data", "bdf", "example1.bdf")
    layout_file_path = joinpath(pkg_dir, "resources", "layouts", "biosemi", "biosemi72.csv")

    if isfile(data_file) && isfile(layout_file_path)
        println("Precompiling EegFun...")

        # 1. Reading *.bdf files and opening databrowser
        dat = EegFun.read_raw_data(data_file)
        layout_file = EegFun.read_layout(layout_file_path)
        dat = EegFun.create_eegfun_data(dat, layout_file)

        # 2. Common preprocessing pipeline
        EegFun.rereference!(dat, :avg)
        EegFun.highpass_filter!(dat, 0.5)
        EegFun.lowpass_filter!(dat, 50)

        # 3. Artifact detection
        EegFun.is_extreme_value!(dat, 100)

        # 4. Epoch extraction & ERP workflow
        epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
        epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.5, 1.0))
        erp = EegFun.average_epochs(epochs)

        # 5. Key plot types (headless — display_plot = false prevents window creation)
        EegFun.plot_databrowser(dat; display_plot = false)
        EegFun.plot_epochs(epochs[1]; display_plot = false)
        EegFun.plot_erp(erp; display_plot = false)
        EegFun.plot_topography(erp; display_plot = false)

        # 6. Predicate generators 
        EegFun.channels()
        EegFun.channels(:Cz)
        EegFun.channels(:Fz, :Cz, :Pz)
        EegFun.channels(1:10)
        EegFun.channels_not(:Fp1)
        EegFun.epochs()
        EegFun.epochs(1:5)
        EegFun.conditions()
        EegFun.conditions(1)
        EegFun.participants()
        EegFun.samples()
        EegFun.times()
        EegFun.times(-0.2, 0.5)
        EegFun.components()

        # 7. Data access utilities 
        EegFun.channel_labels(dat)
        EegFun.channel_labels(erp[1])
        EegFun.channel_labels(epochs[1])
        EegFun.meta_labels(dat)
        EegFun.all_data(dat)

        # 8. Subset operations
        EegFun.subset(erp[1]; channel_selection = EegFun.channels(:Cz), interval_selection = EegFun.times(-0.2, 0.5))
        EegFun.subset(erp; condition_selection = EegFun.conditions(1))

        println("Precompilation complete!")
    else
        println("Skipping precompilation workload (example data files not found)")
        println("  Expected data file: $data_file")
        println("  Expected layout file: $layout_file_path")
    end

end
