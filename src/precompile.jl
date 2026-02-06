import PrecompileTools

PrecompileTools.@compile_workload begin
    println("Precompiling EegFun...")

    # TODO: how much does it really help? Initial tests seem to show it helps a bit
    # TODO: what should we put in here? Most common use cases? Everything?

    # Reading *.bdf files and opening databrowser
    # read raw data
    dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
    # read and preprate layout file
    layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
    # create EegFun data structure (EegFun.ContinuousData)
    dat = EegFun.create_eegfun_data(dat, layout_file)
    # Some minimal preprocessing (average reference and highpass filter)
    EegFun.rereference!(dat, :avg)
    EegFun.highpass_filter!(dat, 0.5)
    EegFun.lowpass_filter!(dat, 50)

    # Artifact detection (common operation)
    EegFun.is_extreme_value!(dat, 100)

    # EPOCHS
    epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
    epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.5, 1.0))

    # EPR's
    erp = EegFun.average_epochs(epochs)

    # Key plot types (with display_plot = false to avoid windows)
    EegFun.plot_databrowser(dat)
    EegFun.plot_epochs(epochs[1])
    EegFun.plot_erp(erp)
    EegFun.plot_topography(erp)

end
