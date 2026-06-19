import PrecompileTools

PrecompileTools.@compile_workload begin

    # Get the package directory to construct absolute paths
    pkg_dir = pkgdir(EegFun)
    data_file = joinpath(pkg_dir, "resources", "data", "bdf", "example1.bdf")
    layout_file_path = joinpath(pkg_dir, "resources", "layouts", "biosemi", "biosemi72.csv")

    if isfile(data_file) && isfile(layout_file_path)
        println("Precompiling EegFun...")

        # 1. Reading raw data files (all formats) and opening databrowser
        # BDF
        dat = EegFun.read_raw_data(data_file)
        layout_file = EegFun.read_layout(layout_file_path)
        dat = EegFun.create_eegfun_data(dat, layout_file)

        # BrainVision
        bv_file = joinpath(pkg_dir, "resources", "data", "brainvision", "example1.vhdr")
        if isfile(bv_file)
            bv_raw = EegFun.read_raw_data(bv_file)
            bv_dat = EegFun.create_eegfun_data(bv_raw)
        end

        # EDF
        edf_file = joinpath(pkg_dir, "resources", "data", "edf", "test.edf")
        if isfile(edf_file)
            edf_raw = EegFun.read_raw_data(edf_file)
            edf_dat = EegFun.create_eegfun_data(edf_raw)
        end

        # FIF (raw and epochs)
        fif_raw_file = joinpath(pkg_dir, "resources", "data", "fif", "test_raw.fif")
        if isfile(fif_raw_file)
            fif_raw = EegFun.read_raw_data(fif_raw_file)
            fif_dat = EegFun.create_eegfun_data(fif_raw)
        end
        fif_epo_file = joinpath(pkg_dir, "resources", "data", "fif", "test_epochs.fif")
        if isfile(fif_epo_file)
            fif_epo = EegFun.read_raw_data(fif_epo_file)
            fif_epo_dat = EegFun.create_eegfun_data(fif_epo)
        end

        # XDF
        xdf_file = joinpath(pkg_dir, "resources", "data", "xdf", "test1.xdf")
        if isfile(xdf_file)
            xdf_raw = EegFun.read_raw_data(xdf_file)
            xdf_dat = EegFun.create_eegfun_data(xdf_raw)
        end

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
        res_db = EegFun.plot_databrowser(dat; display_plot = false)
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, NamedTuple())
        state_con = EegFun._create_browser_state(dat, dat.layout.data.label, res_db.ax, nothing, plot_kwargs)
        EegFun._xforward!(res_db.ax, state_con)
        EegFun._xback!(res_db.ax, state_con)
        EegFun._apply_hp_filter!(state_con)
        EegFun._apply_lp_filter!(state_con)
        EegFun._apply_filters!(state_con)
        EegFun._butterfly_plot!(res_db.ax, state_con)
        EegFun._subset_selected_data(state_con)
        EegFun._clear_all_selected_regions!(res_db.ax, state_con)
        EegFun._get_selected_regions_info(state_con)
        state_con.channels.visible[1] = false
        EegFun._draw(res_db.ax, state_con)

        res_db_epoch = EegFun.plot_databrowser(epochs[1]; display_plot = false)
        state_epoch = EegFun._create_browser_state(epochs[1], epochs[1].layout.data.label, res_db_epoch.ax, nothing, plot_kwargs)
        EegFun._step_epoch_forward(res_db_epoch.ax, state_epoch)
        EegFun._step_epoch_backward(res_db_epoch.ax, state_epoch)
        EegFun._repair_selected_channels!(state_epoch, [epochs[1].layout.data.label[1]], :neighbor_interpolation, res_db_epoch.ax)
        EegFun._repair_selected_channels!(state_epoch, [epochs[1].layout.data.label[1]], :spherical_spline, res_db_epoch.ax)
        EegFun._undo_last_repair!(state_epoch, res_db_epoch.ax)

        res_epochs = EegFun.plot_epochs(epochs[1]; display_plot = false)
        EegFun._handle_shared_navigation!(res_epochs.axes, :up, 0.2)
        EegFun._handle_shared_navigation!(res_epochs.axes, :down, 0.2)
        EegFun._handle_shared_navigation!(res_epochs.axes, :left, 0.2)
        EegFun._handle_shared_navigation!(res_epochs.axes, :right, 0.2)

        selection_state_epochs = EegFun.SharedSelectionState(res_epochs.axes; selection_color = :blue, selection_alpha = 0.3)
        EegFun._start_shared_selection!(first(res_epochs.axes), selection_state_epochs, 0.1)
        EegFun._update_shared_selection!(first(res_epochs.axes), selection_state_epochs, 0.1, 0.5)
        EegFun._finish_shared_selection!(first(res_epochs.axes), selection_state_epochs, 0.5)
        EegFun._get_epochs_selection_bounds(selection_state_epochs, [epochs[1]])
        EegFun._clear_shared_selection!(selection_state_epochs)

        res_erp = EegFun.plot_erp(erp; display_plot = false)
        EegFun._handle_shared_navigation!(res_erp.axes, :up, 0.2)
        EegFun._handle_shared_navigation!(res_erp.axes, :left, 0.2)

        selection_state_erp = EegFun.SharedSelectionState(res_erp.axes; selection_color = :blue, selection_alpha = 0.3)
        EegFun._start_shared_selection!(first(res_erp.axes), selection_state_erp, 0.1)
        EegFun._update_shared_selection!(first(res_erp.axes), selection_state_erp, 0.1, 0.5)
        EegFun._finish_shared_selection!(first(res_erp.axes), selection_state_erp, 0.5)
        EegFun._get_erp_selection_bounds(selection_state_erp, erp)
        EegFun._filter_visible_conditions(erp, Ref(nothing))
        EegFun._average_conditions(erp)
        EegFun._clear_shared_selection!(selection_state_erp)

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

        # 9. Optimized Statistical & Channel Operations
        EegFun.channel_average!(dat, channel_selections = [EegFun.channels([:Fp1, :Fp2])], output_labels = [:Frontal_Avg])
        EegFun.channel_difference!(
            dat,
            channel_selection1 = EegFun.channels(:Fp1),
            channel_selection2 = EegFun.channels(:Fp2),
            channel_out = :vEOG,
        )
        erp[1].condition = 1
        EegFun.condition_average(erp, [[1]])

        erp_fake2 = copy(erp[1])
        erp_fake2.condition = 2
        EegFun.condition_difference([erp[1], erp_fake2], [(1, 2)])

        # 10. ERP Metrics
        EegFun.gfp(erp[1])
        EegFun.global_dissimilarity(erp[1])
        EegFun.lrp(erp[1], erp_fake2; channel_selection = EegFun.channels([:C3, :C4]))
        EegFun.jackknife_average([erp[1], erp_fake2])

        println("Precompilation complete!")
    else
        println("Skipping precompilation workload (example data files not found)")
        println("  Expected data file: $data_file")
        println("  Expected layout file: $layout_file_path")
    end

end
