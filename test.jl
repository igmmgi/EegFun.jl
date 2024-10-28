using BioSemiBDF

# test filter
function test_filter()
  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  filter_data!(dat, "hp", 1, 2)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  filter_data!(dat, "lp", 1, 2)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  filter_data!(epochs, "hp", 1, 2)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  filter_data!(epochs, "lp", 1, 2)

  # test filter
  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = filter_data(dat, "hp", 1, 2)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = filter_data(dat, "lp", 1, 2)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = filter_data(epochs, "hp", 1, 2)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = filter_data(epochs, "lp", 1, 2)
end


function test_rereference()

  # test re-reference
  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  rereference!(dat, dat.layout.label, 1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  rereference!(dat, 1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  rereference!(dat, dat.layout.label, :Fp1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  rereference!(dat, :Fp1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  rereference!(dat, dat.layout.label, "Fp1")

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  rereference!(dat, "Fp1")

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  rereference!(dat, dat.layout.label, [1, 2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  rereference!(dat, [1, 2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  rereference!(dat, dat.layout.label, [:Fp1, :Fp2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  rereference!(dat, [:Fp1, :Fp2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  rereference!(dat, dat.layout.label, ["Fp1", "Fp2"])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  rereference!(dat, ["Fp1", "Fp2"])

  # test re-reference
  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  rereference!(epochs, epochs.layout.label, 1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  rereference!(epochs, 1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  rereference!(epochs, epochs.layout.label, :Fp1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  rereference!(epochs, :Fp1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  rereference!(epochs, epochs.layout.label, "Fp1")

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  rereference!(epochs, "Fp1")

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  rereference!(epochs, epochs.layout.label, [1, 2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  rereference!(epochs, [1, 2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  rereference!(epochs, epochs.layout.label, [:Fp1, :Fp2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  rereference!(epochs, [:Fp1, :Fp2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  rereference!(epochs, epochs.layout.label, ["Fp1", "Fp2"])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  rereference!(epochs, ["Fp1", "Fp2"])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = rereference(dat, dat.layout.label, 1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = rereference(dat, 1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = rereference(dat, dat.layout.label, :Fp1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = rereference(dat, :Fp1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = rereference(dat, dat.layout.label, "Fp1")

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = rereference(dat, "Fp1")

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = rereference(dat, dat.layout.label, [1, 2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = rereference(dat, [1, 2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = rereference(dat, dat.layout.label, [:Fp1, :Fp2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = rereference(dat, [:Fp1, :Fp2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = rereference(dat, dat.layout.label, ["Fp1", "Fp2"])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = rereference(dat, ["Fp1", "Fp2"])

  # test re-reference
  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = rereference(epochs, epochs.layout.label, 1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = rereference(epochs, 1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = rereference(epochs, epochs.layout.label, :Fp1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = rereference(epochs, :Fp1)

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = rereference(epochs, epochs.layout.label, "Fp1")

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = rereference(epochs, "Fp1")

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = rereference(epochs, epochs.layout.label, [1, 2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = rereference(epochs, [1, 2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = rereference(epochs, epochs.layout.label, [:Fp1, :Fp2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = rereference(epochs, [:Fp1, :Fp2])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = rereference(epochs, epochs.layout.label, ["Fp1", "Fp2"])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = rereference(epochs, ["Fp1", "Fp2"])

end


function test_baseline()

  # test baseline
  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  baseline!(dat, dat.layout.label, [0.01 0.02])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  baseline!(dat, [0.01 0.02])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  baseline!(epochs, epochs.layout.label, [-0.5 -0.5])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  baseline!(epochs, [-0.5 -0.5])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  erp = average_epochs(epochs)
  baseline!(erp, erp.layout.label, [-0.5 -0.5])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  erp = average_epochs(epochs)
  baseline!(erp, [-0.5 -0.5])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = baseline(dat, dat.layout.label, [0.01 0.02])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  dat1 = baseline(dat, [0.01 0.02])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = baseline(epochs, epochs.layout.label, [-0.5 -0.5])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  epochs1 = baseline(epochs, [-0.5 -0.5])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  erp = average_epochs(epochs)
  erp1 = baseline(erp, erp.layout.label, [-0.5 -0.5])

  dat = read_bdf("../Flank_C_3.bdf")
  dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
  epochs = extract_epochs(dat, 1, -0.5, 2)
  erp = average_epochs(epochs)
  erp1 = baseline(erp, [-0.5 -0.5])

end

# test_filter()
# test_rereference()
# test_baseline()



