using EegFun

data_file = "my_raw_file.bdf"
layout_file = EegFun.read_layout("my_layout.csv");
EegFun.polar_to_cartesian_xy!(layout_file)
dat = EegFun.read_raw_data(data_file);
dat = EegFun.create_eeg_dataframe(dat, layout_file);

EegFun.plot_databrowser(dat);