using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# basic trigger count from raw file
count = EegFun.trigger_count(dat)

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout_file);

# basic trigger count from EegFun data structure
count = EegFun.trigger_count(dat)

# trigger overview
EegFun.plot_trigger_overview(dat)

# trigger overview with ignored triggers
EegFun.plot_trigger_overview(dat; ignore_triggers = [3, 253])

# trigger timing i.e, when did each trigger occur and interval between triggers
EegFun.plot_trigger_timing(dat)

# trigger timing with ignored triggers (timing interval is updated)
EegFun.plot_trigger_timing(dat; ignore_triggers = [3])
