###############################################################
# Filter functions 
function _apply_filter!(dat::DataFrame, columns, filter)
    for col in names(dat)
        if col in columns
            dat[:, col] .= filtfilt(filter, dat[:, col])
        end
    end
end

function filter_data!(dat::DataFrame, columns, filter_type, freq, order, sample_rate)
    if filter_type == "hp"
        filter = digitalfilter(Highpass(freq, fs = sample_rate), Butterworth(order))
    elseif filter_type == "lp"
        filter = digitalfilter(Lowpass(freq, fs = sample_rate), Butterworth(order))
    end
    _apply_filter!(dat, columns, filter)
end

function filter_data(dat::DataFrame, columns, type, freq, order, sample_rate)
    dat_out = deepcopy(dat)
    filter_data!(dat_out, columns, type, freq, order, sample_rate)
    return dat_out
end

function filter_data!(dat::Union{ContinuousData,ErpData}, type, freq, order)
    filter_data!(dat.data, dat.layout.label, type, freq, order, dat.sample_rate)
end

function filter_data(dat::Union{ContinuousData,ErpData}, type, freq, order)
    dat_out = deepcopy(dat)
    filter_data!(dat_out.data, dat_out.layout.label, type, freq, order, dat_out.sample_rate)
    return dat_out
end


function filter_data!(dat::EpochData, type, freq, order)
    for epoch in eachindex(dat.data)
        filter_data!(dat.data[epoch], dat.layout.label, type, freq, order, dat.sample_rate)
    end
end

function filter_data(dat::EpochData, type, freq, order)
    dat_out = deepcopy(dat)
    filter_data!(dat_out, type, freq, order)
    return dat_out
end


function filter_data!(dat::EpochData, columns, type, freq, order, sample_rate)
    for epoch in eachindex(dat.data)
        filter_data!(dat.data[epoch], columns, type, freq, order, sample_rate)
    end
end

function filter_data(dat::EpochData, columns, type, freq, order, sample_rate)
    dat_out = deepcopy(dat)
    for epoch in eachindex(dat.data)
        filter_data!(dat.data[epoch], columns, type, freq, order, sample_rate)
    end
    return dat_out
end


function filter_data(dat::EpochData, type, freq, order)
    dat_out = deepcopy(dat)
    filter_data!(dat_out, type, freq, order)
    return dat_out
end

function filter_data!(dat::Vector{DataFrame}, columns, type, freq, order, sample_rate)
    for epoch in eachindex(dat)
        filter_data!(dat[epoch], columns, type, freq, order, sample_rate)
    end
end

function filter_data(dat::Vector{DataFrame}, columns, type, freq, order, sample_rate)
    dat_out = deepcopy(dat)
    filter_data!(dat_out, columns, type, freq, order, sample_rate)
    return dat_out
end
