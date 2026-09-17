module EegFunDimensionalDataExt

using EegFun
using DimensionalData

function DimensionalData.DimArray(dat::EegFun.SingleDataFrameEeg)
    A = Matrix(dat)
    t = dat.data.time
    c = EegFun.channel_labels(dat)
    
    # Dimensions: Time × Channel
    return DimArray(A, (Ti(t), Dim{:Channel}(c)))
end

function DimensionalData.DimArray(dat::EegFun.MultiDataFrameEeg)
    A = Array(dat)
    ep = 1:EegFun.n_epochs(dat)
    t = dat.data[1].time
    c = EegFun.channel_labels(dat)
    
    # Dimensions: Epoch × Time × Channel
    return DimArray(A, (Dim{:Epoch}(ep), Ti(t), Dim{:Channel}(c)))
end

end
