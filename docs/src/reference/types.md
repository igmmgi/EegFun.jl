# Core Types

Data structures and type definitions in EegFun.jl.

## Abstract Types

```@docs
EegFun.EegFunData
EegFun.EegData
EegFun.SingleDataFrameEeg
EegFun.MultiDataFrameEeg
EegFun.StatsResult
```

## Data Container Types

```@docs
EegFun.AnalysisInfo
EegFun.ContinuousData
EegFun.ErpData
EegFun.EpochData
EegFun.TimeFreqData
EegFun.TimeFreqEpochData
EegFun.SpectrumData
```

## Layout Types

```@docs
EegFun.Layout
EegFun.Neighbours
```

## Epoch Configuration

```@docs
EegFun.EpochCondition
```

## ICA Types

```@docs
EegFun.IcaPrms
EegFun.InfoIca
```

## Statistics Types

```@docs
EegFun.AnalysisData
EegFun.StatisticalData
EegFun.Cluster
EegFun.ClusterInfo
EegFun.Clusters
EegFun.TestInfo
EegFun.StatMatrix
EegFun.Masks
EegFun.PermutationDistribution
EegFun.PermutationResult
EegFun.AnalyticResult
```

## Decoding Types

```@docs
EegFun.DecodingParameters
EegFun.DecodedData
```

## RSA Types

```@docs
EegFun.NoiseCeiling
EegFun.RsaData
```

## Pipeline Configuration Types

```@docs
EegFun.FilterSection
EegFun.FilterConfig
EegFun.EogConfig
EegFun.EegConfig
EegFun.IcaConfig
EegFun.PreprocessConfig
```

## See Also

- [Data structures explanation](../explanations/data-structures.md)
