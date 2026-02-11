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

## Processing Info Types

```@docs
EegFun.FilterInfo
EegFun.BaselineInfo
EegFun.EpochInfo
EegFun.TriggerInfo
EegFun.Interval
EegFun.PlotLayout
```

## Artifact & Rejection Types

```@docs
EegFun.Rejection
EegFun.ArtifactInfo
EegFun.ArtifactComponents
EegFun.ZScoreRejectionInfo
EegFun.EpochRejectionInfo
EegFun.EpochRepairInfo
EegFun.ChannelRepairInfo
EegFun.ContinuousRepairInfo
EegFun.TemporalCluster
```

## GUI & State Types

```@docs
EegFun.EpochRejectionState
EegFun.SharedSelectionState
EegFun.TopoSelectionState
EegFun.PlotHelpInfo
```

## ERP Measurement Types

```@docs
EegFun.ErpMeasurementsResult
EegFun.ValidationResult
```

## Decoding & Statistics Result Types

```@docs
EegFun.DecodingStatisticsResult
```

## Batch & Pipeline Types

```@docs
EegFun.BatchConfig
EegFun.BatchResult
EegFun.AnalysisSettings
EegFun.ConfigParameter
EegFun.PipelineTemplateOptions
```

## Constants

```@docs
EegFun.VALID_MEASUREMENT_TYPES
EegFun.MEASUREMENT_TYPE_LABELS
EegFun._LATENCY_MEASUREMENT_TYPES
```

## See Also

- [Data structures explanation](../explanations/data-structures.md)
