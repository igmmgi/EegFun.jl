```@meta
CollapsedDocStrings = true
```

# Core Types

Data structures and type definitions in `EegFun.jl`.

## Abstract & Utility Types

```@docs
EegFun.EegFunData
EegFun.EegData
EegFun.SingleDataFrameEeg
EegFun.MultiDataFrameEeg
EegFun.EegFunError
EegFun.Interval
```

## Data Container & Metadata Types

```@docs
EegFun.AnalysisInfo
EegFun.BaselineInfo
EegFun.FilterInfo
EegFun.FilterSection
EegFun.TriggerInfo
EegFun.ContinuousData
EegFun.ErpData
EegFun.EpochData
EegFun.TimeFreqData
EegFun.TimeFreqEpochData
EegFun.SpectrumData
```

## Layout & Plotting Types

```@docs
EegFun.Layout
EegFun.Neighbors
EegFun.PlotLayout
```

## Epoch & RIDE Types

```@docs
EegFun.EpochCondition
EegFun.RideComponent
EegFun.RideResult
```

## ICA Types

```@docs
EegFun.IcaPrms
EegFun.InfoIca
EegFun.ArtifactComponents
```

## Statistics Types

```@docs
EegFun.StatisticalData
EegFun.StatsResult
EegFun.PermutationResult
EegFun.AnalyticResult
EegFun.TFStatisticalData
EegFun.TFStatsResult
EegFun.TFClusterPermutationResult
EegFun.TFAnalyticResult
```

## Decoding & RSA Types

```@docs
EegFun.DecodingParameters
EegFun.DecodedData
EegFun.DecodingStatisticsResult
EegFun.NoiseCeiling
EegFun.RsaData
```

## ERP Measurement Types

```@docs
EegFun.ErpMeasurementsResult
EegFun.ValidationResult
```

## Artifact, Rejection & Repair Types

```@docs
EegFun.ArtifactInfo
EegFun.Rejection
EegFun.EpochRejectionInfo
EegFun.EpochRejectionState
EegFun.ZScoreRejectionInfo
EegFun.ContinuousRepairInfo
EegFun.EpochRepairInfo
EegFun.ChannelRepairInfo
EegFun.AutoRejectConfig
EegFun.AutoRejectResult
```

## Batch, Configuration & Pipeline Types

```@docs
EegFun.AnalysisSettings
EegFun.BatchResult
EegFun.CleanLineConfig
EegFun.EegConfig
EegFun.EogConfig
EegFun.FilesConfig
EegFun.FilterConfig
EegFun.IcaConfig
EegFun.InputFilesConfig
EegFun.OutputFilesConfig
EegFun.PipelineConfig
EegFun.PipelineTemplateOptions
EegFun.PreprocessConfig
EegFun.ResampleConfig
```

## Constants

```@docs
EegFun.VALID_MEASUREMENT_TYPES
EegFun.MEASUREMENT_TYPE_LABELS
```

## See Also

- [Data structures](../explanations/data-structures.md)
