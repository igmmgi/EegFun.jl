```@meta
CollapsedDocStrings = true
```

# Core Types

Data structures and type definitions in EegFun.jl.

## Abstract Types

```@docs
EegFun.EegFunData
EegFun.EegData
EegFun.SingleDataFrameEeg
EegFun.MultiDataFrameEeg
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
EegFun.StatsResult
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

## ERP Measurement Types

```@docs
EegFun.ErpMeasurementsResult
EegFun.ValidationResult
```

## Artifact & Rejection Types

```@docs
EegFun.ArtifactInfo
EegFun.EpochRejectionInfo
EegFun.ContinuousRepairInfo
EegFun.ZScoreRejectionInfo
EegFun.EpochRepairInfo
EegFun.ChannelRepairInfo
```

## Batch & Pipeline Types

```@docs
EegFun.AnalysisSettings
EegFun.PipelineTemplateOptions
```

## Decoding & Statistics Result Types

```@docs
EegFun.DecodingStatisticsResult
```

## Constants

```@docs
EegFun.VALID_MEASUREMENT_TYPES
EegFun.MEASUREMENT_TYPE_LABELS
```

## See Also

- [Data structures](../explanations/data-structures.md)
