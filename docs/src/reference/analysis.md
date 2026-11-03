```@meta
CollapsedDocStrings = true
```

# Analysis Functions

All public data access, I/O, preprocessing, ERP, RIDE, time-frequency, statistical, decoding, RSA, and pipeline functions in `EegFun.jl`.

## Index

```@index
Pages = ["analysis.md"]
```

## Functions

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Filter = t -> let s = applicable(nameof, t) ? string(nameof(t)) : string(t); !startswith(s, "plot_") && !startswith(s, "signal_example_") && s != "simulate_erp" && s != "theme_eegfun" && s != "add_topo_rois!" end
```
