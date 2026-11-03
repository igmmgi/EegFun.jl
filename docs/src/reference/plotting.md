```@meta
CollapsedDocStrings = true
```

# Plotting & Interactive Demos

All public visualization functions, themes, topography helpers, and interactive teaching demos in `EegFun.jl`.

## Index

```@index
Pages = ["plotting.md"]
```

## Plotting Functions

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function]
Filter = t -> let s = applicable(nameof, t) ? string(nameof(t)) : string(t); startswith(s, "plot_") end
```

## Themes & Topography Overlays

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function]
Filter = t -> let s = applicable(nameof, t) ? string(nameof(t)) : string(t); s in ("theme_eegfun", "add_topo_rois!") end
```

## Interactive Teaching Demos

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function]
Filter = t -> let s = applicable(nameof, t) ? string(nameof(t)) : string(t); startswith(s, "signal_example_") || s == "simulate_erp" end
```
