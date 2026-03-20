```@meta
CollapsedDocStrings = true
```

# Analysis Functions

All public analysis and processing functions in EegFun.jl.

## Functions

```@autodocs
Modules = [EegFun]
Order = [:function, :macro]
Filter = t -> let s = string(t); !startswith(s, "_") && !startswith(s, "plot_") && !startswith(s, "Base.") && !endswith(s, ".show") && !endswith(s, ".copy") && !endswith(s, ".convert") && !endswith(s, ".getproperty") end
```
