#let config = (
  skip-default-titlepage: true,
  header-mode: "none",  // Optional: hide header if needed, but default is fine
)

// Custom Cover Page (rendered before the #show rule and Documenter setup)
#align(center + horizon)[
  #grid(
    columns: (1fr, 1fr),
    align: horizon,
    [
      #text(size: 4em, weight: "bold", fill: rgb("d42b60"))[EegFun.jl]\
      #v(0.5em)
      #text(size: 2.2em, weight: "bold", fill: luma(80))[EEG/ERP analysis in Julia]
    ],
    [
      #align(right)[
        // We use absolute path to the build directory's public assets because Typst 
        // resolves paths relative to the generated .typ file in the build directory.
        // DocumenterTypst copies `src/assets` to `build_pdf/assets`.
        #image("public/EegFunLogo.png", width: 85%)
      ]
    ]
  )
  
  #v(5em)
  
  #grid(
    columns: (1fr, 1fr, 1fr, 1fr),
    gutter: 2.5em,
    align: left + top,
    [
      #text(weight: "bold", size: 1.3em, fill: rgb("d42b60"))[Raw data readers]\
      #v(0.5em)
      #text(size: 1.1em, fill: luma(60))[Biosemi, BrainVision Analyser, EDF, BIDS, EEGLAB, Fieldtrip ...]
    ],
    [
      #text(weight: "bold", size: 1.3em, fill: rgb("d42b60"))[Data Processing]\
      #v(0.5em)
      #text(size: 1.1em, fill: luma(60))[Filtering, referencing, artifact detection/correction, ICA, epoching ...]
    ],
    [
      #text(weight: "bold", size: 1.3em, fill: rgb("d42b60"))[Analysis]\
      #v(0.5em)
      #text(size: 1.1em, fill: luma(60))[Epochs, time-frequency analysis, ERP measurements, statistical testing]
    ],
    [
      #text(weight: "bold", size: 1.3em, fill: rgb("d42b60"))[Visualization]\
      #v(0.5em)
      #text(size: 1.1em, fill: luma(60))[Interactive plots, topographic maps, and publication-quality figures with Makie.jl]
    ]
  )
]

#pagebreak()
