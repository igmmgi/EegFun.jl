# Why Julia for EEG Research?

Most EEG analysis today happens in MATLAB (EEGLAB, FieldTrip, ERPLAB) or Python
(MNE-Python). These are well-established packages with many years of
development, large user communities, and broad feature sets — they will
typically offer more functionality than EegFun.jl. Each ecosystem has clear
strengths, and researchers should choose the tool that best fits their needs.
This page explains why EegFun.jl is built on Julia and where the trade-offs
lie.

## Strengths

**Fully open source — language and libraries.**
Julia itself and its entire package ecosystem are open source. MATLAB-based
toolboxes such as EEGLAB, FieldTrip, and ERPLAB are open source, but they
still require a proprietary MATLAB licence to run. GNU Octave provides a free
alternative, but compatibility with these toolboxes is not 100%, and
performance and graphics support can differ significantly.

**One language from prototype to production.**
Scientific computing often hits a "two-language problem": researchers prototype
in a high-level language and then rewrite performance-critical code in C or
Fortran. In practice this means R code calling Rcpp/C++, Python relying on C
extensions (NumPy, SciPy, MNE's compiled backends), and MATLAB using MEX
files. Julia is designed to avoid this split — user-level code compiles to
efficient native code, so there is no need to drop into a second language for
speed. In EegFun.jl, there is no hidden C or Fortran layer; the analysis code
you read is the code that runs.

**Performance without boilerplate.**
Julia's JIT compiler generates machine code specialised to the types you
actually use. For numerically intensive EEG workflows — filtering, time-
frequency decomposition, permutation statistics — this often matches C/Fortran
performance without manual memory management or type annotations.

**Multiple dispatch for clean extensibility.**
Julia's type system and multiple dispatch make it straightforward to extend
existing functions to new data types. In EegFun.jl this means, for example,
the same `lowpass_filter!` function works on continuous data, epoched data, and ERP
averages with no wrapper boilerplate.

**Built-in package manager and reproducibility.**
Every Julia project records its exact dependency versions in `Manifest.toml`.
Collaborators can reproduce your environment with a single `Pkg.instantiate()`
call, without conda environments or Docker containers.

## Trade-offs

**Time to first execution (TTFX).**
Julia compiles code the first time it runs in a session, which means the first
call to a function can take noticeably longer than subsequent calls. Python
users familiar with Numba will recognise this pattern — both compile code to
fast machine code on first use — but in Julia the compilation applies to the
entire language rather than individual decorated functions. This "compilation
latency" has improved substantially in recent Julia versions (1.9+ introduced
package images, and 1.12 brings further improvements), but it remains more
noticeable than starting a MATLAB or Python session. In practice, most users
start a Julia session once and keep it running.

**Smaller ecosystem.**
The Julia ecosystem for EEG is young compared to MATLAB and
Python. Packages like EEGLAB (2004), FieldTrip (2011), and MNE-Python (2014)
have over a decade of development, thousands of users, and extensive
functionality that EegFun.jl does not yet match. There are also fewer Stack
Overflow answers and fewer colleagues down the hall who can help debug your
code. This is improving steadily (and is part of the motivation for
EegFun.jl), but it is an honest limitation.

**Tooling is still maturing.**
IDE support, debugging, and profiling tools exist (VS Code with the Julia
extension is quite capable), but they are not as mature or polished as the
equivalent MATLAB IDE or established Python tool-chains (PyCharm, Jupyter
ecosystem). The debugger, in particular, can be slow on large codebases.

**Package loading time.**
Large packages with many dependencies take time to load (`using EegFun` is not
instant). The very first `using EegFun` after installation is particularly
slow because Julia precompiles (partly) the package and all of its dependencies — this
is a one-off cost. Subsequent `using EegFun` calls in new sessions are
considerably faster, but still slower than `import mne` in Python.

> [!TIP]
> If you are coming from MATLAB or Python, the
> [MATLAB–Python–Julia cheat sheet](https://cheatsheets.quantecon.org/) is a
> helpful reference for translating familiar idioms.
