# eegfun.jl

Documentation for eegfun.jl

## Overview

eegfun is a Julia package for EEG data analysis and processing.

## Installation

```julia
using Pkg
Pkg.add("eegfun")
```

## Quick Start

```julia
using eegfun

# Load EEG data
data = read_biosemi_data("your_data.bdf")

# Create layout
layout = read_layout("biosemi64.csv")

# Preprocess data
preprocess_eeg_data("config.toml")
```

## API Reference

### Core Types

```@autodocs
Modules = [eegfun]
Filter = t -> typeof(t) <: Type
```

### Functions

```@autodocs
Modules = [eegfun]
Filter = f -> typeof(f) <: Function
```

## Index

```@index
```
