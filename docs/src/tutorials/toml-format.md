# TOML Format

EegFun.jl uses [TOML](https://toml.io) (Tom's Obvious Minimal Language) for
its pipeline configuration files and epoch condition files. TOML is designed
to be easy to read and write — if you can read an INI file, you already know
most of what you need.

## Key–Value Pairs

The most basic building block. Values can be strings, numbers, or booleans:

```toml
directory = "./preprocessed_files"
epoch_start = -0.2
apply = true
```

Strings use double quotes. Numbers and booleans (`true` / `false`) are
unquoted.

## Comments

Lines starting with `#` are comments and are ignored:

```toml
# High-pass filter settings
freq = 0.1   # cutoff in Hz
```

## Tables (Sections)

Square brackets define a **table** (a group of related settings). Everything
below a table header belongs to that table, until the next header:

```toml
[files.input]
directory = "."
layout_file = "biosemi72.csv"

[files.output]
directory = "./preprocessed_files"
save_ica_data = true
```

Dots create **nested** tables. The example above is equivalent to:

```toml
[files]

  [files.input]
  directory = "."

  [files.output]
  directory = "./preprocessed_files"
```

## Arrays

Square brackets around values create arrays (lists):

```toml
raw_data_files = ["participant_01.bdf", "participant_02.bdf"]
```

Nested arrays are arrays of arrays:

```toml
# Each inner array lists alternative channel names (first match wins)
vEOG_channels = [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]]
```

## Arrays of Tables

Double square brackets `[[...]]` define an **array of tables** — a list of
identically structured sections. This is how epoch conditions are specified:

```toml
[epochs]

[[epochs.conditions]]
name = "standard"
trigger_sequences = [[1]]
reference_index = 1

[[epochs.conditions]]
name = "deviant"
trigger_sequences = [[2]]
reference_index = 1

[[epochs.conditions]]
name = "novel"
trigger_sequences = [[3, 4]]
reference_index = 1
timing_pairs = [[1, 2]]
min_interval = 0.1
max_interval = 2.0
```

Each `[[epochs.conditions]]` block adds one entry to the `conditions` list
inside the `epochs` table.

## Quick Reference

| Feature | Syntax | Example |
| --- | --- | --- |
| String | `"..."` | `"avg"` |
| Integer | bare number | `1` |
| Float | decimal number | `0.1` |
| Boolean | `true` / `false` | `true` |
| Array | `[...]` | `[1, 2, 3]` |
| Table | `[name]` | `[preprocess]` |
| Nested table | `[a.b]` | `[preprocess.eeg]` |
| Array of tables | `[[name]]` | `[[epochs.conditions]]` |
| Comment | `# ...` | `# cutoff in Hz` |

## Further Reading

- [TOML homepage](https://toml.io/en/) — introduction and full language reference
- [Pipeline configuration](../tutorials/workflows/pipeline_templates.md) — using TOML with EegFun.jl's preprocessing pipeline
