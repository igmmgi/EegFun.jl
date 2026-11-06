This demo showcases the interactive widget controls, time scrubbers, and keyboard-driven control panels built into EegFun's core plotting viewers when using `GLMakie`.

### Included Interactive Modes

| Viewer / Mode | Interactive Controls |
| --- | --- |
| `plot_topography(...; scrubber = true)` | 2D scalp map with time slider, `<` / `Play` / `>` playback, and `Mean` button |
| `plot_topography_3d(...; scrubber = true)` | Rotatable 3D scalp mesh with synchronized real-time playback scrubber |
| `plot_channel_spectrum` | Live `X-Axis` (`Linear` / `Log10`) and `Y-Axis` (`Linear` / `Log10` / `dB`) dropdown menus |
| `plot_trigger_timing` | Dual sliders for scrolling through recording time and resizing the visible window |
| `plot_erp` / `plot_epochs` | Press `c` for live Baseline & Condition Control Panel, `l` to toggle legends, `s` for Save Figure dialog, `i` for help, or drag-select + right-click to launch topographic popups |

