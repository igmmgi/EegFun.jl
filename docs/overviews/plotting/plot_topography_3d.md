This demo demonstrates interactive 3D scalp topography visualization (`plot_topography_3d`) for continuous, epoched, and averaged ERP data.

### 2D vs. 3D Topographic Maps

Unlike 2D `plot_topography` (which projects electrodes onto planar `x2, y2` coordinates), `plot_topography_3d` maps interpolated scalp potentials onto a 3D head mesh using spherical layout coordinates (`inc`, `azi`).

> **Note**: Use `GLMakie` when working interactively so you can rotate, pan, and zoom the 3D scene.

### Key Options

- **Data selection**: Select latencies or windows via `interval_selection = EegFun.times(t1, t2)` or predicates via `sample_selection`.
- **Camera orientation**: Set initial viewing angles in degrees with `camera_azimuth` and `camera_elevation`.
- **Styling**: Shares aesthetics with `plot_topography`, including `colormap`, `colorrange`, `colorbar_plot`, `point_plot`, `point_markersize`, `label_plot`, and `label_fontsize`.
- **Supported data types**: Works directly with `ContinuousData`, `EpochData` (specifying epoch index), and `ErpData`.

