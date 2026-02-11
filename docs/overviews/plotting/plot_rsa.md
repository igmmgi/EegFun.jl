This demo shows how to visualise Representational Similarity Analysis (RSA) results.

### What is RSA?

- **RDM** (Representational Dissimilarity Matrix) captures pairwise differences between neural representations
- **RDM Timecourse** shows how representational structure evolves over time
- **Model Correlations** compare observed RDMs against theoretical models

### Key Functions

| Function | Purpose |
| --- | --- |
| `plot_rdm_heatmap(rsa_result)` | Display an RDM as a heatmap |
| `plot_rdm_timecourse(rsa_result)` | Dissimilarity over time |
| `plot_model_correlations(rsa_result)` | Model comparison over time |

### Key Parameters

| Parameter | Function | Description |
| --- | --- | --- |
| `time_point` | `plot_rdm_heatmap` | Specific time (seconds or index); `nothing` = average |
| `condition_pairs` | `plot_rdm_timecourse` | Which pairs to plot; `:all` = all |
| `show_colorbar` | All | Show/hide colour bar |
| `colormap` | All | Colour map for the heatmap |

### What You'll Learn

1. Plotting RDM heatmaps at specific time points or averaged
2. Visualising dissimilarity timecourses for selected condition pairs
3. Comparing neural RDMs against theoretical model predictions
