## Overview

This demo demonstrates different re-referencing schemes for EEG data.

### Why Re-reference?

The reference electrode determines the voltage baseline:
- **Recording reference**: Often active during recording (Cz, mastoids, etc.)
- **Analysis reference**: Should match research question and conventions

### Common Reference Schemes

**Average Reference:**
- Subtract mean of all electrodes
- Assumes zero average scalp potential
- Most common for high-density EEG

**Mastoid Reference:**
- Reference to average of left/right mastoids
- Common in ERP research
- Reduces activity from reference sites

**Linked Mastoids:**
- Average of earlobes/mastoids
- Traditional ERP standard

**Cz Reference:**
- Single central electrode
- Useful for motor/sensory studies

### Best Practices

- Apply same reference across conditions/participants
- Re-reference before interpolating bad channels
- Consider effect on topographies
- Document reference choice in publications
