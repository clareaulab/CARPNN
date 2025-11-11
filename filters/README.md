## This directory contains BindCraft-like filters for Boltz output

### `reference_BindCraft_default_filters`
This the original filter file from BindCraft in which the Boltz filters are all based on. It's for reference purpose only and does not work as a filter

### `Boltz1_default_filters.json`
This is the Boltz filter criteria we adopted based on the BindCraft default filters. Note the i_pae value in BindCraft is normalized (i_pae = pae_interaction / 31) and we use the equivalent unnormalized pAE interaction as filter.

### `Boltz1_unrelaxed_filters.json`
If the predicted structures are not relaxed, the interface dG could be higher than usual. As such we relax this condition in this filter

### `Boltz1_default_filters_Lowered_UnsatHBonds.json`
Empirically we observe often times binders would fail this biophysical filter but have other wise good metrics so we slightly relaxed this filter

