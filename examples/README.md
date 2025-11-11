## This directory contains example input/output frome ach step of the workflow

Please see the notebooks in the `notebooks` directory for detail on how the files are created. But breifly:

`00_bindcraft_outputs`: examples binder outputs from 2 BindCraft runs

`01_bindcraft_boltz_refolded`: binders from `00_bindcraft_outputs` but refolded and new lead binders are identified

`02_lead_candidate_mutagenized`: Mutagenesis of the lead binder identified from `01_bindcraft_boltz_refolded`

`03_lead_mutagenesis_boltz_refolded`: binders from `02_lead_candidate_mutagenized` but refolded and filters were applied

`04_final_selection`: The final list of binders identified from `03_lead_mutagenesis_boltz_refolded` to be experimentally tested.

Note that these are 10 random binders we grabbed from our past BindCraft campaigns and none have been experimentally tested.