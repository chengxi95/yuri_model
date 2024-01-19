## How to generate simulation data

`yuri_model_data_generate_fixed_initial.m` is used to generate simulation data, simply run the first code cell it will generate simulation data and saved it into a .mat file and there are several options in the script.

1. *total_trial_num* controls how many trials you want to simulate
2. *stimulation_type* controals which type of stimulation will be presented (1:blue triangle, 2:red circle, 3:green circle, 4:yellow triangle)
3. *VA_off, MD_off, pPFC_off, PV_off*, legion test, whether to turn off certain parts of the brain, the simulation will be for healthy subject if all set to 0.
4. *initialization*, to make sure the cells are the same for different types of simulation, we use the same external driving current for different types of simulation, that external driving current is saved in `initialization.mat`. **If set *initialization* is set to 0 (default), you need to load `initialization.mat` (simply double click the file) before running the cell.** If you want to have different initializations, set *initialization* to 1 and uncomment *clear all* in line 3.


## Testing and plotting script
`tuned_cell_spike_density_plot.m` is used to plot the spike density function for all the tuned cells. Simply load the saved simulation data (by double click the .mat file) then run `tuned_cell_spike_density_plot.m`
