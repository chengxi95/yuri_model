## How to generate simulation data

`yuri_model_data_generate_fixed_initial.m` is used to generate simulation data, simply run the first code cell it will generate simulation data and saved it into a .mat file and there are several options in the script.

1. *total_trial_num* controls how many trials you want to simulate
2. *stimulation_type* controals which type of stimulation will be presented (1:blue triangle, 2:red circle, 3:green circle, 4:yellow triangle)
3. *VA_off, MD_off, pPFC_off, PV_off*, legion test, whether to turn off certain parts of the brain, the simulation will be for healthy subject if all set to 0.
4. *initialization*, to make sure the cells are the same for different types of simulation, we use the same external driving current for different types of simulation, that external driving current is saved in `initialization.mat`. **If set *initialization* is set to 0 (default), you need to load `initialization.mat` (simply double click the file) before running the cell.** If you want to have different initializations, set *initialization* to 1 and uncomment *clear all* in line 3.


## Testing and plotting script
- `yuri_model_data_generate_fixed_initial.m` also contains script to generate rater plot for certain neuron across all the trials or certain trial across all the neurons. You can change *neuron_id* in cell 3 and *trial_num* in cell 4 and run the corresponding cell to get the raster plot.

- `tuned_cell_spike_density_plot.m` is used to plot the spike density function for all the tuned cells. Simply load the saved simulation data (by double click the saved .mat file) then run `tuned_cell_spike_density_plot.m`

- `print_success_trial_ratio.m` is used to get the ratio of a successful trial in the saved data, a trial is considered successful if number of peak in testing is greater than the number of peak in
the baseline by at least 4 (default). And only if 20 cells (default) fires within 1ms (default) will be considered a peak. Currently this
script is only used on PFC ensemble. The rule of thumb is that if the stimulation type is blue triangle or red circle, the success_ratio for full_PFC_shape_ensemble is large while for full_PFC_ori_ensemble is low and vice versa. Simply laod the saved simulation data (by double click the saved .mat file) then run `print_success_trial_ratio`

- `calculate_si_all_cells.m` is used to get the selectivity index (SI) of all the cells in the saved data. The SI is defined in the paper. By default, it assumes the stimulatation is presented at 200000 bin (2s), the baseline is from 1s-2s and the testing is from 2s-3s. The window size is 0.1s and it only uses two rules (shape or orientation, n = 2). Choose the set of data you want to use by commenting out other bc_files, gc_files, rc_files, yt_files and run the script. 
