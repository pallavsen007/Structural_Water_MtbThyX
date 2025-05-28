### Import the necessary libraries
import MDAnalysis as mda
from MDAnalysis.analysis.waterdynamics import MeanSquareDisplacement as MSD
from MDAnalysis.analysis.waterdynamics import SurvivalProbability as SP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress

## Load the universe with trajectory file and the structure file
u = mda.Universe("./first_sp.gro", "./traj_sp.xtc")
u_msd = mda.Universe("./first_nojump.gro", "./msd_nojump.xtc")

## Define the selection for the water molecule around the substrate (dUMP)
select = "resname SOL and around 3.4 (resid 250 and name O4)"   ### for the SP analysis, we select the CW1/CW2 by defining an area (O4 for CW1 and O2P for CW2)
selection_msd = "resid 373"  ### for the MSD analysis, we select the specific water molecule

## Validate selection
selected_atoms = u.select_atoms(select)
print(f"Number of selected atoms: {len(selected_atoms)}")
if len(selected_atoms) == 0:
    raise ValueError("Selection returned no atoms. Check your selection string.")

## Validate trajectory range
n_frames = len(u.trajectory)
print(f"Number of frames in trajectory: {n_frames}")

## Perform MSD analysis
MSD_analysis = MSD(u_msd, selection_msd, 0, n_frames, 21)
MSD_analysis.run()

## Perform Survival Probability analysis
sp = SP(u, select, verbose=True)
tau_timeseries = sp.tau_timeseries
sp_timeseries = sp.sp_timeseries

## Calculate MSD/SP ratio
meansd = np.array(MSD_analysis.timeseries)
survival = np.array(sp_timeseries)
MSD_vs_SP = meansd/survival

### If one wants to calculate the diffusion coefficient for individual trajectories, then do comment out the suceeding lines
#y = MSD_vs_SP 
#x = np.arange(len(y))
# Perform linear regression
#slope, intercept, r_value, p_value, std_err = linregress(x, y)
#plt.plot(x,y)
#plt.show()
### Find out the diffusion coefficient
#D = slope/2
#print(D/100)  ### depending on the frames for your conversion


# Convert to DataFrame with column name "Value"
df = pd.DataFrame(MSD_vs_SP, columns=['value'])
print(df)

# Write to text files with control over formatting
df.to_csv('MSD_vs_SP_values_F.txt', sep='\t', index=True, float_format='%.6f')

### One can store the results individually in a text file as shown above
### for individual runs and then plot the average and with error associated with it

