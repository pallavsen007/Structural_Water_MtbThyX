## Import the necessary libraries
import MDAnalysis as mda
from MDAnalysis.analysis.waterdynamics import WaterOrientationalRelaxation as WOR
import numpy as np
import pandas as pd

## Load the universe with trajectory file and the structure file
u = mda.Universe("./npt.gro", "./traj.xtc")

## Define the selection for the water molecule around the substrate (dUMP)
select = "resname SOL and around 3.4 (resid 250 and name O4)"  ### O4 for CW1 and O2P for CW2 

## Validate selection
selected_atoms = u.select_atoms(select)
print(f"Number of selected atoms: {len(selected_atoms)}")
if len(selected_atoms) == 0:
    raise ValueError("Selection returned no atoms. Check your selection string.")

## Validate trajectory range
n_frames = len(u.trajectory)
print(f"Number of frames in trajectory: {n_frames}")

## Initialize the Water Orientation Relaxation analysis
WOR_analysis = WOR(u, select, 0, n_frames, 21)
WOR_analysis.run()

## Collect the timeseries data
time_series = list(WOR_analysis.timeseries)
WOR_dip = np.array([entry[2] for entry in time_series])

## Convert to DataFrame with column name "Value"
df_1 = pd.DataFrame(WOR_dip, columns=['value'])
print(df_1)

## Write to text files with control over formatting
df_1.to_csv('dipole_value.txt', sep='\t', index=True, float_format='%.6f')

### One can store the results individually in a text file as shown above
### for individual runs and then plot the average and with error associated with it


