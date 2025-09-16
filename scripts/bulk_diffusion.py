### Import all the necessary libraries

import MDAnalysis as mda
from MDAnalysis.analysis.waterdynamics import MeanSquareDisplacement as MSD
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import linregress


# Load the universe i.e the trajectory file and the 
u = mda.Universe("./first.gro", "./nojump.xtc")

### Selects all the Oxygen atoms in the water box
select = "name OW"

# Validate selection
selected_atoms = u.select_atoms(select)
print(f"Number of selected atoms: {len(selected_atoms)}")
if len(selected_atoms) == 0:
    raise ValueError("Selection returned no atoms. Check your selection string.")

# Validate trajectory range
n_frames = len(u.trajectory)
print(f"Number of frames in trajectory: {n_frames}")

## Run the analyis here
MSD_analysis = MSD(u, select,0,n_frames,11)
MSD_analysis.run()
msd = MSD_analysis.timeseries

y = msd
x = np.arange(len(y))

# Perform linear regression
slope, intercept, r_value, p_value, std_err = linregress(x, y)

### Plot and check  
plt.plot(x,y)
plt.show()

### Print the self diffusion coefficient (D) of the water
print(y)
D = slope/6
print(D) 

# Convert to DataFrame with column name "Value"
df = pd.DataFrame(msd, columns=['value'])
print(df)

# Write to text files with control over formatting
df.to_csv('bulk_values.txt', sep='\t', index=True, float_format='%.6f')
