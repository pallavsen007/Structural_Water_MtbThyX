### Import the necessary libraries
import MDAnalysis as mda
from MDAnalysis.analysis.waterdynamics import MeanSquareDisplacement as MSD
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import linregress


## Load the universe with trajectory file and the structure file
u = mda.Universe("./first.gro", "./nojump.xtc")

## Define the selection for any water molecule from the water box 
select = "resid 1000"

## Validate trajectory range
n_frames = len(u.trajectory)
print(f"Number of frames in trajectory: {n_frames}")



## Perform MSD analysis
MSD_analysis = MSD(u, select,0,n_frames,11)
MSD_analysis.run()
msd = MSD_analysis.timeseries

y = msd
x = np.arange(len(y))

# Perform linear regression
slope, intercept, r_value, p_value, std_err = linregress(x, y)
plt.plot(x,y)
plt.show()

D = slope/2
print(D/100)

# Convert to DataFrame with column name "Value"
df = pd.DataFrame(msd, columns=['value'])
print(df)

# Write to text files with control over formatting
df.to_csv('bulk_values.txt', sep='\t', index=True, float_format='%.6f')
