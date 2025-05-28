## Import the necessary libraries
import MDAnalysis as mda
from MDAnalysis.analysis.waterdynamics import SurvivalProbability as SP

## Load the universe with trajectory file and the structure file
u = mda.Universe("./first_sp.gro", "./traj_sp.xtc")


## Define the selection for the water molecule around the substrate (dUMP)
select = "resname SOL and around 3.4 (resid 250 and name O4)" ### O4 for CW1 and O2P for CW2 

## Validate selection
selected_atoms = u.select_atoms(select)
print(f"Number of selected atoms: {len(selected_atoms)}")
if len(selected_atoms) == 0:
    raise ValueError("Selection returned no atoms. Check your selection string.")

## Validate trajectory range
n_frames = len(u.trajectory)
print(f"Number of frames in trajectory: {n_frames}")

## Initialize the SurvivalProbability analysis
sp = SP(u, select, verbose=True)
sp.run()

# Extract the tau_timeseries and sp_timeseries
tau_timeseries = sp.tau_timeseries
sp_timeseries = sp.sp_timeseries

# Print the results in the console
for tau, sp in zip(tau_timeseries, sp_timeseries):
    print("{time} {sp}".format(time=tau, sp=sp))

### One can store the results individually in a text file
### for individual runs and then plot the average and with error associated with it
