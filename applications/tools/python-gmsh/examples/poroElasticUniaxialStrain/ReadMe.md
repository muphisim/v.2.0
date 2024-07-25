Test case for poroelastic model. There are two primary modes of loading:

1) Step change in normal stress at x=0
2) Step change in pore pressure at x=0

Superposing these two modes is equivalent to loading by a fluid (Mode3). 

These loading comditions can be configured using the model.py file to create the appropriate input.inp file.

Each simualation is then started from inside their corresponding folder using ./run.sh from terminal. The ./run.sh file also cleans the outputs before the run. 


