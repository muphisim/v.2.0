import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import os
import sys

plt.close("all")
mpl.rcParams.update({"font.size":30})

plt.figure(figsize=(20,15))
data = pd.read_csv(os.path.join(sys.argv[1],"logIter.txt"))

plt.plot(data["Iteration"].values,data["maxDisp"].values,"-", label="Maximum displacement", color='k', linewidth=3)
plt.xlabel("Iteration")
plt.ylabel("Maximum Displacement")
plt.savefig("maxDisplacement.png")
plt.close()
