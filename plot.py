import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("output.csv")
plt.scatter(df["T"], df["prob"])
plt.show()