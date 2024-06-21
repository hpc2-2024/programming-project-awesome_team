import numpy as np

n = 5 

A = np.diag(4 * np.ones(n**2)) + np.diag(-1 * np.ones(n**2-1),1) + np.diag(-1 * np.ones(n**2-1),-1) + np.diag(-1 * np.ones(n**2-n),n) + np.diag(-1 * np.ones(n**2-n),-n)  #+ np.diag(4 * np.ones(n**2))
