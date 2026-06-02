import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import pandas as pd


# Paramètres physiques
gamma = 1.4 ;
R     = 287.0 ;
end_tuyere = 3 ;

# Fonction de A(x)

def A(x) :
    return 1.0 + 2.2 * (x -1.5)**2 ;

def target(x, A_star) :
    return A(x)/A_star ;

# Fonction à résoudre analytiquement pour avoir M

def area_ratio(M, gamma):
    return (1/M) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * M**2))**((gamma+1)/(2*(gamma-1))) ;

# Recherche de M suivant la variation de x dans la tuyere
x_values = np.linspace(0, end_tuyere, 500) ;  # 500 points entre 0 et 3
mach_Value = [] ;
A_star = 1 ;
for x in x_values:
    if x < 1.5:  # convergent → subsonique
        M = brentq(lambda M: area_ratio(M, gamma) - target(x,A_star), 0.01, 1.0) ; # on résoud d'abord sur le domaine subsonique [0.01,1.0]
    else:  # divergent → supersonique
        M = brentq(lambda M: area_ratio(M, gamma) - target(x,A_star), 1.0, 5.0) ; # on résoud ensuite sur le domaine supersonique ]1, end_tuyere]
    mach_Value.append(M) ;



# Charge le dernier CSV de ta simulation
df = pd.read_csv("../cmake-build-debug/resultat_226.csv")
df['mach'] = df['vitesse'] / np.sqrt(1.4 * df['pression'] / df['rho'])

mach_interpolé = np.interp(df['x'], x_values, mach_Value)
erreur_mach = np.sqrt(np.mean((df['mach'] - mach_interpolé)**2))
print(f"Erreur L2 Mach : {erreur_mach:.4f}")

fig, ax1 = plt.subplots()

ax1.plot(x_values, mach_Value, label='analytique')
ax1.plot(df['x'], df['mach'], label=f"numérique (L2 = {erreur_mach:.4f})")
ax1.set_xlabel('Position x (m)')
ax1.set_ylabel('Mach number M')
ax1.legend()
ax1.grid()

plt.show()
