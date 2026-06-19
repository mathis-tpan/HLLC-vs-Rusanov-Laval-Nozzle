import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib.animation import FuncAnimation
import numpy as np

# 1. Recherche des fichiers CSV générés par le C++

# Chemin vers results/ depuis scripts/
DOSSIER_CSV = os.path.join(os.path.dirname(__file__), '..', 'results')

fichiers = sorted(
    [f for f in os.listdir(DOSSIER_CSV)
     if f.startswith('resultat_') and f.endswith('.csv')],
    key=lambda x: int(x.split('_')[1].split('.')[0])
)
if not fichiers:
    print("Aucun fichier resultat_X.csv trouvé!")
    exit()

# 2. Configuration de la fenêtre graphique
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 10), sharex=True)


# 3. La fonction magique qui dessine chaque image (frame)
def update(frame_id):
    # On récupère le fichier correspondant à l'image actuelle
    f = fichiers[frame_id]
    data = pd.read_csv(os.path.join(DOSSIER_CSV, f))
    data['mach'] = data['vitesse'] / np.sqrt(1.4 * data['pression'] / data['rho'])
    # On nettoie les axes pour redessiner proprement
    ax1.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()

    # Graphique 1 : Densité
    ax1.plot(data['x'], data['rho'], color='royalblue', lw=2, label='Densité (rho)')
    ax1.set_ylabel('Densité')
    ax1.set_ylim(-0.05, 1.5)
    ax1.legend(loc='upper right')
    ax1.set_title(f"Tuyère de De Laval - Image {frame_id}/{len(fichiers) - 1}", fontsize=12, fontweight='bold')
    ax1.grid(True)

    # Graphique 2 : Vitesse
    ax2.plot(data['x'], data['vitesse'], color='crimson', lw=2, label='Vitesse (u)')
    ax2.set_ylabel('Vitesse')
    ax2.set_ylim(-1, 700)
    ax2.legend(loc='upper right')
    ax2.grid(True)

    # Graphique 3 : Pression
    ax3.plot(data['x'], data['pression'], color='forestgreen', lw=2, label='Pression (P)')
    ax3.set_ylabel('Pression')
    ax3.set_xlabel('Position dans le tube (x)')
    ax3.set_ylim(0, 150000)
    ax3.legend(loc='upper right')
    ax3.grid(True)

    # Graphique 4 : Mach
    ax4.plot(data['x'], data['mach'], color='purple', lw=2, label='Mach')
    ax4.set_ylabel('Mach')
    ax4.set_xlabel('Position dans le tube (x)')
    ax4.set_ylim(0, 5)
    ax4.legend(loc='upper right')
    ax4.grid(True)

# 4. Création de l'animation globale
print("Préparation du film...")
ani = FuncAnimation(fig, update, frames=len(fichiers), repeat=False)

# 5. ENREGISTREMENT DU FICHIER VIDEO (GIF)
print("Enregistrement de la vidéo en cours (cela peut prendre quelques secondes)...")
ani.save(str(DOSSIER_CSV / 'Laval.gif'), writer='pillow', fps=240)
print("🏁 Vidéo enregistrée avec succès sous le nom 'Laval.gif' !")

# 6. Affichage à l'écran
plt.show()