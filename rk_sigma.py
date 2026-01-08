import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

# Paramètres du problème
S = 200  # Nombre d'équations couplées
T = 1000  # Nombre d'itérations
h = 0.1  # Pas de temps
K = 10 #nombre de CI



# Fonction représentant le système d'équations différentielles couplées
def system(t, N):
    dNdt = np.zeros(S)
    for i in range(S):
        cof = 0
        for j in range(S):
            if j == i:
                cof = cof + 0
            else:
                cof = cof + alpha[i,j]*N[j]
        dNdt[i] = N[i]*(1-N[i]) - N[i]*cof 

    return dNdt

# Méthode de Runge-Kutta d'ordre 4 pour le système d'équations différentielles couplées
def runge_kutta(t, N, h):
    """
    Applique une étape de la méthode de Runge-Kutta d'ordre 4.

    t : temps actuel
    y : vecteur des variables y_i à l'instant t
    h : pas de temps

    Retourne le vecteur y à l'instant t + h.
    """
    k1 = h * system(t, N)
    k2 = h * system(t + 0.5 * h, N + 0.5 * k1)
    k3 = h * system(t + 0.5 * h, N + 0.5 * k2)
    k4 = h * system(t + h, N + k3)

    return N + (k1 + 2*k2 + 2*k3 + k4) / 6

def une_CI():
    # Conditions initiales
    N_init = np.random.rand(S)  # Initialiser les N_i avec des valeurs aléatoires entre 0 et 1
    t_init = 0  # Temps initial

    # Stockage des résultats
    N = N_init.copy()
    t = t_init

    # Intégration numérique sur T itérations
    for iteration in range(T):
        N = runge_kutta(t, N, h)
        t += h
    return N

#lis = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1, 1.6, 2, 4, 6, 8, 10]
lis = [0.5]
var = []
#s = [0.5, 0.9]
for s in tqdm(lis):
    alpha = np.random.normal(loc=4/S, scale= s/np.sqrt(S), size=(S, S))
    for i in range(S):
        for j in range(i):
            alpha[i,j] = alpha[j,i]
    Eq = []
    for i in tqdm(range(K)):
        Eq.append(une_CI())


    Variance = []
    for i in range(S):
        l = [Eq[j][i] for j in range(K)]
        Variance.append(np.std(l))

    #print("Intégration terminée.")
    #print(s)
    var.append(np.mean(Variance))
'''
lis1 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7]
var1 = [6.380919229242113e-08, 1.5114978426244188e-06, 0.00010615007394414178, 0.0010651338308433644, 0.0031826107934183057, 0.005115017304584216, 0.007893472704669866, 0.023676408446747316, 0.020060778200496406]

lis2 = [0.75, 0.8, 0.9, 1, 1.1, 1.2]
var2 = [0.039393853542964756, 0.13856986697984733, 0.25016234118544634, 0.3145663848593625, 0.560358917076399, 1.1512382735498894]


a, b = np.polyfit(lis1, var1, 1)
c, d, e, f = np.polyfit(lis2, var2, 3)
l1 = np.linspace(0, 0.8, 10000)
l2 = np.linspace(0.7, 1.22, 10000)
k1 = [a*i + b for i in l1]
k2 = [c*i*i*i + d*i*i + e*i +f for i in l2]
x = 0.7439
y = 0.017
''' 
#La partie du dessus sert à déterminer le comportement limite

plt.tick_params(axis='both', which='major', labelsize=16)
plt.plot(lis, var, marker ='x', color = 'black', linestyle='None', markersize=8, label = r'$\sqrt{\mathbb{V}~(\sigma)}$')
#plt.plot(l1, k1, '--', label = r'fit polynomial $\sigma < \sigma_c$')
#plt.plot(l2, k2, '--', label = r'fit polynomial $\sigma > \sigma_c$')
plt.plot(x,y, marker = '.', color = 'blue', markersize = 20, label = r'$\sigma_c = 0.74$')
plt.text(x -0.05, 0.09, r'$\left(\sqrt{\mathbb{V}~(\sigma_c)},~\sigma_c\right)$', ha='center', fontsize=18, color='blue')
plt.xlabel(r'$\sigma$', fontsize = 20)
plt.ylabel(r'$\sqrt{\mathbb{V}~(\sigma)}$', fontsize = 20)
plt.legend(fontsize = 20)
plt.show()
