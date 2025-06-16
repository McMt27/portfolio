# Exercice 1 : Bornage des racines réelles
def inter(P):
    """Retourne un entier b tel que toutes les racines réelles de P soient dans [-b, b]."""
    deg = P[0]
    coeffs = P[1:]
    a_n = coeffs[-1]
    if a_n == 0:
        raise ValueError("Le coefficient dominant ne peut pas être nul.")
    max_ratio = max(abs(coeffs[i] / a_n) for i in range(len(coeffs) - 1))
    return int(max(1, 1 + max_ratio))

# Exercice 2 : Dérivation d’un polynôme
def deriv(P):
    """Retourne le polynôme dérivé de P."""
    deg = P[0]
    coeffs = P[1:]
    if deg == 0:
        return [0, 0]
    der = [deg - 1] + [coeffs[i+1] * (i + 1) for i in range(deg)]
    return der

# Exercice 3 : Division euclidienne de polynômes
def euclide(P, Q):
    """Retourne le reste de la division euclidienne de P par Q."""
    deg_P, coeffs_P = P[0], P[1:]
    deg_Q, coeffs_Q = Q[0], Q[1:]

    # Convertir en polynôme normal (coeffs uniquement) pour manipulations
    R = coeffs_P[:]
    while len(R) >= len(coeffs_Q):
        d = len(R) - len(coeffs_Q)
        lead_coeff = R[-1] / coeffs_Q[-1]
        for i in range(len(coeffs_Q)):
            R[d + i] -= lead_coeff * coeffs_Q[i]
        while R and abs(R[-1]) < 1e-12:
            R.pop()
    return [len(R)-1] + R if R else [0, 0]

# Exercice 4 : Suite de Sturm
def poly_neg(P):
    """Retourne l’opposé du polynôme P."""
    return [P[0]] + [-c for c in P[1:]]

def sturm(P):
    """Retourne la suite de Sturm du polynôme P."""
    S = [P, deriv(P)]
    while S[-1][1:] != [0]:
        r = euclide(S[-2], S[-1])
        S.append(poly_neg(r))
    return S

# Exercice 5 : Nombre de racines réelles distinctes
def sign(x):
    """Retourne le signe de x."""
    if x > 0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0

def count_sign_changes(lst, x):
    """Compte les changements de signe dans la suite de polynômes à la valeur x."""
    values = []
    for P in lst:
        val = sum(P[i+1] * (x ** i) for i in range(P[0] + 1))
        values.append(sign(val))
    changes = 0
    for i in range(len(values) - 1):
        if values[i] != 0 and values[i] * values[i+1] < 0:
            changes += 1
    return changes

def NBracines(P):
    """Retourne le nombre de racines réelles distinctes de P."""
    suite = sturm(P)
    b = inter(P)
    V_a = count_sign_changes(suite, -b)
    V_b = count_sign_changes(suite, b)
    return V_a - V_b

# Exercice 6 : Méthode de dichotomie pour une racine
def eval_poly(P, x):
    """Évalue le polynôme P en x."""
    return sum(P[i+1] * x ** i for i in range(P[0] + 1))

def solu(P, epsilon=0.1):
    """Retourne une solution approchée à 0.1 près de P(x)=0 par dichotomie."""
    b = inter(P)
    a, c = -b, b
    if eval_poly(P, a) * eval_poly(P, c) > 0:
        return None  # Aucune racine réelle dans l'intervalle
    while (c - a) > epsilon:
        m = (a + c) / 2
        if eval_poly(P, a) * eval_poly(P, m) <= 0:
            c = m
        else:
            a = m
    return (a + c) / 2