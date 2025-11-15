import numpy as np
import matplotlib.pyplot as plt

# ---- Inputs (ring-load case numbers you used) ----
a = 3.25            # in (effective clamp radius / short side)
E = 10_000_000.0    # psi
nu = 0.33
q = 150.0           # psi
FOS = 2
q_eff = q * FOS     # psi
K = 2.719618e-6     # Roark coefficient (ring-load case)
delta_allow = -0.001  # in

# ---- Model ----
def D(t):      # flexural rigidity
    return E * t**3 / (12.0 * (1.0 - nu**2))

def delta(t):  # deflection (negative = downward)
    return -(q_eff * a**4 / D(t)) * K

# Thickness where |δ| = 0.001 in (solve δ = -C/t^3)
Cconst = (q_eff * a**4 * K) * 12.0 * (1.0 - nu**2) / E
t_req  = (abs(Cconst) / abs(delta_allow))**(1/3)
y_req  = delta(t_req)

# ---- Plot window focused around the intersection ----
x_left  = min(0.02, 0.5 * t_req)   # a little left of t_req
x_right = max(0.25, 6.0 * t_req)   # show some right of t_req
t = np.linspace(x_left, x_right, 600)
y = delta(t)

plt.figure(figsize=(9,5.2))
plt.plot(t, y, lw=2.5, label="Deflection vs Thickness")
plt.axhline(delta_allow, ls="--", lw=2, label="Max Allowable Deflection (-0.001 in)")
plt.scatter([t_req], [y_req], s=70, marker="x", label=f'Intersection t = {t_req:.3f}"')

plt.title("Intersection view (inches): shows t where δ = -0.001 in")
plt.xlabel("Plate Thickness (inches)")
plt.ylabel("Deflection (inches)")
plt.grid(True, alpha=0.3)
plt.legend(loc="best")

# y-limits for a clear view
ymin = min(y.min(), delta_allow) - 0.0005
plt.ylim(ymin, 0.0002)

plt.tight_layout()
plt.savefig("curve_intersection_ringK.png", dpi=160)
plt.show()

print(f"Intersection thickness t_req ≈ {t_req:.4f} in  (δ = {y_req:.6f} in)")