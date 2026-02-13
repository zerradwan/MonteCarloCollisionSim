import numpy as np
from scipy import integrate
from numpy.linalg import inv, det, LinAlgError

def unit_vector(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("Cannot normalize zero vector")
    return v / n

def collision_plane_basis(rel_vel):
    n = unit_vector(rel_vel)
    if abs(n[0]) < 0.9:
        a = np.array([1.0, 0.0, 0.0])
    else:
        a = np.array([0.0, 1.0, 0.0])
    u = a - np.dot(a, n) * n
    u = u / np.linalg.norm(u)
    v = np.cross(n, u)
    v = v / np.linalg.norm(v)
    return np.column_stack((u, v))

def project_covariance_to_plane(C3, B):
    return B.T @ C3 @ B


def pc_monte(mu, C, R, N = 1000000):
    samples = np.random.multivariate_normal(mu, C, size=N)
    dist = samples[:, 0]**2 + samples[:, 1]**2
    inside = np.sum(dist <= R**2)
    pc = inside/N
    return pc

def compute_pc_monte(
    primary_pos, secondary_pos,
    primary_cov3, secondary_cov3,
    rel_vel,
    primary_radius_m, secondary_radius_m,
    N=200000
):
    primary_pos = np.asarray(primary_pos, float)
    secondary_pos = np.asarray(secondary_pos, float)
    rel_pos3 = secondary_pos - primary_pos

    B = collision_plane_basis(rel_vel)

    C1_2 = project_covariance_to_plane(primary_cov3, B)
    C2_2 = project_covariance_to_plane(secondary_cov3, B)
    C2 = C1_2 + C2_2

    mu2 = B.T @ rel_pos3

    R_km = (primary_radius_m + secondary_radius_m) / 1000.0

    pc = pc_monte(mu2, C2, R_km, N=N)

    return {
        "Pc": pc,
        "mu_plane_km": mu2,
        "C_plane_km2": C2,
        "R_km": R_km,
    }

'''
if __name__ == "__main__":

    primary_pos = np.array([7000.000, 0.000, 0.000])       
    secondary_pos = np.array([7000.002, 0.003, 0.000])     

    primary_cov3 = np.diag([2e-6, 2e-6, 2e-6])   
    secondary_cov3 = np.diag([3e-6, 3e-6, 3e-6]) 

    rel_vel = np.array([-0.001, 7.5, 0.0])       

    primary_radius_m = 5.0   
    secondary_radius_m = 1.0 

    out = compute_pc_monte(
        primary_pos,
        secondary_pos,
        primary_cov3,
        secondary_cov3,
        rel_vel,
        primary_radius_m,
        secondary_radius_m
    )

    print("==== RESULTS ====")
    print("Probability of Collision Pc:", out["Pc"])
    print("Collision-plane mean (km):", out["mu_plane_km"])
    print("Collision-plane covariance (km^2):\n", out["C_plane_km2"])
    print("Combined hard-body radius R (km):", out["R_km"])
'''