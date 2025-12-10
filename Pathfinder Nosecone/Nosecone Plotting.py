import numpy as np
import matplotlib.pyplot as plt

def rotated_profile(filename_points='surface_points.csv',
                    filename_meridian='meridian_curve.csv',
                    y_min=1e-6, y_max=12-1e-6,
                    Ny=32, Ntheta=32):
    """
    Generate (x, y, z) coordinates for the surface formed by rotating
    x(y) = -3 * ((y - 12) / 12)^0.6 about the y-axis, for 0 < y < 12.

    Saves:
      - filename_points: CSV of all surface points (x,y,z)
      - filename_meridian: CSV of the generating curve (x(y), y, 0)
    """

    # Radius function r(y) = x(y).
    # For fractional power with negative base, use sign * |base|^0.6 to stay in reals.
    def r_of_y(y):
        base = (y - 12.0) / 12.0  # in (-1, 0) over (0, 12)
        return -3.0 * np.sign(base) * (np.abs(base) ** 0.6)

    # Sample parameters
    y = np.linspace(y_min, y_max, Ny)
    theta = np.linspace(0.0, 2.0 * np.pi, Ntheta, endpoint=False)

    # Meridian curve
    r = r_of_y(y)            # can be negative (that’s fine)
    x_meridian = r
    z_meridian = np.zeros_like(y)

    # Save meridian curve as (x, y, z=0)
    meridian = np.column_stack([x_meridian, y, z_meridian])
    np.savetxt(filename_meridian, meridian, delimiter=',')

    # Create surface of revolution
    Y, Theta = np.meshgrid(y, theta, indexing='ij')
    R = r_of_y(Y)

    X = R * np.cos(Theta)
    Z = R * np.sin(Theta)

    # Save all surface points
    surface_points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    np.savetxt(filename_points, surface_points, delimiter=',')

    # Basic summaries
    r_abs = np.abs(r)
    print(f'y range: [{y_min:.6f}, {y_max:.6f}] (Ny={Ny})')
    print(f'theta samples: {Ntheta}')
    print(f'min |radius|: {r_abs.min():.6f}, max |radius|: {r_abs.max():.6f}')
    print(f'points saved: {surface_points.shape[0]} -> {filename_points}')
    print(f'meridian saved: {meridian.shape[0]} -> {filename_meridian}')

    # Plot
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(121)
    ax.plot(x_meridian, y, 'k', label='x(y)')
    ax.plot(-x_meridian, y, 'k:', alpha=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Generating Curve')
    ax.axis('equal')
    ax.grid(True, alpha=0.3)

    ax3d = fig.add_subplot(122, projection='3d')
    # Thin wireframe for speed; switch to plot_surface for solid
    step_y = max(1, Ny // 80)
    step_t = max(1, Ntheta // 120)
    ax3d.plot_wireframe(X[::step_y, ::step_t], Y[::step_y, ::step_t], Z[::step_y, ::step_t],
                        rstride=1, cstride=1, color='tab:blue', linewidth=0.5, alpha=0.9)
    ax3d.set_xlabel('X')
    ax3d.set_ylabel('Y')
    ax3d.set_zlabel('Z')
    ax3d.set_title('Surface of Revolution')

    # Equal aspect-ish
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0
    mid_x = (X.max()+X.min()) * 0.5
    mid_y = (Y.max()+Y.min()) * 0.5
    mid_z = (Z.max()+Z.min()) * 0.5
    ax3d.set_xlim(mid_x - max_range, mid_x + max_range)
    ax3d.set_ylim(mid_y - max_range, mid_y + max_range)
    ax3d.set_zlim(mid_z - max_range, mid_z + max_range)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    rotated_profile(
        filename_points='surface_points.csv',
        filename_meridian='meridian_curve.csv',
        y_min=1e-6, y_max=12-1e-6,  # avoid the exact endpoint at y=12 where base=0
        Ny=32, Ntheta=32
    )