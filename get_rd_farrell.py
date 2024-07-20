# Diffusion Theory
def get_rd_farrell(mua, musp, n):
    """
    Calculates the diffuse reflectance (Rd) from
    semi-infinite homogeneous medium.

    Args:
      mua: Absorption coefficient (cm^-1 or mm^-1).
      musp: Reduced scattering coefficient (cm^-1 or mm^-1), musp = mus(1-g),
            where mus is the scattering coefficient, and
            g is the anisotropy of scattering.
      n: Ratio of refractive index of object medium to outside medium.
      
    Returns:
      Rd: Diffusion coefficient (same units as mua and musp).

    Example:
       print(get_rd_farrell([.001, .01],[1, 2],1.4))
       
    Farrell TJ, MS Patterson, B Wilson, Medical Physics 19:879-888, 1992
    """

    mua=np.array(mua)
    musp=np.array(musp)
    
    ri = 0.6681 + 0.0636 * n + 0.7099 / n - 1.4399 / n**2
    A = (1 + ri) / (1 - ri)

    zo = 1 / (mua + musp)
    D = zo / 3
    delta = np.sqrt(D / mua)
    mueff = 1 / delta
    ap = musp / (mua + musp)

    Rd = ap * np.exp(-mueff * zo) / 2 * (1 + np.exp(-4/3 * A * np.sqrt(3 * (1 - ap))))
    return Rd

