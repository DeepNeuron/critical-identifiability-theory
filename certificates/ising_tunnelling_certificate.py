"""Reproducible numerical checks for the open Ising chain in Critical Identifiability Theory.

Run with Python 3, numpy and scipy. Frequencies use J=hbar=1.
"""
import numpy as np
from scipy.linalg import eigvalsh, svdvals
from scipy.optimize import brentq


def block(n, g, beta):
    d = np.full(n, float(g))
    d[0] *= beta
    return np.diag(d) - np.diag(np.ones(n-1), -1)


def slope(n, g):
    return 2*g**n / np.linalg.norm(g**np.arange(n))


def spin_gap(n, g, beta):
    # Work in the X product basis. ZZ flips neighboring bits, preserving parity.
    e = []
    for parity in (0, 1):
        basis = [i for i in range(2**n) if i.bit_count()%2 == parity]
        index = {s:i for i,s in enumerate(basis)}
        H = np.zeros((len(basis), len(basis)))
        for i,s in enumerate(basis):
            H[i,i] = -g*sum((beta if j == 0 else 1)*(1-2*((s>>j)&1))
                            for j in range(n))
            for j in range(n-1):
                H[i,index[s^(3<<j)]] = -1
        e.append(eigvalsh(H,subset_by_index=[0,0])[0])
    return e[1]-e[0]


def run():
    for n in (3, 5, 8):
        for g in (.62, 1., 1.15):
            B = block(n,g,0)
            zero = g**np.arange(n-1,-1,-1)
            assert np.linalg.norm(B@zero) < 1e-12
            bulk = 2*np.sqrt(1+g*g-2*g*np.cos(np.pi/n))
            assert np.isclose(2*svdvals(B)[-2], bulk)
            assert np.isclose(spin_gap(n,g,1e-4)/1e-4,slope(n,g),rtol=1e-6)
            for beta in (-.4, .4, 1.):
                assert np.isclose(spin_gap(n,g,beta),np.sign(beta)*2*svdvals(block(n,g,beta))[-1])
    print("Parity-resolved spin spectra, zero modes, slope and bulk gaps verified.")
    x=.8
    for n in (16,64,256):
        g=np.exp(-x/n)
        ratio=n*slope(n,g)**2/(4*(2*x/np.expm1(2*x)))
        print(f'N={n:3d} critical-window QFI / limit = {ratio:.7f}')
    b=.7
    y=brentq(lambda z:z*np.tan(z)-b*b,0,np.pi/2-1e-12)
    for n in (16,64,256):
        beta=b/np.sqrt(n)
        gap=2*svdvals(block(n,1.,beta))[-1]
        k=2*np.arcsin(gap/4)
        assert np.isclose(np.tan(n*k)*np.tan(k/2),beta**2/(2-beta**2))
        print(f'N={n:3d} N*gap/(2*y) = {n*gap/(2*y):.7f}')
    # Exact sum of the odd tunnelling sectors, compared to the zero-mode calculation.
    from math import comb
    n,g=8,.92
    q=g**n
    series=sum(comb(2*m,m)/4**m*q**(2*m+1) for m in range(50))
    assert np.isclose(2*np.sqrt(1-g*g)*series,slope(n,g))
    for g in (.72,.90):
        n,beta=35,.6
        approx=2*g**n*abs(beta)*(1-g*g)/np.sqrt(1-g*g+beta**2*g*g)
        assert np.isclose(2*svdvals(block(n,g,beta))[-1],approx,rtol=.02,atol=0)
    from scipy.integrate import quad
    from scipy.special import kv
    for z in (2.,8.):
        cut = quad(lambda x: np.exp(-z*x)*np.sqrt(x*x-1),1,np.inf)[0]
        assert np.isclose(cut,kv(1,z)/z)
        n=1000
        k=np.pi*np.arange(1,n)/n
        finite=1/n + (2/n)*np.sum(np.cos(k/2)**2*np.exp(-z*np.sin(k/2)))
        limit=(4/np.pi)*quad(lambda x: np.sqrt(1-x*x)*np.exp(-z*x),0,1)[0]
        assert abs(finite-limit)<3e-6
    print("Boundary spectral measure and Borel cut integral verified.")
    print("All checks passed.")


if __name__=='__main__':
    run()
