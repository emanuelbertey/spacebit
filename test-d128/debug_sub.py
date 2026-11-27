#!/usr/bin/env python3
"""Simula en Python el algoritmo de subanillos generado en poly_qbig_sub.c."""
D = 128
S = 4
K = D // S
q = 134215681


def root(order):
    e = (q - 1) // order
    for g in range(2, 1000000):
        r = pow(g, e, q)
        if pow(r, order, q) == 1 and pow(r, order // 2, q) != 1:
            return r


eta = root(2 * K)
zeta = eta * eta % q


def bitrev(a):
    n = len(a)
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j ^= bit
        if i < j:
            a[i], a[j] = a[j], a[i]


def ct(a, w, n):
    a = list(a)
    bitrev(a)
    length = 2
    while length <= n:
        wl = pow(w, n // length, q)
        for i in range(0, n, length):
            ww = 1
            for j in range(length // 2):
                u = a[i + j]
                v = a[i + j + length // 2] * ww % q
                a[i + j] = (u + v) % q
                a[i + j + length // 2] = (u - v) % q
                ww = ww * wl % q
        length <<= 1
    return a


def ntt_sub(a):
    out = [0] * D
    for r in range(S):
        pw = 1
        tmp = []
        for qi in range(K):
            tmp.append(a[qi * S + r] * pw % q)
            pw = pw * eta % q
        tmp = ct(tmp, zeta, K)
        for i in range(K):
            out[r * K + i] = tmp[i]
    return out


def intt_sub(A):
    out = [0] * D
    zinv = pow(zeta, q - 2, q)
    etainv = pow(eta, q - 2, q)
    kinv = pow(K, q - 2, q)
    for r in range(S):
        tmp = [A[r * K + i] % q for i in range(K)]
        tmp = ct(tmp, zinv, K)
        pw = 1
        for qi in range(K):
            tmp[qi] = tmp[qi] * kinv % q * pw % q
            pw = pw * etainv % q
        for qi in range(K):
            out[qi * S + r] = tmp[qi]
    return out


def mult_sub(A, B):
    out = [0] * D
    for i in range(K):
        beta = pow(eta, 2 * i + 1, q)
        ca = [A[r * K + i] % q for r in range(S)]
        cb = [B[r * K + i] % q for r in range(S)]
        co = [0] * S
        for j in range(S):
            for k in range(S):
                prod = ca[j] * cb[k] % q
                if j + k < S:
                    co[j + k] = (co[j + k] + prod) % q
                else:
                    co[j + k - S] = (co[j + k - S] + beta * prod) % q
        for r in range(S):
            out[r * K + i] = co[r]
    return out


a = [(i * 99991 + 17) % q for i in range(D)]
b = [(i * 77771 + 41) % q for i in range(D)]
print("eta orden 2K:", pow(eta, 2 * K, q) == 1 and pow(eta, K, q) != 1)
print("roundtrip sub:", intt_sub(ntt_sub(a)) == a)
ref = [0] * D
for i in range(D):
    for j in range(D):
        p = a[i] * b[j] % q
        k = i + j
        ref[k if k < D else k - D] = (
            ref[k if k < D else k - D] + (p if k < D else q - p)) % q
got = intt_sub(mult_sub(ntt_sub(a), ntt_sub(b)))
print("convolucion sub:", got == ref,
      [(i, got[i], ref[i]) for i in range(D) if got[i] != ref[i]][:3])
