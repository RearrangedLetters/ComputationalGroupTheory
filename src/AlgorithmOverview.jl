""" Overview over the implemented algorithms and open things.

Symbols:
    • (I) : implemented
    • (P) : partially implemented
    • (N) : not implemented
    • (O) : optimization known
    • (T) : tested
    • (A) : indirectly available, i.e. with slight modification/extension

Chapter 01 & 02 - Orbits and Stabilizers
    • Orbit (I)
    • Transversal
        ◀ Orbit-Transversal (I)
        ◀ Orbit-Factored-Transversal (I)
        ◀ Factored-Transversal-Reconstruction (I)
    • Normal Closure (I)
    • Commutator Subgroup (I)
    • Pseudorandom Elements aka Product Replacement (I)

Chapter 03 - Stabilizer Chains
    • Sift - write element as product of coset representatives (P)
    • Membership Test (I): We have sift(g) == e iff g ∈ G
    • Schreier Sims (I)
    • Order (I)
    • Derived Series (A)
    • Lower Central Series (A)
    • Elements in same coset (A)
    • Determine permutation action on the coset of a subgroup (A)
    • Pointwise Stabilizer (A)
    • Enumerate g (A)
    • Homomorphisms

Chapter 04 - Backtracking
    • Plain Backtrack (P)
    • Setwise Stabilizer (N)
    • Find conjugating element (N)
    • Centralizer (N)
    • Normalizer (N)

"""