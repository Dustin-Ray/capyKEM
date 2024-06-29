## capyKEM - research into quantum-resistant algorithm design

This repo is a pure rust, no-std interpretation of [FIPS 203 (draft)](https://csrc.nist.gov/pubs/fips/203/ipd) which leverages a module learning-with-errors based construction aimed towards providing a secure means of key exchange when faced with a potential quantum adversary.

THIS LIBRARY IS A DRAFT AND IS NOT SAFE FOR USE. It exists purely for acedemic exeperimentation.

### Current working items:

- [ ] document in style of FIPS
- [ ] support other two parameter sets
- [ ] parameterize sample_poly_cbd over eta
- [x] condense encoding/decoding to single function
- [x] parameterize encoding/decoding over d
- [ ] implement API-level functions
- [ ] replace usage of ```Vec``` with ```hybridarray```


## Notions beyond IND-CCA:

Schmieg proves [here](https://eprint.iacr.org/2024/523) that misbinding properties can occur due to the way private keys are serialized and fixed by using a single seed to generate the private key, and thus ML-KEM-768 (generalized to other variants as well) is not MAL-BIND-K-CT or MAL-BIND-K-PK secure. This conclusion is drawn from this [paper](https://eprint.iacr.org/2023/1933) which introduces the MAL-BIND security notions which extend beyond IND-CCA.

NIST is now [proposing](https://groups.google.com/a/list.nist.gov/g/pqc-forum/c/5CT4NC_6zRI/m/lpifFrpWAwAJ?utm_medium=email&utm_source=footer) the following modification to the FIPS 203 IPD:

> "We propose ML-KEM uses a single 32-byte seed as decapsulation key, from which rho, sigma, and z are expanded.
> 
> This is smaller and simpler. Simpler, because we do not need to think about decapsulation key formatting or validation. In particular, it ensures that ML-KEM is MAL-BIND-K-CT and MAL-BIND-K-PK"
> 

## Proposed update to FIPS 203 IPD:

The IPD currently specifies that key expansion is unpacked before decapsulation, also precomputing $A$:

```
"packed"          unpack        ready-to-use                      keygen
decaps key:     -------->       decaps key:           <---------  64 byte seed
s, ek, H(ek), z                 s, ek, H(ek), z, A                d, z
```

NISTs proposal is to simplify this process:

```
"packed"       unpack=keygen        ready-to-use
decaps key:       -------->         decaps key:
g                                   s, ek, H(ek), z, A
```
With the caveat: 

> To be clear, we do *not* propose that FIPS 203 specifies this two-step approach. It merely should not preclude it.

## Implementation Routes:

Two possible means of applying the proposed update are suggested:

### Solution 1: brief:

1. Rename K-PKE.KeyGen to K-PKE.ExpandPrivate; remove lines 1, 2; and add rho and sigma as arguments.
2. Rename ML-KEM.KeyGen to ML-KEM.ExpandPrivate; add a 32-byte g as argument; call K-PKE.ExpandPrivate instead of K-PKE.KeyGen on line 2 passing rho and sigma; and replace line 1 by: $(\rho, \sigma, z) = J(g)$
3. Define a new ML-KEM.KeyGen as:

```
g <$- B^32
ek, dk = ML-KEM.ExpandPrivate(g)
return (ek, g)
```

4. Change ML-KEM.Decaps to take a 32-byte g as argument instead of dk. Add before line 1: 

```
_, dk = ML-KEM.ExpandPrivate(g)
```

### Solution 2: tight integration and readability, closer to implementations in practice:

1. Rename K-PKE.KeyGen to K-PKE.ExpandPrivate; remove lines 1, 2, 20, and 21; add rho and sigma as arguments; and return (^A, ^t, ^s).

2. Add a new function ML-KEM.UnpackPrivate that takes a 32-byte seed dk as argument, and acts as follows:

SHAKE-256:
$(\rho, \sigma, z) = J(dk)$    
```
^A, ^t, ^s = K-PKE.ExpandPrivate(rho, sigma)
ek = ByteEncode_12(^t) || rho
return (^A, ^t, ^s, ek, H(ek), z)
```

3. Change ML-KEM.KeyGen to:

```
dk <$- B^32
(^A, ^t, ^s, ek, h, z) = ML-KEM.UnpackPrivate(dk)
ek = ByteEncode_12(^t) || rho
return (ek, dk)
```

4. Change K-PKE.Encrypt to accept ^A, and ^t directly instead of ek_PKE, removing lines 2–8.
6. Change K-PKE.Decrypt to accept ^s as argument directly instead of dk_PKE, removing line 5.
7. Change ML-KEM.Decaps to accept the shortened 32-byte dk. Replace lines 1-4 by
```
(^A, ^t, ^s, ek, h, z) = ML-KEM.UnpackPrivate(dk)
```
and pass ^s instead of dk_PKE to K-PKE.Decrypt (line 7); and ^A, ^t instead of ek_PKE to K-PKE.Encrypt (line 8).

8. Change ML-KEM.Encaps to include the lines 2–8 removed from K-PKE.Encrypt before the call to K-PKE.Encrypt, to which ^A and ^t are passed instead of ek_PKE.

## (Plausible) Post-Quantum Secure Cryptosystem
Our larger [cryptographic algorithm library](https://github.com/drcapybara/capyCRYPT) pairs ML-KEM-768 to a SHA3-sponge construction for a quantum-safe public-key cryptosystem. It offers theoretic quantum-security through the use of the KEM and sponge primitives, which are both based on problems conjectured to be hard to solve for a quantum adversary. This design seeds the SHA-3 sponge with the secret shared through the KEM + a session nonce, which then faciliates high-performance symmetric encryption/decryption of arbitrary-length messages.

Our construction is non-standard, has not been subject to peer review, and lacks any formal audit. Our [ML-KEM library](https://github.com/drcapybara/capyKEM) itself is a work in progress and only supports the recommended NIST-II security parameter-set of 768. Furthermore, the current FIPS 203 IPD is, (as the name indicates), a draft, and final details about secure implementation may be subject to change. Our design currently exists in this library purely as an academic curiosity. Use it at your own risk, we provide no guarantee of security, reliability, or efficiency.

## Acknowledgements
Our [KEM](https://github.com/drcapybara/capyKEM) is inspired by the excellent ML-KEM articles and [go implementation](https://pkg.go.dev/filippo.io/mlkem768) by Filippo Valsorda and the always wonderful rust-crypto implementation by the great Tony Arcieri [here](https://crates.io/crates/ml-kem).
