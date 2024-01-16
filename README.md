# ml-kem-stuff

A toy implementation of ML-KEM, aka Kyber, based on the current [FIPS 203 draft](https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.203.ipd.pdf).

This implementation is INCOMPLETE and **INSECURE**, the latter being explicitly out of scope. The current status is "it looks like it works" - it's able to derive a shared secret, but I have no idea if it's doing it correctly. ~~I'll do some more rigorous testing later.~~ It passes [these test vectors](https://github.com/post-quantum-cryptography/KAT/tree/main/MLKEM) at least.

This implementation is *NOT* constant-time, currently.

I wrote this to get my head around ML-KEM, and it's honestly a lot simpler than I was expecting. I'd also read [this](https://words.filippo.io/dispatches/kyber-math/) beforehand, which helped a lot. I was particularly surprised by how simple the implementation of the NTT is, given how conceptually mind-bending it is (at least, it is for me).

You may also be interested in [GiacomoPope/kyber-py](https://github.com/GiacomoPope/kyber-py), which is basically a more polished version of this repo. However, I avoided looking at any other implementations before embarking on my own, because I wanted to challenge myself and make sure I wasn't skipping any details.

## TODO:

- Proper input validation.

- ~~Actually test it against some test vectors.~~ ML-KEM-768 KAT passed!

- ~~Find out if the minor spec bugs I encountered are real and/or already reported.~~ (answer: yes to both!)