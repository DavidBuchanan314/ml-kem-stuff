# https://words.filippo.io/dispatches/kyber-math/
# https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.203.ipd.pdf

from typing import List, Tuple
import hashlib
import os
from shakestream import ShakeStream
from functools import reduce

# ML-KEM-768 params:
N = 256
Q = 3329
K = 3
ETA1 = 2
ETA2 = 2
DU = 10
DV = 4

def bitrev7(n: int) -> int:
	return int(f"{n:07b}"[::-1], 2)  # gross but it works

# 17 is primitive 256th root of unity mod Q
ZETA = [pow(17, bitrev7(k), Q) for k in range(128)] # used in ntt and ntt_inv
GAMMA = [pow(17, 2*bitrev7(k)+1, Q) for k in range(128)] # used in ntt_mul

# can be reused for NTT representatives
def poly256_add(a: List[int], b: List[int]) -> List[int]:
	return [(x + y) % Q for x, y in zip(a, b)]

def poly256_sub(a: List[int], b: List[int]) -> List[int]:
	return [(x - y) % Q for x, y in zip(a, b)]

# naive O(n^2) multiplication algorithm for testing/comparison purposes.
# this is not used for the main impl.
def poly256_slow_mul(a: List[int], b: List[int]) -> List[int]:
	c = [0] * 511

	# textbook multiplication, without carry
	for i in range(256):
		for j in range(256):
			c[i+j] = (c[i+j] + a[j] * b[i]) % Q

	# now for reduction mod X^256 + 1
	for i in range(255):
		c[i] = (c[i] - c[i+256]) % Q
		# we could explicitly zero c[i+256] here, but there's no need...
	
	# ...because we're about to truncate c
	return c[:256]


# by the way, this is O(n logn)
def ntt(f_in: List[int]) -> List[int]:
	f_out = f_in.copy()
	k = 1
	for log2len in range(7, 0, -1):
		length = 2**log2len
		for start in range(0, 256, 2 * length):
			zeta = ZETA[k]
			k += 1
			for j in range(start, start + length):
				t = (zeta * f_out[j + length]) % Q
				f_out[j + length] = (f_out[j] - t) % Q
				f_out[j] = (f_out[j] + t) % Q
	return f_out


# so is this
def ntt_inv(f_in: List[int]) -> List[int]:
	f_out = f_in.copy()
	k = 127
	for log2len in range(1, 8):
		length = 2**log2len
		for start in range(0, 256, 2 * length):
			zeta = ZETA[k]
			k -= 1
			for j in range(start, start + length):
				t = f_out[j]
				f_out[j] = (t + f_out[j + length]) % Q
				f_out[j + length] = (zeta * (f_out[j + length] - t)) % Q

	for i in range(256):
		f_out[i] = (f_out[i] * 3303) % Q  # 3303 == pow(128, -1, Q)

	return f_out

ntt_add = poly256_add  # it's just elementwise addition

# and this is just O(n)
def ntt_mul(a: List[int], b: List[int]) -> List[int]:
	c = []
	for i in range(128):
		a0, a1 = a[2 * i: 2 * i + 2]
		b0, b1 = b[2 * i: 2 * i + 2]
		c.append((a0 * b0 + a1 * b1 * GAMMA[i]) % Q)
		c.append((a0 * b1 + a1 * b0) % Q)
	return c


# crypto functions

def mlkem_prf(eta: int, data: bytes, b: int) -> bytes:
	return hashlib.shake_256(data + bytes([b])).digest(64 * eta)

def mlkem_xof(data: bytes, i: int, j: int) -> ShakeStream:
	return ShakeStream(hashlib.shake_128(data + bytes([i, j])).digest)

def mlkem_hash_H(data: bytes) -> bytes:
	return hashlib.sha3_256(data).digest()

def mlkem_hash_J(data: bytes) -> bytes:
	return hashlib.shake_256(data).digest(32)

def mlkem_hash_G(data: bytes) -> bytes:
	return hashlib.sha3_512(data).digest()


# encode/decode logic

def bits_to_bytes(bits: List[int]) -> bytes:
	assert(len(bits) % 8 == 0)
	return bytes(
		sum(bits[i + j] << j for j in range(8))
		for i in range(0, len(bits), 8)
	)

def bytes_to_bits(data: bytes) -> List[int]:
	bits = []
	for word in data:
		for i in range(8):
			bits.append((word >> i) & 1)
	return bits

def byte_encode(d: int, f: List[int]) -> bytes:
	assert(len(f) == 256)
	bits = []
	for a in f:
		for i in range(d):
			bits.append((a >> i) & 1)
	return bits_to_bytes(bits)

def byte_decode(d: int, data: bytes) -> List[int]:
	bits = bytes_to_bits(data)
	return [sum(bits[i * d + j] << j for j in range(d)) for i in range(256)]

def compress(d: int, x: List[int]) -> List[int]:
	return [(((n * 2**d) + Q // 2 ) // Q) % (2**d) for n in x]

def decompress(d: int, x: List[int]) -> List[int]:
	return [(((n * Q) + 2**(d-1) ) // 2**d) % Q for n in x]


# sampling

def sample_ntt(xof: ShakeStream):
	res = []
	while len(res) < 256:
		a, b, c = xof.read(3)
		d1 = ((b & 0xf) << 8) | a
		d2 = c << 4 | b >> 4
		if d1 < Q:
			res.append(d1)
		if d2 < Q and len(res) < 256:
			res.append(d2)
	return res


def sample_poly_cbd(eta: int, data: bytes) -> List[int]:
	assert(len(data) == 64 * eta)
	bits = bytes_to_bits(data)
	f = []
	for i in range(256):
		x = sum(bits[2*i*eta+j] for j in range(eta))
		y = sum(bits[2*i*eta+eta+j] for j in range(eta))
		f.append((x - y) % Q)
	return f


# K-PKE

def kpke_keygen(seed: bytes=None) -> Tuple[bytes, bytes]:
	d = os.urandom(32) if seed is None else seed
	ghash = mlkem_hash_G(d)
	rho, sigma = ghash[:32], ghash[32:]

	ahat = []
	for i in range(K):
		row = []
		for j in range(K):
			row.append(sample_ntt(mlkem_xof(rho, i, j)))
		ahat.append(row)
	
	shat = [
		ntt(sample_poly_cbd(ETA1, mlkem_prf(ETA1, sigma, i)))
		for i in range(K)
	]
	ehat = [
		ntt(sample_poly_cbd(ETA1, mlkem_prf(ETA1, sigma, i+K)))
		for i in range(K)
	]
	that = [ # t = a * s + e
		reduce(ntt_add, [
			ntt_mul(ahat[j][i], shat[j])
			for j in range(K)
		] + [ehat[i]])
		for i in range(K)
	]
	ek_pke = b"".join(byte_encode(12, s) for s in that) + rho
	dk_pke = b"".join(byte_encode(12, s) for s in shat)
	return ek_pke, dk_pke


def kpke_encrypt(ek_pke: bytes, m: bytes, r: bytes) -> bytes:
	that = [byte_decode(12, ek_pke[i*128*K:(i+1)*128*K]) for i in range(K)]
	rho = ek_pke[-32:]

	# this is identical to as in kpke_keygen
	ahat = []
	for i in range(K):
		row = []
		for j in range(K):
			row.append(sample_ntt(mlkem_xof(rho, i, j)))
		ahat.append(row)
	
	rhat = [
		ntt(sample_poly_cbd(ETA1, mlkem_prf(ETA1, r, i)))
		for i in range(K)
	]
	e1 = [
		sample_poly_cbd(ETA2, mlkem_prf(ETA2, r, i+K))
		for i in range(K)
	]
	e2 = sample_poly_cbd(ETA2, mlkem_prf(ETA2, r, 2*K))

	u = [ # u = ntt-1(AT*r)+e1
		poly256_add(ntt_inv(reduce(ntt_add, [
			ntt_mul(ahat[i][j], rhat[j]) # note that i,j are reversed here
			for j in range(K)
		])), e1[i])
		for i in range(K)
	]
	mu = decompress(1, byte_decode(1, m))
	v = poly256_add(ntt_inv(reduce(ntt_add, [
		ntt_mul(that[i], rhat[i])
		for i in range(K)
	])), poly256_add(e2, mu))

	c1 = b"".join(byte_encode(DU, compress(DU, u[i])) for i in range(K))
	c2 = byte_encode(DV, compress(DV, v))
	return c1 + c2


def kpke_decrypt(dk_pke: bytes, c: bytes) -> bytes:
	c1 = c[:32*DU*K]
	c2 = c[32*DU*K:]
	u = [
		decompress(DU, byte_decode(DU, c1[i*32*DU:(i+1)*32*DU]))
		for i in range(K)
	]
	v = decompress(DV, byte_decode(DV, c2))
	shat = [byte_decode(12, dk_pke[i*384:(i+1)*384]) for i in range(K)]
	# NOTE: the comment in FIPS203 seems wrong here?
	# it says "NTT−1 and NTT invoked k times", but I think NTT−1 is only invoked once.
	w = poly256_sub(v, ntt_inv(reduce(ntt_add, [
		ntt_mul(shat[i], ntt(u[i]))
		for i in range(K)
	])))
	m = byte_encode(1, compress(1, w))
	return m


# KEM time

def mlkem_keygen(seed1=None, seed2=None):
	z = os.urandom(32) if seed1 is None else seed1
	ek_pke, dk_pke = kpke_keygen(seed2)
	ek = ek_pke
	dk = dk_pke + ek + mlkem_hash_H(ek) + z
	return ek, dk


def mlkem_encaps(ek: bytes, seed=None) -> Tuple[bytes, bytes]:
	# TODO !!!! input validation !!!!!!!
	m = os.urandom(32) if seed is None else seed
	ghash = mlkem_hash_G(m + mlkem_hash_H(ek))
	k = ghash[:32]
	r = ghash[32:]
	c = kpke_encrypt(ek, m, r)
	return k, c


def mlkem_decaps(c: bytes, dk: bytes) -> bytes:
	# TODO !!!! input validation !!!!!!!
	dk_pke = dk[:384*K]
	ek_pke = dk[384*K : 768*K + 32]
	h = dk[768*K + 32 : 768*K + 64]
	z = dk[768*K + 64 : 768*K + 96]
	mdash = kpke_decrypt(dk_pke, c)
	ghash = mlkem_hash_G(mdash + h)
	kdash = ghash[:32]
	rdash = ghash[32:]
	# NOTE: J() has unnecessary second argument in the spec???
	kbar = mlkem_hash_J(z + c)
	cdash = kpke_encrypt(ek_pke, mdash, rdash)
	if cdash != c:
		# I suppose this branch ought to be constant-time, but that's already out the window with this impl
		#print("did not match") # XXX: what does implicit reject mean? I suppose it guarantees it fails in a not-attacker-controlled way?
		return kbar
	return kdash


if __name__ == "__main__":
	a = list(range(256))
	b = list(range(1024, 1024+256))

	ntt_res = ntt_inv(ntt_add(ntt(a), ntt(b)))
	poly_res = poly256_add(a, b)

	assert(ntt_res == poly_res)

	ntt_prod = ntt_inv(ntt_mul(ntt(a), ntt(b)))
	poly_prod = poly256_slow_mul(a, b)

	assert(ntt_prod == poly_prod)


	ek_pke, dk_pke = kpke_keygen(b"SEED"*8)

	msg = b"This is a demonstration message."
	ct = kpke_encrypt(ek_pke, msg, b"RAND"*8)
	pt = kpke_decrypt(dk_pke, ct)
	print(pt)
	assert(pt == msg)


	ek, dk = mlkem_keygen()
	k1, c = mlkem_encaps(ek)
	print("encapsulated:", k1.hex())

	k2 = mlkem_decaps(c, dk)
	print("decapsulated:", k2.hex())

	assert(k1 == k2)