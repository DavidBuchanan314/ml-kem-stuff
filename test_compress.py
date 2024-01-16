from fractions import Fraction
from mlkem import compress, decompress, Q

def round_half_up(n: Fraction) -> int:
	return (n + Fraction(1, 2)).__floor__()

def compress_proper(d: int, x: int) -> int:
	return round_half_up(Fraction(2**d * x, Q)) % (2**d)

def decompress_proper(d: int, y: int) -> int:
	return round_half_up(Fraction(Q * y, 2**d)) % Q

# https://github.com/pq-crystals/kyber/blob/b628ba78711bc28327dc7d2d5c074a00f061884e/ref/poly.c#L191-L210
def compress_kyber_ref_poly_tomsg(t: int) -> int:
	t <<= 1
	t += 1665
	t *= 80635
	t >>= 28
	t &= 1
	return t

# https://github.com/pq-crystals/kyber/blob/b628ba78711bc28327dc7d2d5c074a00f061884e/ref/poly.c#L32-L36
def compress_kyber_ref_poly_128(t: int) -> int:
	t <<= 4
	t += 1665
	t *= 80635
	t >>= 28
	t &= 0xf
	return t

# https://github.com/pq-crystals/kyber/blob/b628ba78711bc28327dc7d2d5c074a00f061884e/ref/poly.c#L52-L56
def compress_kyber_ref_poly_160(t: int) -> int:
	t <<= 5
	t += 1664
	t *= 40318
	t >>= 27
	t &= 0x1f
	return t

# https://github.com/pq-crystals/kyber/blob/b628ba78711bc28327dc7d2d5c074a00f061884e/ref/polyvec.c#L28-L33
def compress_kyber_ref_polyvec_128(t: int) -> int:
	t <<= 11
	t += 1664
	t *= 645084
	t >>= 31
	t &= 0x7ff
	return t

# https://github.com/pq-crystals/kyber/blob/b628ba78711bc28327dc7d2d5c074a00f061884e/ref/polyvec.c#L58-L63
def compress_kyber_ref_polyvec_160(t: int) -> int:
	t <<= 10
	t += 1665
	t *= 1290167
	t >>= 32
	t &= 0x3ff
	return t

for d in range(12):
	for n in range(Q):
		assert(compress_proper(d, n) == compress(d, [n])[0])

for d in range(12):
	for n in range(2**d):
		assert(decompress_proper(d, n) == decompress(d, [n])[0])

for n in range(Q):
	assert(compress_proper(1, n) == compress_kyber_ref_poly_tomsg(n))
	assert(compress_proper(4, n) == compress_kyber_ref_poly_128(n))
	assert(compress_proper(5, n) == compress_kyber_ref_poly_160(n))
	assert(compress_proper(11, n) == compress_kyber_ref_polyvec_128(n))
	assert(compress_proper(10, n) == compress_kyber_ref_polyvec_160(n))