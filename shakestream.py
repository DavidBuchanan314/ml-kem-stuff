# this is kinda gross. it exists to work around API limitations.
class ShakeStream:
	def __init__(self, digestfn) -> None:
		# digestfn is anything we can call repeatedly with different lengths
		self.digest = digestfn
		self.buf = self.digest(32) # arbitrary starting length
		self.offset = 0
	
	def read(self, n: int) -> bytes:
		# double the buffer size until we have enough
		while self.offset + n > len(self.buf):
			self.buf = self.digest(len(self.buf) * 2)
		res = self.buf[self.offset:self.offset + n]
		self.offset += n
		return res


if __name__ == "__main__":
	from hashlib import shake_128

	a = ShakeStream(shake_128(b"hello").digest)
	foo = a.read(17) + a.read(5) + a.read(57) + a.read(1432) + a.read(48)
	bar = shake_128(b"hello").digest(17 + 5 + 57 + 1432 + 48)
	assert(foo == bar)