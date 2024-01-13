/* eslint-disable @typescript-eslint/no-unused-vars */

import { strict as assert } from 'assert';

// ML-KEM-768 params:
const N = 256;
const Q = 3329;
const K = 3;
const ETA1 = 2;
const ETA2 = 2;
const DU = 10;
const DV = 4;

// ring of integers mod Q
class Zq {
	readonly n: number;

	constructor(n: number) {
		assert(n >= 0);
		assert(n < Q);
		assert(Number.isInteger(n));
		this.n = n;
	}

	// converts numbers to Zq instances, leaves existing Zq instances untouched
	static cast(n: Zq | number): Zq {
		if (n instanceof Zq) return n;
		return new Zq(n);
	}

	add(other: Zq | number): Zq {
		return new Zq((this.n + Zq.cast(other).n) % Q);
	}

	sub(other: Zq | number): Zq {
		return new Zq((Q + this.n - Zq.cast(other).n) % Q);
	}

	// Zq elements should be 12-bit values, yielding 24-bit multiplication results, prior to reduction.
	// These values should all fit in JS's fp64 numbers, without precision issues.
	mul(other: Zq | number): Zq {
		return new Zq((this.n * Zq.cast(other).n) % Q);
	}

	// XXX: this is NOT constant-time, but it's also only used in LUTs.
	// Implements the square-and-multiply algorithm
	pow(exp: Zq | number): Zq {
		const expq = Zq.cast(exp);
		let val = new Zq(1); // "multiplicative identity element" (!)
		for (let i=11; i>=0; i--) { // read bits MSB to LSB. 11 is floor(log2(Q))
			val = val.mul(val); // square
			if ((expq.n >> i) & 1) {
				val = val.mul(this); // conditional multiply
			}
		}
		return val;
	}
}

// Polynomial ring Zq[X]/(X^n + 1)
class Rq {
	readonly f: Zq[];

	constructor(f: Zq[]) {
		assert(f.length === N);
		this.f = f;
	}

	static cast(n: Rq | Zq[] | number[]): Rq {
		if (n instanceof Rq) return n;
		return new Rq(n.map(Zq.cast));
	}

	add(other: Rq): Rq {
		return new Rq(this.f.map((e, i) => e.add(other.f[i])))
	}

	sub(other: Rq): Rq {
		return new Rq(this.f.map((e, i) => e.sub(other.f[i])))
	}

	// this is the slow textbook multiplication algorithm
	// it exists for testing purposes and isn't used in the
	// ML-KEM impl proper.
	mul(other: Rq): Rq {
		const c: Zq[] = Array(511).fill(new Zq(0)); // Zq(0) is additive identity element
		
		// multiplication loop
		for (let i=0; i<256; i++) {
			for (let j=0; j<256; j++) {
				c[i+j] = c[i+j].add(this.f[j].mul(other.f[i]));
			}
		}

		// reduction loop
		for (let i=0; i<255; i++) {
			c[i] = c[i].sub(c[i+256]);
			// c[i+256] is implicitly zeroed
		}

		return new Rq(c.slice(0, 256));
	}

	ntt(): Tq {
		const f = [...this.f];
		let k = 1;
		for (let len=128; len>=2; len/=2) {
			for (let start=0; start<256; start += 2*len) {
				const zeta = ZETA[k++];
				for (let j=start; j<start+len; j++) {
					const t = zeta.mul(f[j + len]);
					f[j + len] = f[j].sub(t);
					f[j] = f[j].add(t);
				}
			}
		}
		return new Tq(f);
	}
}

// NTT-transformed Rq
class Tq {
	readonly f: Zq[]; // stored as a 1d array, but represents an array of pairs

	constructor(f: Zq[]) {
		assert(f.length === N);
		this.f = f;
	}

	static cast(n: Tq | Zq[] | number[]): Tq {
		if (n instanceof Tq) return n;
		return new Tq(n.map(Zq.cast));
	}

	add(other: Tq): Tq {
		return new Tq(this.f.map((e, i) => e.add(other.f[i])))
	}

	sub(other: Tq): Tq {
		return new Tq(this.f.map((e, i) => e.sub(other.f[i])))
	}

	// in the NTT domain, multiplication is relatively trivial
	mul(other: Tq): Tq {
		const c: Zq[] = [];
		for (let i=0; i<256; i += 2) {
			const [a0, a1] = this.f.slice(i, i+2);
			const [b0, b1] = other.f.slice(i, i+2);
			c.push(a0.mul(b0).add(a1.mul(b1).mul(GAMMA[i/2])));
			c.push(a0.mul(b1).add(a1.mul(b0)));
		}
		return new Tq(c);
	}

	nttInv(): Rq {
		const f = [...this.f];
		let k = 127;
		for (let len=2; len<=128; len*=2) {
			for (let start=0; start<256; start += 2*len) {
				const zeta = ZETA[k--];
				for (let j=start; j<start+len; j++) {
					const t = f[j];
					f[j] = t.add(f[j + len]);
					f[j + len] = zeta.mul(f[j + len].sub(t));
				}
			}
		}
	
		for (let i=0; i<256; i++) {
			f[i] = f[i].mul(3303); // 3303 === powQ(128, Q-2)
		}
	
		return new Rq(f);
	}
}


// doesn't matter if this is slow, it gets rolled into LUTs anyway
function bitRev7(n: number): number {
	assert((n >= 0) && (n < (1<<7)));

	// 8-bit rev bithack
	n = ((n & 0b10101010) >> 1) | ((n & 0b01010101) << 1);
	n = ((n & 0b11001100) >> 2) | ((n & 0b00110011) << 2);
	n = ((n & 0b11110000) >> 4) | ((n & 0b00001111) << 4);

	return n >> 1; // shift off 0th bit (giving us a 7-bit rev overall)
}

// indexing these LUTs won't be a sidechannel because the index is not secret
const ZETA: Zq[] = [];
const GAMMA: Zq[] = [];
const UNITY_ROOT = new Zq(17);

// initialise ZETA and GAMMA LUTs
for (let k=0; k<128; k++) {
	ZETA.push(UNITY_ROOT.pow(bitRev7(k)));
	GAMMA.push(UNITY_ROOT.pow(2*bitRev7(k)+1));
}


// conversions

function bitsToBytes(bits: number[]): Uint8Array {
	assert((bits.length % 8) === 0);
	const bytes: number[] = [];
	for (let i=0; i<bits.length; i+=8) {
		let byte = 0;
		for (let j=0; j<8; j++) {
			byte |= bits[i+j] << j;
		}
		bytes.push(byte);
	}
	return new Uint8Array(bytes);
}

function bytesToBits(bytes: Uint8Array): number[] {
	const bits: number[] = [];
	for (const byte of bytes) {
		for (let i=0; i<8; i++) {
			bits.push((byte >> i) & 1);
		}
	}
	return bits;
}

function byteEncode(d: number, f: number[]): Uint8Array {
	assert(f.length === N);
	const bits: number[] = [];
	for (const a of f) {
		for (let i=0; i<d; i++) {
			bits.push((a >> i) & 1);
		}
	}
	return bitsToBytes(bits);
}

function byteDecode(d: number, bytes: Uint8Array): number[] {
	assert(bytes.length === d*N/8);
	const bits = bytesToBits(bytes);
	const values: number[] = [];
	for (let i=0; i<bits.length; i += d) {
		let value = 0;
		for (let j=0; j<d; j++) {
			value |= bits[i + j] << j
		}
		values.push(value);
	}
	return values;
}

const mynumbers = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255];
const foo = Rq.cast(mynumbers);
const bar = foo.add(foo);

console.log(foo.mul(bar));

const foontt = foo.ntt();
const barntt = bar.ntt();
const nttres = foontt.mul(barntt);
console.log(nttres.nttInv())

/*
console.log(nttres);
console.log(nttres.nttInv());
const mybytes = byteEncode(12, mynumbers);
console.log(mybytes);
console.log(byteDecode(12, mybytes))

console.log(new Zq(123));*/