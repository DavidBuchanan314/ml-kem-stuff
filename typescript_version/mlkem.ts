/* eslint-disable @typescript-eslint/no-unused-vars */

import { strict as assert } from 'assert';
import { HashOptions, createHash, randomBytes } from 'crypto';


// "global" ML-KEM params:
const N = 256;
const Q = 3329;

type MLKemParams = {
	k: number,
	eta1: number,
	eta2: number,
	du: number,
	dv: number,
};

// the params that actually vary
const MLKEM_PARAMS_512:  MLKemParams = { k: 2, eta1: 3, eta2: 2, du: 10, dv: 4 };
const MLKEM_PARAMS_768:  MLKemParams = { k: 3, eta1: 2, eta2: 2, du: 10, dv: 4 };
const MLKEM_PARAMS_1024: MLKemParams = { k: 4, eta1: 2, eta2: 2, du: 11, dv: 5 };


// doesn't matter if this is slow, it gets rolled into LUTs anyway
function bitRev7(n: number): number {
	assert((n >= 0) && (n < (1<<7)));
	return parseInt(n.toString(2).split("").reverse().join("").padEnd(7, "0"), 2);
}


// hashes (FIPS 203, Section 4.1)

function hashPRF(eta: number, data: Uint8Array, b: number): Uint8Array {
	const hash = createHash("shake256", {outputLength: 64 * eta});
	hash.update(data);
	hash.update(new Uint8Array([b]));
	return hash.digest();
}

type XOFStream = {
	read: (length: number) => Uint8Array
};

function hashXOF(data: Uint8Array, i: number, j: number): XOFStream {
	const hash = createHash("shake128");
	hash.update(data);
	hash.update(new Uint8Array([i, j]));
	let buffer = hash.copy().digest();
	let offset = 0;
	return {
		read: (length: number) => {
			// since there isn't a "proper" streaming API, we double the buffer length
			// each time we run out of data, and extract from the current offset.
			while (buffer.length < offset + length) {
				// `as HashOptions` is a workaround, hash.copy() isn't typed properly
				buffer = hash.copy({outputLength: buffer.length * 2} as HashOptions).digest();
			}
			const res = buffer.subarray(offset, offset += length);
			// APIs that can return less than the requested number of bytes are scary
			assert(res.length == length);
			return res;
		}
	}
}

function hashH(data: Uint8Array): Uint8Array {
	const hash = createHash("sha3-256");
	hash.update(data);
	return hash.digest();
}

function hashJ(data: Uint8Array): Uint8Array {
	const hash = createHash("shake256", {outputLength: 32});
	hash.update(data);
	return hash.digest();
}

function hashG(data: Uint8Array): Uint8Array {
	const hash = createHash("sha3-512");
	hash.update(data);
	return hash.digest();
}


// conversions

// FIPS 203 Algorithm 2
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

// FIPS 203 Algorithm 3
function bytesToBits(bytes: Uint8Array): number[] {
	const bits: number[] = [];
	for (const byte of bytes) {
		for (let i=0; i<8; i++) {
			bits.push((byte >> i) & 1);
		}
	}
	return bits;
}

// FIPS 203 Algorithm 4
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

// FIPS 203 Algorithm 5
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

// indexing these LUTs won't be a sidechannel because the index is not secret
const ZETA: Zq[] = [];
const GAMMA: Zq[] = [];
const UNITY_ROOT = new Zq(17);

// initialise ZETA and GAMMA LUTs
for (let k=0; k<128; k++) {
	ZETA.push(UNITY_ROOT.pow(bitRev7(k)));
	GAMMA.push(UNITY_ROOT.pow(2*bitRev7(k)+1));
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

	// Sample from the distribution Dη(Rq)
	// FIPS 203 Algorithm 7
	static sampleCBD(eta: number, data: Uint8Array): Rq {
		assert(data.length == 64 * eta);
		const bits = bytesToBits(data);
		const f: Zq[] = [];
		for (let i=0; i<256*2; i+=2) {
			const x = new Zq(bits.slice((i+0)*eta, (i+1)*eta).reduce((a,b)=>a+b));
			const y = new Zq(bits.slice((i+1)*eta, (i+2)*eta).reduce((a,b)=>a+b));
			f.push(x.sub(y));
		}
		return new Rq(f);
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

	// FIPS 203 Algorithm 8
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

	// sample a uniform random element of Tq, from an XOF output stream
	// FIPS 203, Algorithm 6
	static sample(b: XOFStream): Tq {
		const a: Zq[] = [];
		while (a.length < 256) {
			const buf = b.read(3);
			const d1 = buf[0] | ((buf[1] & 0x0f) << 8);
			const d2 = (buf[1] >> 4) | (buf[2] << 4);
			if (d1 < Q) {
				a.push(new Zq(d1));
			}
			if (d2 < Q && a.length < 256) {
				a.push(new Zq(d2));
			}
		}
		return new Tq(a);
	}

	add(other: Tq): Tq {
		return new Tq(this.f.map((e, i) => e.add(other.f[i])))
	}

	sub(other: Tq): Tq {
		return new Tq(this.f.map((e, i) => e.sub(other.f[i])))
	}

	// in the NTT domain, multiplication is relatively trivial
	// FIPS 203 Algorithm 10, 11
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

	// FIPS 203 Algorithm 9
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
			f[i] = f[i].mul(3303); // 3303 === 128 ^ -1 mod Q
		}
	
		return new Rq(f);
	}

	toBytes(): Uint8Array {
		return byteEncode(12, this.f.map(x=>x.n));
	}
}

function KPKEKeyGen(
	params: MLKemParams,
	seed: Uint8Array=new Uint8Array(randomBytes(32))
): [Uint8Array, Uint8Array]
{
	const ghash = hashG(seed);
	const rho = ghash.subarray(0, 32);
	const sigma = ghash.subarray(32, 64);

	const a: Tq[][] = [];
	const s: Tq[] = [];
	const e: Tq[] = [];
	const t: Tq[] = [];

	for (let i=0; i<params.k; i++) {
		a.push([]);
		for (let j=0; j<params.k; j++) {
			a[i].push(Tq.sample(hashXOF(rho, i, j)));
		}

		s.push(Rq.sampleCBD(params.eta1, hashPRF(params.eta1, sigma, i)).ntt());
		e.push(Rq.sampleCBD(params.eta1, hashPRF(params.eta1, sigma, params.k + i)).ntt());
	}

	for (let i=0; i<params.k; i++) {
		let acc = e[i];
		for (let j=0; j<params.k; j++) {
			acc = acc.add(a[j][i].mul(s[j]));
		}
		t.push(acc);
	}

	const ekpke = Buffer.concat([Buffer.concat(t.map(x=>x.toBytes())), rho]);
	const dkpke = Buffer.concat(s.map(x=>x.toBytes()))
	
	return [new Uint8Array(ekpke), new Uint8Array(dkpke)];
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

const blah = hashXOF(new Uint8Array([0x41]), 0, 0);
console.log(blah.read(32));
console.log(blah.read(5));
console.log(blah.read(72));

console.log("foo")
console.log(KPKEKeyGen(MLKEM_PARAMS_768))
console.log("abc")
console.log(KPKEKeyGen(MLKEM_PARAMS_768))