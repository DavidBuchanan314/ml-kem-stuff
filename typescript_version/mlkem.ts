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

// doesn't matter if this is slow, it gets rolled into LUTs anyway
function bitRev7(n: number): number {
	assert((n >= 0) && (n < (1<<7)));

	// 8-bit rev bithack
	n = ((n & 0xaa) >> 1) | ((n & 0x55) << 1);
	n = ((n & 0xcc) >> 2) | ((n & 0x33) << 2);
	n = ((n & 0xf0) >> 4) | ((n & 0x0f) << 4);

	return n >> 1; // shift off 0th bit (giving us a 7-bit rev overall)
}

// raise a number to a power, mod Q, using square-and-multiply
// this is NOT constant-time, but it's also only used in LUTs
function powQ(base: number, exp: number): number {
	// TODO: range checks on args?
	// `base` and `exp` should be 12-bit values, yielding 24-bit multiplication results.
	// these values should all fit in JS's fp64 numbers, without precision issues
	let val = 1; // "multiplicative identity element" (!)
	for (let i=11; i>=0; i--) { // read bits MSB to LSB
		val = (val * val) % Q; // square
		if ((exp >> i) & 1) {
			val = (val * base) % Q; // conditional multiply
		}
	}
	return val;
}

// indexing these LUTs won't be a sidechannel because the index is not secret
const ZETA: number[] = [];
const GAMMA: number[] = [];

// initialise ZETA and GAMMA LUTs
for (let k=0; k<128; k++) {
	ZETA.push(powQ(17, bitRev7(k)));
	GAMMA.push(powQ(17, 2*bitRev7(k)+1));
}

// I checked these against my python impl
//console.log(JSON.stringify(ZETA));
//console.log(JSON.stringify(GAMMA));


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


// NTT stuff

function ntt(f_in: number[]): number[] {
	assert(f_in.length === N);
	const f = [...f_in];
	let k = 1;
	for (let len=128; len>=2; len/=2) {
		for (let start=0; start<256; start += 2*len) {
			const zeta = ZETA[k++];
			for (let j=start; j<start+len; j++) {
				const t = (zeta * f[j + len]) % Q
				f[j + len] = (Q + f[j] - t) % Q
				f[j] = (f[j] + t) % Q
			}
		}
	}
	return f;
}

function nttInv(f_in: number[]): number[] {
	assert(f_in.length === N);
	const f = [...f_in];
	let k = 127;
	for (let len=2; len<=128; len*=2) {
		for (let start=0; start<256; start += 2*len) {
			const zeta = ZETA[k--];
			for (let j=start; j<start+len; j++) {
				const t = f[j]
				f[j] = (t + f[j + len]) % Q
				f[j + len] = (zeta * (Q + f[j + len] - t)) % Q
			}
		}
	}

	for (let i=0; i<256; i++) {
		f[i] = (f[i] * 3303) % Q  // 3303 == powQ(128, Q-2)
	}

	return f;
}



const mynumbers = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255];
const nttres = ntt(mynumbers);
console.log(JSON.stringify(nttres));
console.log(JSON.stringify(nttInv(nttres)));
const mybytes = byteEncode(12, mynumbers);
console.log(mybytes);
console.log(byteDecode(12, mybytes))