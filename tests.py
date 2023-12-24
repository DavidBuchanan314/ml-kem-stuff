from mlkem import *
from tqdm import tqdm

def get_tests():
	values = {}
	for line in open("KAT/MLKEM/kat_MLKEM_768.rsp"):
		k, v = line.strip().split(" = ")
		if k == "count":
			if values:
				yield values
				values = {}
			values[k] = int(v)
			continue
		values[k] = bytes.fromhex(v)
	yield values

for testcase in tqdm(get_tests()):
	#print(testcase["count"])

	ek, dk = mlkem_keygen(seed1=testcase["z"], seed2=testcase["d"])
	assert(ek == testcase["pk"])
	assert(dk == testcase["sk"])

	k, c = mlkem_encaps(ek, seed=testcase["msg"])
	assert(k == testcase["ss"])
	assert(c == testcase["ct"])

	k2 = mlkem_decaps(c, dk)
	assert(k2 == k)

	k3 = mlkem_decaps(testcase["ct_n"], dk)
	assert(k3 == testcase["ss_n"])