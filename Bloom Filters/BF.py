from sklearn.utils import murmurhash3_32
from math import log10, ceil, floor
from bitarray import bitarray
import random
import csv
import pandas as pd
import sys
import string
import matplotlib.pyplot as plt

seeds = [i for i in range(500)]

def hash_func(m, hash_tbl_size, seed = 0):
    return murmurhash3_32(m, seed, positive=True) % hash_tbl_size

class BloomFilter():
    def __init__(self, n, fp_rate) -> None:
        '''
        fp_rate = 0.618**(R/N)
        R: size of bitmap
        N: expected number of elements to hash
        k: number of hash functions to use
        '''
        R = ceil(log10(fp_rate)/log10(0.618)*n)
        tmp = 1
        while(tmp < R):
            tmp *= 2
        self.R = tmp if tmp - R <= R - (tmp//2) else tmp//2
        self.k = ceil(self.R/n*0.69314718)
        self.bit_arr = bitarray(self.R)
        self.bit_arr.setall(0)


    def insert(self, key):
        for seed in seeds[:self.k]:
            # hx = hash_func(key, self.R, seed)
            hx = murmurhash3_32(key, seed=seed, positive=True) % self.R
            self.bit_arr[hx] = 1

    def test(self, key):
        flag = True
        for seed in seeds[:self.k]:
            # hx = hash_func(key, self.R, seed)
            hx = murmurhash3_32(key, seed=seed, positive=True) % self.R
            if self.bit_arr[hx] == 0:
                flag = False
                break
        return flag


# generate datasets
random.seed(42)
membership_set = random.sample(range(10000, 100000), 10000)
test_set = random.sample(membership_set, 1000) + random.sample(set(range(10000, 100000)) - set(membership_set), 1000)

fp_rates = [0.01, 0.001, 0.0001]
for fp_rate in fp_rates:
    print("==================================")
    print("Theoretical fp_rate: ", fp_rate)
    bf = BloomFilter(len(membership_set), fp_rate)
    for elem in membership_set:
        bf.insert(elem)

    fp_count = 0
    for elem in test_set:
        test_res = bf.test(elem)
        if test_res and elem not in membership_set:
            fp_count += 1
    
    print("Real fp_rate:        ", round(fp_count/len(test_set), 5))

class BloomFilterURL():
    def __init__(self, n, R) -> None:
        '''
        k = 0.7*R/N
        R: size of bitmap
        N: expected number of elements to hash, N=377871 in this dataset
        k: number of hash functions to use
        '''
        self.R = R
        self.k = round(0.7*self.R/n)
        self.bit_arr = bitarray(self.R)
        self.bit_arr.setall(0)


    def insert(self, key):
        for seed in seeds[:self.k]:
            # hx = hash_func(key, self.R, seed)
            hx = murmurhash3_32(key, seed=seed, positive=True) % self.R
            self.bit_arr[hx] = 1

    def test(self, key):
        flag = True
        for seed in seeds[:self.k]:
            # hx = hash_func(key, self.R, seed)
            hx = murmurhash3_32(key, seed=seed, positive=True) % self.R
            if self.bit_arr[hx] == 0:
                flag = False
                break
        return flag

data = pd.read_csv("user-ct-test-collection-01.txt", sep="\t")
urllist = data.ClickURL.dropna().unique()


bf_url_list = [BloomFilterURL(len(urllist), 2**i) for i in range(17, 25)]

hash_tbl = dict()

count = 0
for url in urllist:
    for bf in bf_url_list:
        bf.insert(url)
    hash_tbl[url] = 1
    
    count += 1
    if count % 10000 == 0:
        print("progress: ", count)

bf_size = [sys.getsizeof(bf.bit_arr) for bf in bf_url_list]

hash_tbl_size = sys.getsizeof(hash_tbl)
print("size_hashtbl: " , sys.getsizeof(hash_tbl))
print("size_bf: ", bf_size)

test_set = random.sample(range(len(urllist)), 1000)

# generating random letters
letters = string.ascii_letters + string.digits
random_url_set = []

for i in range(1000):
    random_url_set.append(''.join(random.choice(letters) for i in range(10)))

bf_fp_count = [0]*len(bf_url_list)
for random_url in random_url_set:
    for i, bf in enumerate(bf_url_list):
        res = bf.test(random_url)
        if res and random_url in random_url_set:
            bf_fp_count[i] += 1
bf_fp_rate = [round(x/2000, 4) for x in bf_fp_count]
plt.xlabel('Memory comsumption [bytes]')
plt.ylabel('False Positive Rate')
plt.plot(bf_size, bf_fp_rate)
plt.show()