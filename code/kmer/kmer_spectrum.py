from pybloomfilter import BloomFilter
# from pybloom_live import BloomFilter
import pysam
import psutil
import os
import numpy as np

from matplotlib import pyplot as plt

import datetime


def get_kmers(read, ksize):
    kmers = []

    for i in range(len(read) - ksize + 1):
        kmers.append(read[i: i + ksize])
    return kmers


def plot_spectrum(stat_dict):
    frequency = []
    density = []
    for key in stat_dict.keys():
        frequency.append(key)
        density.append(stat_dict[key])

    fig = plt.figure(figsize=(17, 8))
    ax1 = fig.add_subplot(111)
    plt.xlim(0, 60, 2)
    new_ticks = range(0, 60, 5)
    plt.xticks(new_ticks)
    ax1.scatter(frequency, density, s=5)
    plt.xlabel('pos')
    plt.ylabel('qual')
    plt.show()


def get_hash(s, mod, p=131):
    # p:131或13331

    # 规定 idx(x) = x−′a′+1
    def idx(x):
        return ord(x) - ord('A') + 1

    hash_num = np.zeros((len(s)))
    hash_num[0] = idx(s[0])

    # hash[i]=(hash[i−1])∗p+idx(s[i]) % mod
    for i in range(1, len(s)):
        hash_num[i] = hash_num[i - 1] * p + idx(s[i])
    return int((hash_num[len(s) - 1]) % mod)


def get_usage():
    info = psutil.virtual_memory()
    print('内存使用：', psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024, "MB")
    print('总内存：', info.total / 1024 / 1024, 'MB')
    print('内存占比：', info.percent)
    print('打印本机cpu占用率： ' + str(psutil.cpu_percent(0)) + '%')

def get_spectrum(input_file, size=31):
    bams = pysam.AlignmentFile(input_file, 'rb')

    bloom_filter = BloomFilter(capacity=9999999999, error_rate=0.1)
    # print(bloom_filter.bitarray)
    # print(bloom_filter.num_bits)

    # 统计每一个kmer出现的次数，即多样性
    hash_dict = {}
    cnt = 0
    for r in bams:
        cnt += 1
        print(cnt)
        # if cnt % 100000 == 0:
        #     print(cnt)
        if cnt == 100000:
            break
        read = r.query_sequence
        kmers = get_kmers(read, size)
        # 将kmer加入到bloom中
        for kmer in kmers:
            # a= int(hash(kmer) % 3)
            # get_hash(kmer, 3)
            is_in = kmer in bloom_filter
            if is_in is True:
                # 排除假阳性
                if kmer in hash_dict:
                    hash_dict[kmer] += 1
                else:
                    hash_dict[kmer] = 1
            else:
                bloom_filter.add(kmer)
                hash_dict[kmer] = 1

    unique_kmer = []
    stat_dict = {}
    # 统计多样性为相同值的所有kmer个数,并删除只有一个kmer
    for key in hash_dict.keys():
        multiplicity = hash_dict[key]
        if multiplicity == 1:
            unique_kmer.append(key)
        elif multiplicity not in stat_dict.keys():
            stat_dict[multiplicity] = 1
        else:
            stat_dict[multiplicity] += multiplicity

    for i in range(len(unique_kmer)):
        hash_dict.pop(unique_kmer[i])

    get_usage()


    return stat_dict



    # print(bloom_filter.bitarray)


if __name__ == '__main__':

    start_time = datetime.datetime.now()

    stat_dict = get_spectrum('../../bams/sampleChr20_sorted.bam')
    print(stat_dict)
    # plot_spectrum(stat_dict)


    end_time = datetime.datetime.now()
    print((end_time - start_time).seconds)
