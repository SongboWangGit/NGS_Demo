from pybloomfilter import BloomFilter
# from pybloom_live import BloomFilter
import pysam
import psutil
import os
from matplotlib import pyplot as plt
import multiprocessing
import datetime
import numpy as np
import time
def get_kmers(read, ksize):
    kmers = []

    for i in range(len(read) - ksize + 1):
        kmers.append(read[i: i + ksize])
    return kmers

def plot_spectrum(frequency, density):
    fig = plt.figure(figsize=(17, 8))
    ax1 = fig.add_subplot(111)
    plt.xlim(0, 60, 2)
    new_ticks = range(0, 60, 5)
    plt.xticks(new_ticks)
    ax1.scatter(frequency, density, s=5)
    plt.xlabel('pos')
    plt.ylabel('qual')
    plt.show()


def get_spectrum(input_file, size=31):
    bams = pysam.AlignmentFile(input_file, 'rb')

    bloom_filter = BloomFilter(capacity=999999999, error_rate=0.1)
    # print(bloom_filter.bitarray)
    # print(bloom_filter.num_bits)

    # 统计每一个kmer出现的次数，即多样性
    hash_dict = {}
    cnt = 0
    for r in bams:
        cnt += 1
        print(cnt)
        if cnt == 200000:
            break
        read = r.query_sequence
        kmers = get_kmers(read, size)
        # print(kmers)
        for kmer in kmers:

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

    # 删除只有一个kmer
    unique_kmer = []
    for key in hash_dict.keys():
        if hash_dict[key] == 1:
            unique_kmer.append(key)

    for i in range(len(unique_kmer)):
        hash_dict.pop(unique_kmer[i])
    # print(hash_dict)

    # 统计多样性为相同值的所有kmer个数
    stat_dict = {}

    for key in hash_dict.keys():

        multiplicity = hash_dict[key]
        if multiplicity not in stat_dict.keys():
            stat_dict[multiplicity] = 1
        else:
            stat_dict[multiplicity] += multiplicity

    frequency = []
    density = []
    for key in stat_dict.keys():
        frequency.append(key)
        density.append(stat_dict[key])




    return stat_dict, frequency, density





    # print(bloom_filter.bitarray)


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


def master_node(in_file, kmer_pool, ksize, slave_num):
    try:
        bams = pysam.AlignmentFile(in_file, 'rb')
        cnt = 0
        for r in bams:
            cnt += 1
            print(cnt)
            if cnt % 10 == 0:
                for i in range(slave_num):
                    while len(kmer_pool[i]) != 0:
                        print("sleep")
                        time.sleep(0)
            if cnt == 1000:
                break
            read = r.query_sequence
            kmers = get_kmers(read, ksize)

            for kmer in kmers:
                # print(kmer)
                # print(kmer_pool)
                hash_num = get_hash(kmer, slave_num)

                tmp_list = kmer_pool[hash_num]
                tmp_list.append(kmer)
                kmer_pool[hash_num] = tmp_list

                # print(type(kmer_pool[hash_num]))
            # break
        print("Master End")
    except Exception as ex:
        print(ex)


def merge_dict(x, y):
    for k, v in y.items():
        if k in x.keys():
            x[k] += v
        else:
            x[k] = v
    return x


def slave_node(kmer_pool, slave_num, stat_pool):
    local_filter = BloomFilter(capacity=9999999999, error_rate=0.1)

    local_dict = {}
    try:
        print('Slave: ', slave_num,  ' start')
        time_to_stop = 0
        while True:
            self_kmres = kmer_pool[slave_num]
            while len(self_kmres) is not 0:
                time_to_stop = 0
                kmer = self_kmres[0]
                # print('Slave: ', slave_num, " del: ", kmer)
                del self_kmres[0]

                # 加入到Bloom filter中
                is_in = kmer in local_filter
                if is_in is True:
                    # 排除假阳性
                    if kmer in local_dict:
                        local_dict[kmer] += 1
                    else:
                        local_dict[kmer] = 1
                else:
                    local_filter.add(kmer)
                    local_dict[kmer] = 1

                kmer_pool[slave_num] = self_kmres

            time.sleep(0.001)
            time_to_stop += 1

            # 停止时间到
            if time_to_stop == 50:
                # 删除只有一个kmer
                unique_kmer = []
                local_stat = {}

                # 去除唯一的kmer,统计多样性为相同值的所有kmer个数
                for key in local_dict.keys():
                    multiplicity = local_dict[key]
                    if multiplicity == 1:
                        unique_kmer.append(key)
                    elif multiplicity not in local_stat.keys():
                        local_stat[multiplicity] = 1
                    else:
                        local_stat[multiplicity] += multiplicity

                for i in range(len(unique_kmer)):
                    local_dict.pop(unique_kmer[i])
                stat_pool = merge_dict(stat_pool, local_stat)
                break
        print('Slave: ', slave_num, ' end------------------')
    except Exception as ex:
        print('Slave: ', slave_num, ex)


def Start(in_file, slave_num=5, ksize=31):
    process_pool = multiprocessing.Pool(slave_num + 2)
    kmer_pool = multiprocessing.Manager().dict()
    stat_pool = multiprocessing.Manager().dict()
    for i in range(slave_num):
        kmer_pool[i] = []

    # master_node(in_file, kmer_pool, ksize, slave_num)
    # slave_node(kmer_pool, 0)

    process_pool.apply_async(master_node, args=(in_file, kmer_pool, ksize, slave_num))

    for i in range(slave_num):
        process_pool.apply_async(slave_node, args=(kmer_pool, i, stat_pool))
    process_pool.close()
    process_pool.join()
    #
    print(kmer_pool)
    print(stat_pool)

    info = psutil.virtual_memory()
    print('内存使用：', psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024, "MB")
    print('总内存：', info.total / 1024 / 1024, 'MB')
    # for i in range(SLAVES):
    #     process_pool.apply_async(produceT, args=(i, contentPool, seLen, inFile))

if __name__ == '__main__':
    start_time = datetime.datetime.now()
    in_file = '../../bams/sampleChr20_sorted.bam'
    Start(in_file)

    end_time = datetime.datetime.now()
    print((end_time - start_time).seconds)

    # start_time = datetime.datetime.now()
    #
    # stat_dict, frequency, density = get_spectrum('../../bams/sampleChr20_sorted.bam')
    # print(stat_dict)
    # plot_spectrum(frequency, density)
    #
    #
    # end_time = datetime.datetime.now()
    # print((end_time - start_time).seconds)
