import random
import numpy as np
import math
import datetime
import sys
import time
import multiprocessing
import scipy.stats
import getopt
def SaveAsList(path):

    f = open(path, 'r')

    cnt = 0
    sequences = ''

    # 遍历fa文件，存储到字典中
    for line in f:
        if cnt < 100000:
            if line.startswith(">"):
                pass
            else:
                line = line.replace('N', '')
                sequences = sequences + line.rstrip("\n")
                cnt = cnt + 1
    # 字典的本地存储
    dict = open('../genome/chr20_list.txt', 'w')
    dict.write(str(sequences))
    return (sequences)

def SaveAsList2(path):

    f = open(path, 'r')

    cnt = 0
    sequences = ''

    # 遍历fa文件，存储到字典中
    for line in f:

        if line.startswith(">"):
            pass
        elif 'N' in line:
            pass
        else:
            cnt += 1
            print(cnt)
            sequences = sequences + line.rstrip("\n")
    # 字典的本地存储
    dict = open('../genome/chr20.txt', 'w')
    dict.write(str(sequences))
    return (sequences)

def SaveAsDict(path):

    f = open(path, 'r')

    cnt = 0
    sequences = {}
    # 遍历fa文件，存储到字典中
    for line in f:
        if cnt < 100000:
            if line.startswith(">"):
                # name = line.rstrip("\n")
                name = 'chr20'
                sequences[name] = ""
            else:
                line = line.replace('N', '')
                sequences[name] = sequences[name] + line.rstrip("\n")
                cnt = cnt + 1
    # 字典的本地存储
    dict = open('../genome/chr20_dict.txt', 'w')
    dict.write(str(sequences))
    return (sequences['chr20'])

def FetchFastqc(path1, path2=''):

    f = open(path1, 'r')
    outf = open('../genome/chr20_real_2.txt', 'w')
    cnt = 0
    # 遍历fa文件，存储到字典中
    for line in f:
        line = line.strip('\n')
        cnt += 1
        if cnt % 4 == 0:
            cnt = 0
            continue
        if cnt % 2 == 0:
            outf.write(line)

            # break

    outf.close()
    f.close()


def GetRead(fileName, seLen, readLen):
    # 随机选取长度为READ_LEN的片段
    index = random.randint(0, int(64428000 / 2))
    file = open(fileName, 'r')

    file.seek(index)
    read1 = file.read(seLen)

    # 方便第一个补成N
    file.seek(index + readLen - seLen - 1)
    read2 = file.read(seLen + 1)

    file.close()
    return read1, read2


def InverseBp(bp):
    if bp == 'A':
        return 'T'
    elif bp == 'T':
        return 'A'
    elif bp == 'C':
        return 'G'
    elif bp == 'G':
        return 'C'

def CountBps(l):
    bpDict = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0}

    for key in bpDict.keys():
        bpDict[key] = np.shape(np.where(l == key))[1]
    # bpDict['A'] = np.shape(np.where(l == 'A'))[1]

    # for i in range(len(l)):
    #     bpDict[l[i]] += 1
    return bpDict

# 计算质量
def calQual(errP, phredVal=33):
    if errP == 0:
        return chr(74)
    mostQual = int(-10 * math.log(errP, 10) + phredVal)
    if mostQual >= (74 + 70) / 2:
        mostQual = 74
    elif mostQual >= (70 + 65) / 2:
        mostQual = 70
    elif mostQual >= (65 + 60) / 2:
        mostQual = 65
    elif mostQual >= (60 + 55) / 2:
        mostQual = 60
    elif mostQual >= (55 + 45) / 2:
        mostQual = 55
    elif mostQual >= (45 + 41) / 2:
        mostQual = 45
    elif mostQual >= 41:
        mostQual = 41
    else:
        mostQual = mostQual

    # Hist({74: 43, 70: 29, 45: 28, 55: 16, 60: 14, 65: 13, 41: 7})
    # # 0.000079 0.00019        .......   0.15
    return chr(mostQual)


# 酶的活性曲线模型
# def GetEnzymeCurve(seLen):
#     xs = np.linspace(0.7, 3.0, seLen)
#     roll = 1 - np.exp(-10 * xs)
#     return np.array(roll)
# def GetEnzymeCurve(seLen):
#     xs = np.linspace(-0.1, 0.1, 6)
#     roll = np.power(xs, 3) + 0.999
#     roll = np.pad(roll, ((0, seLen)), 'constant', constant_values=(1))
#     return np.array(roll)
# def GetEnzymeCurve(seLen):
#     roll = np.array([0.9987, 0.9989, 0.99905, 0.9998, 0.9999, 0.9999])
#     roll = np.pad(roll, ((0, seLen)), 'constant', constant_values=(1))
#     return roll
def GetEnzymeCurve(seLen):
    roll = np.array([0.998, 0.999, 0.9992, 0.9993, 0.9994, 0.9995])
    roll = np.pad(roll, ((0, seLen)), 'constant', constant_values=(1))
    return roll
# 计算遗留的荧光总数目
def CalSumWrongLight(dict):
    sum = 0
    for key, value in dict.items():
        sum += value
    return sum

# def CountBps()
def SequenceRead1(read, seLen, factor=1):
    CLUSTER = 1000

    # read = list(GetRead('../genome/hs_ref_GRCh38.p12_chr20.fa', 150))
    read.append('N')
    read = np.array(read)

    curInserIndex = np.zeros(CLUSTER, dtype=np.int16) - 1
    lastInserIndex = np.zeros(CLUSTER, dtype=np.int16) - 1

    wrongHist = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0}

    resultBps = []
    resultQual = []
    resultErr = []
    # 生成酶的活性模型曲线
    enzymeCurve = GetEnzymeCurve(seLen)
    # print(enzymeCurve)
    # 总过需要length轮就可以结束, i也是现在应该插入的位置
    for i in range(seLen):
        # print('------------', i)
        curInserIndex = lastInserIndex.copy()

        # # 1. 酶活性的影响,活性影响是否补上
        # # 酶的活性， 即酶影响错误的概率
        enzyme = round(1 - enzymeCurve[i], 3)
        enzymeRoll = np.random.randint(1, 1001, CLUSTER) / 1001
        lowDNA = np.shape(np.where(enzymeRoll < enzyme))[1]


        # 2.1 保护因子(叠氮基团)没有去掉概率(导致不能补上当前位置的bp)
        roll = np.random.randint(1, int(10000001 / factor), CLUSTER)
        # print(roll)
        rightDNA = np.where(roll > 5)
        curInserIndex[rightDNA] += 1

        # 修正index，防止越界
        curInserIndex[curInserIndex >= seLen] = seLen

        # # 3. 计算错误率
        # 所有当前轮补上碱基的位置
        inserIndex = np.where((curInserIndex - lastInserIndex) >= 1)

        # # 当前这轮补上的碱基
        curInserBp = read[curInserIndex[inserIndex]]

        # #找到数量最多的碱基
        counterBps = CountBps(curInserBp)
        # print(counterBps)
        mostNum = -1
        for key, value in counterBps.items():
            if mostNum <= value:
                mostNum = value
                mostBps = key
        resultBps.append(mostBps)

        # 计算错误率
        errP = 1 - ((mostNum + wrongHist[mostBps] - lowDNA) / (np.shape(inserIndex)[1] + CalSumWrongLight(wrongHist)))
        # errP = 1 - ((mostNum + wrongHist[mostBps] - lowDNA) / (np.shape(inserIndex)[1] + CalSumWrongLight(wrongHist) - counterBps['N']))
        resultErr.append(errP)
        # 修改factor
        factor += errP * i * 40
        mostQual = calQual(errP)

        # print(mostBps, errP, mostQual)
        resultQual.append(mostQual)

        # # 4. 荧光基团的影响
        # 荧光基团没有切除，导致下一轮的荧光数目增大
        # roll = np.random.randint(1, int(10000001 / factor), CLUSTER)
        # wrongDNA = np.where(roll <= 5)
        # wrongLight = read[curInserIndex[wrongDNA]]
        # # 统计错误的碱基和数目（返回类型：({'G': 2, 'T': 2})）
        # wrongHist = CountBps(wrongLight)

        # 5 保护因子(叠氮基团)过早去掉(导致多补上一个)
        roll = np.random.randint(1, int(10000001 / factor), CLUSTER)
        rightDNA = np.where(roll <= 5)
        curInserIndex[rightDNA] += 1

        lastInserIndex = curInserIndex
        # break
    # print(read)
    # print(resultBps)
    # print(resultErr)
    # print(resultQual)
    # fig = plt.figure(figsize=(17, 8))
    # ax1 = fig.add_subplot(121)
    # ax2 = fig.add_subplot(122)
    # x = list(range(1, 151))
    # # ax1.scatter(x, resultQual, s=5)
    # ax2.scatter(x, resultErr, s=5)
    # plt.xlabel('pos')
    # plt.ylabel('qual')
    # plt.show()
    return resultBps, resultQual, factor


def SequenceRead2(read, seLen, factor=1):
    CLUSTER = 1000

    if factor > 1:
        factor = 5
    elif factor > 2:
        factor = 10

    read[0] = 'N'
    read = np.array(read)

    curInserIndex = np.zeros(CLUSTER, dtype=np.int16) - 1
    lastInserIndex = np.ones(CLUSTER, dtype=np.int16) * (seLen + 1)

    wrongHist = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0}
    resultBps = []
    resultQual = []
    resultErr = []

    # 生成酶的活性模型曲线
    enzymeCurve = GetEnzymeCurve(seLen)
    # print(enzymeCurve)
    # 总过需要length轮就可以结束, i也是现在应该插入的位置
    for i in range(seLen):
        # print('------------', i)
        curInserIndex = lastInserIndex.copy()

        # # 1. 酶活性的影响,活性影响是否补上
        # # 酶的活性， 即酶影响错误的概率
        enzyme = round(1 - enzymeCurve[i], 3)
        enzymeRoll = np.random.randint(1, 1001, CLUSTER) / 1001
        lowDNA = np.shape(np.where(enzymeRoll < enzyme))[1]

        # 2.1 保护因子(叠氮基团)没有去掉概率(导致不能补上当前位置的bp)
        roll = np.random.randint(1, int(10000001 / factor), CLUSTER)
        # print(roll)
        rightDNA = np.where(roll > 5)
        curInserIndex[rightDNA] -= 1

        # 修正index，防止越界
        curInserIndex[curInserIndex <= 0] = 0

        # # 3. 计算错误率
        # 所有当前轮补上碱基的位置
        inserIndex = np.where((lastInserIndex - curInserIndex) >= 1)

        # # 当前这轮补上的碱基
        curInserBp = read[curInserIndex[inserIndex]]
        # print(curInserIndex[inserIndex].values)

        # #找到数量最多的碱基
        counterBps = CountBps(curInserBp)
        # print(counterBps)
        mostNum = -1
        for key, value in counterBps.items():
            if mostNum <= value:
                mostNum = value
                mostBps = key

        resultBps.append(InverseBp(mostBps))

        # 计算错误率
        errP = 1 - ((mostNum + wrongHist[mostBps] - lowDNA) / (np.shape(inserIndex)[1] + CalSumWrongLight(wrongHist)))
        resultErr.append(errP)
        factor += errP * i * 40

        mostQual = calQual(errP)

        # print(mostBps, errP, mostQual)
        resultQual.append(mostQual)

        # # 4. 荧光基团的影响
        # 荧光基团没有切除，导致下一轮的荧光数目增大
        # roll = np.random.randint(1, int(10000001 / factor), CLUSTER)
        # wrongDNA = np.where(roll <= 2)
        # wrongLight = read[curInserIndex[wrongDNA]]
        # # 统计错误的碱基和数目（返回类型：({'G': 2, 'T': 2})）
        # wrongHist = CountBps(wrongLight)

        # 5 保护因子(叠氮基团)过早去掉(导致多补上一个)
        roll = np.random.randint(1, int(10000001 / factor), CLUSTER)
        rightDNA = np.where(roll <= 2)
        curInserIndex[rightDNA] -= 1


        lastInserIndex = curInserIndex
        # break
    # print(read)
    # print(resultBps)
    # print(resultErr)
    # print(resultQual)
    # fig = plt.figure(figsize=(17, 8))
    # ax1 = fig.add_subplot(121)
    # ax2 = fig.add_subplot(122)
    # x = list(range(1, 151))
    # # ax1.scatter(x, resultQual, s=5)
    # ax2.scatter(x, resultErr, s=5)
    # plt.xlabel('pos')
    # plt.ylabel('qual')
    # plt.show()
    return resultBps, resultQual, factor

def consumeT(contentPool, seLen, outFile1, outFile2):

    fout1 = open(outFile1, mode='w')
    fout2 = open(outFile2, mode='w')
    cnt = 1
    flag = 0
    while True:

        if cnt == 10:
            fout1.close()
            fout2.close()
            break

        # mutexLock.acquire()
        while len(contentPool) is not 0:
            flag = 1
            # print('Writing', contentPool[0][0])
            if contentPool[0][0] % 1000 == 0:
                print("Writing: ", contentPool[0][0])

            fout1.write('@Chr20_%s_:%s   LENGTH=%s\n' % (str(contentPool[0][0]), str(1), str(seLen)))
            fout1.write(''.join(contentPool[0][1]))
            fout1.write('\n')
            fout1.write('+Chr20_%s_:%s   LENGTH=%s\n' % (str(contentPool[0][0]), str(1), str(seLen)))
            fout1.write(''.join(contentPool[0][2]))
            fout1.write('\n')

            fout2.write('@Chr20_%s_:%s   LENGTH=%s\n' % (str(contentPool[0][0]), str(1), str(seLen)))
            fout2.write(''.join(contentPool[0][3]))
            fout2.write('\n')
            fout2.write('+Chr20_%s_:%s   LENGTH=%s\n' % (str(contentPool[0][0]), str(1), str(seLen)))
            fout2.write(''.join(contentPool[0][4]))
            fout2.write('\n')

            del contentPool[0]
        # mutexLock.release()

        if flag == 0:
            time.sleep(0.1)
            cnt = cnt + 1
        else:
            flag = 0
            cnt = 1

def produceT(i, contentPool, seLen, inFile):
    if i % 1000 == 0:
        print("Sequencing: ", i)

    roll = np.random.randint(1, 10) / 10
    readLen = int(scipy.stats.norm.ppf(roll, 500, 50))
    read1, read2 = GetRead(inFile, seLen, readLen)

    resultBps1, resultQual1, factor1 = SequenceRead1(list(read1), seLen)
    resultBps2, resultQual2, factor2 = SequenceRead2(list(read2), seLen, factor1)
    curRead = [i, resultBps1, resultQual1, resultBps2, resultQual2]
    # print('Sequencing: ', curRead[0])

    contentPool.append(curRead)


def StartMethod(inFile='../genome/chr20.txt', readsNum=5000, seLen=150, outFile='../genome/readTT'):
    startTime = datetime.datetime.now()

    outFile1 = outFile + '1.fastq'
    outFile2 = outFile + '2.fastq'
    processPool = multiprocessing.Pool(8)
    contentPool = multiprocessing.Manager().list()

    processPool.apply_async(consumeT, args=(contentPool, seLen, outFile1, outFile2))

    for i in range(readsNum):
        processPool.apply_async(produceT, args=(i, contentPool, seLen, inFile))

    processPool.close()
    processPool.join()

    print('All subprocesses done.')
    endTime = datetime.datetime.now()
    print('time: ', (endTime - startTime).seconds)



if __name__ == '__main__':
    if len(sys.argv) == 1:
        StartMethod()
    else:
        inFile = '../genome/chr20_list.txt'
        outFile = '../genome/read'
        readsNum = 5000
        seLen = 150

        try:
            opts, args = getopt.getopt(sys.argv[1:], "hi:o:n:l:")
        except getopt.GetoptError:
            print('test.py -i <inputfile> -o <outputfile> -l <sequence length> -n <number of reads>')
            sys.exit(2)
        for opt, arg in opts:
            if opt == '-h':
                print('pairM_NP_MP.py -i <inputfile> -o <outputfile> -l <sequence length> -n <number of reads>')
                print("-i   The whole path of input file")
                print("-o   The path and name of output file")
                print("-l   The length of reads")
                print("-n   Totle numbers of reads")
                sys.exit()
            elif opt in ("-i"):
                inFile = arg
            elif opt in ("-o"):
                outFile = arg
            elif opt in ("-n"):
                readsNum = int(arg)
            elif opt in ("-l"):
                seLen = int(arg)
        StartMethod(inFile, readsNum, seLen, outFile)


