import random
import numpy as np
import hist_demo
import math
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
import distribution_demo
import thinkplot
import datetime
import fcntl
import threading
import time
import multiprocessing

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

def GetRead(fileName, length):
    # 随机选取长度为READ_LEN的片段
    index = random.randint(0, 10000)
    file = open(fileName, 'r')
    file.seek(index)
    read = file.read(length)
    file.close()
    return read

    # return sequences[index : index + length]

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
    # 0.000079 0.00019        .......   0.15
    # if mostQual > 74:
    #     mostQual = 74
    # elif mostQual > 70:
    #     mostQual = 70
    # elif mostQual > 65:
    #     mostQual = 65
    # elif mostQual > 60:
    #     mostQual = 60
    # elif mostQual > 55:
    #     mostQual = 55
    # elif mostQual > 45:
    #     mostQual = 45
    # elif mostQual > 41:
    #     mostQual = 41
    # return mostQual
    return chr(mostQual)


# 酶的活性曲线模型
def getEnzymeCurve():
    expoModel = distribution_demo.ExpoCDF()
    xs, roll = expoModel.RenderExpoCdf(10, 0.05, 3.0, 150)
    roll = np.array(roll)

    return roll


# 计算遗留的荧光总数目
def calSumWrongLight(dict):
    sum = 0
    for key, value in dict.Items():
        sum += value
    return sum


def Sequence(read, length):
    CLUSTER = 1000

    # read = list(GetRead('../genome/hs_ref_GRCh38.p12_chr20.fa', 150))
    read.append('N')
    read = np.array(read)


    lastInserIndex = pd.Series(np.zeros(CLUSTER, dtype=np.int16) - 1)
    lastLightIndex = pd.Series(np.zeros(CLUSTER, dtype=np.int16))
    curInserIndex = []

    wrongHist = hist_demo.Hist()
    resultBps = []
    resultQual = []
    resultErr = []

    # 生成酶的活性模型曲线
    enzymeCurve = getEnzymeCurve()

    # 总过需要length轮就可以结束, i也是现在应该插入的位置
    for i in range(len(read) - 1):
        # print('------------', i)
        curInserIndex = lastInserIndex.copy()

        # # # 1. 酶活性的影响,活性影响是否补上
        # # 酶的活性， 即酶影响错误的概率
        # enzyme = round(1 - enzymeCurve[i], 3)
        # enzymeRoll = pd.Series(np.random.randint(1, 1001, CLUSTER))

        # # 2.1 保护因子(叠氮基团)没有去掉概率(导致不能补上当前位置的bp)
        roll = pd.Series(np.random.randint(1, 1000001, CLUSTER))
        # print(roll)
        rightDNA = roll[roll > 2].index
        # # 检查酶的活性
        # curEnzymeRoll = enzymeRoll[rightDNA]
        # rightDNA = curEnzymeRoll[curEnzymeRoll > enzyme].index
        curInserIndex[rightDNA] += 1

        # 修正index，防止越界
        curInserIndex[curInserIndex >= length] = length

        # # 3. 计算错误率
        # 所有当前轮补上碱基的位置
        inserIndex = (curInserIndex - lastInserIndex)
        inserIndex = inserIndex[inserIndex >= 1].index

        # # 当前这轮补上的碱基
        curInserBp = read[curInserIndex[inserIndex].values]
        # print(curInserIndex[inserIndex].values)

        # #找到数量最多的碱基
        counterBps = Counter(curInserBp)
        # print(counterBps)
        mostNum = -1
        for key, value in counterBps.items():
            if mostNum <= value:
                mostNum = value
                mostBps = key
        resultBps.append(mostBps)

        # 计算错误率
        # print(wrongHist[mostBps])
        errP = 1 - ((mostNum + wrongHist[mostBps]) / (inserIndex.size + calSumWrongLight(wrongHist) - counterBps['N']))
        resultErr.append(errP)
        mostQual = calQual(errP)
        # print(mostBps, errP, mostQual)
        resultQual.append(mostQual)

        # # 4. 荧光基团的影响
        # 荧光基团没有切除，导致下一轮的荧光数目增大
        roll = pd.Series(np.random.randint(1, 10000001, CLUSTER))
        wrongDNA = roll[roll <= 2].index
        wrongLight = read[curInserIndex[wrongDNA].values]
        # 统计错误的碱基和数目（返回类型：({'G': 2, 'T': 2})）
        wrongHist = hist_demo.Hist(wrongLight)

        # 5 保护因子(叠氮基团)过早去掉(导致多补上一个)
        roll = pd.Series(np.random.randint(1, 10000001, CLUSTER))
        rightDNA = roll[roll <= 2].index
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
    return resultBps, resultQual


# 将结果写入fq文件
def WriteFile(fout, resultBps, resultQual, name, num, length):
    fout.write('@' + name + ':' + str(num) + ' LENGTH=%s\n' % str(length))

    fout.write(''.join(resultBps))

    fout.write('\n')

    fout.write('+' + name + ':' + str(num) + ' LENGTH=%s\n' % str(length))
    fout.write(''.join(resultQual))
    fout.write('\n')


def Generate(readLength):
    startTime = datetime.datetime.now()

    out_fq = '../genome/out.fastq'
    fout = open(out_fq, mode='a')

    for i in range(1, 50):
        nowTime = (datetime.datetime.now() - startTime).seconds
        read = GetRead('../genome/chr20_list.txt',  readLength)
        print(i, nowTime)

        resultBps, resultQual = Sequence(list(read), readLength)

        WriteFile(fout, resultBps, resultQual, 'Chr20_%s ' % i, 1, readLength)

    fout.close()
    endTime = datetime.datetime.now()

    print('sum: ', (endTime - startTime).seconds)


def consumeT(contentPool):
    global READLENGTH
    out_fq = '../genome/out.fastq'
    fout = open(out_fq, mode='w')

    cnt = 1
    flag = 0
    while True:

        if cnt == 25:
            fout.close()
            break

        # mutexLock.acquire()
        while len(contentPool) is not 0:
            flag = 1
            print('Writing', contentPool[0][0])

            fout.write('@Chr20_%s_:%s   LENGTH=%s\n' % (str(contentPool[0][0]), str(1), str(READLENGTH)))
            fout.write(''.join(contentPool[0][1]))
            fout.write('\n')
            fout.write('+Chr20_%s_:%s   LENGTH=%s\n' % (str(contentPool[0][0]), str(1), str(READLENGTH)))
            fout.write(''.join(contentPool[0][2]))
            fout.write('\n')

            del contentPool[0]
        # mutexLock.release()

        if flag == 0:
            time.sleep(1)
            cnt = cnt + 1
        else:
            flag = 0
            cnt = 1



def produceT(i, contentPool):
    global READLENGTH

    read = GetRead('../genome/chr20_list.txt', READLENGTH)

    resultBps, resultQual = Sequence(list(read), READLENGTH)
    # print(resultQual)
    curRead = [i, resultBps, resultQual]
    print(curRead)

    contentPool.append(curRead)




if __name__ == '__main__':

    # Generate(150)
    startTime = datetime.datetime.now()

    READLENGTH= 150
    mutexLock = multiprocessing.Lock()
    processPool = multiprocessing.Pool(5)
    contentPool = multiprocessing.Manager().list()

    processPool.apply_async(consumeT, args=(contentPool, ))

    for i in range(5000):
        processPool.apply_async(produceT, args=(i, contentPool))

    processPool.close()
    processPool.join()

    print('All subprocesses done.')
    endTime = datetime.datetime.now()
    print('sum: ', (endTime - startTime).seconds)
