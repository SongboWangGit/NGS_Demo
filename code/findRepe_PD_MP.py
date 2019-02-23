import pysam
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import numpy as np
import os
import sys
import getopt
import multiprocessing
import shutil
import datetime


def FindRepe(read, qual, STARTCUT, ENDCUT, READ, repeNum):
    REPELEN = 1
    words = ['A', 'T', 'C', 'G']
    strLen = len(read)
    repeDict = {}

    # 重复的个数
    for i in range(1, REPELEN + 1):
        # 找到words的全排列
        allAlign = itertools.product(words, repeat=i)

        for each in allAlign:

            each = ''.join(list(each))

            startIndex = 0
            index = read.find(each, startIndex)
            indexLists = []

            # 找到所有出现的位置，返回位置的列标
            while startIndex <= strLen - i and index != -1:
                indexLists.append(index)
                startIndex = index + i
                index = read.find(each, startIndex)

            # 存放连续index的列表
            continueLists = []

            # 检查列表，是否有连续的
            if len(indexLists) > 1:
                # 存储当前连续index
                continueList = []
                continueList.append(indexLists[0] + STARTCUT)

                # 逐个检查当前项和前一项
                for m in range(1, len(indexLists)):
                    if indexLists[m] == indexLists[m - 1] + i:
                        continueList.append(indexLists[m] + STARTCUT)
                    # 不再连续时
                    else:
                        if len(continueList) > 1:
                            continueLists.append(continueList)
                        continueList = []
                        continueList.append(indexLists[m] + STARTCUT)

            if len(continueLists) >= 1:
                # print(each, continueLists)
                repeDict[each] = continueLists

    for key in repeDict.keys():
        continueLists = repeDict[key]
        for indexs in continueLists:
            if len(indexs) == int(repeNum):
                WriteToFile(key, indexs, qual, READ, repeNum)

    return repeDict


def WriteToFile(key, indexs, qual, READ, repeNum):
    outFile = open("../repe/repe_%s/repe%s_txt/%s.txt" % (repeNum, READ, key), 'a')
    startIndex = indexs[0]
    endIndex = indexs[len(indexs) - 1] + len(key) - 1
    # print(key, indexs[0], indexs[len(indexs) - 1] + len(key) - 1,  (qual[startIndex:endIndex]))
    if endIndex < 144 and startIndex >= 5:
        outFile.write('%s, %s, ' % (startIndex, endIndex))
        outFile.write(str(qual[startIndex - 5:endIndex + 1 + 5])[1:-1])
        outFile.write('\n')


def saveRepePlot(qualSeries, repeDict, cnt, READ):
    # # 绘制初始质量图
    fig = plt.figure(figsize=(13, 6.5))
    ax1 = fig.add_subplot(111)
    ax1.scatter(qualSeries.index, qualSeries.values, s=5)

    allRepeList = []

    # 找到重复区域的质量的index
    for key in repeDict.keys():
        continueLists = repeDict[key]

        for lists in continueLists:
            # 得到所有重复序列的index
            cutSeries = qualSeries[lists[0]:lists[len(lists) - 1] + len(key)]
            allRepeList.extend(cutSeries.index)

    allRepeList.sort()

    qualRepeSeries = qualSeries[allRepeList]

    ax1.scatter(qualRepeSeries.index, qualRepeSeries.values, color='red', s=5)
    fig.savefig('../repe/repe_img%s/test%s.jpg' % (READ, cnt))
    plt.close()


def CalAveQual(READ, repeNum):

    rootdir = '../repe/repe_%s/repe%s_txt' % (repeNum, READ)
    files = os.listdir(rootdir)  # 列出文件夹下所有的目录与文件

    for file in files:
        # 计算文件全部路径
        path = os.path.join(rootdir, file)
        # 是文件而不是文件夹时
        if os.path.isfile(path):

            names = ['s', 'e']
            for i in range(int(repeNum) + 5 + 5):
                names.append('q')

            qualInfoDF = pd.DataFrame(pd.read_table(path, sep=',',
                                           names=names))

            # 计算全部的平均
            aveQual = qualInfoDF.iloc[:, 2:].mean().tolist()
            SaveAveImg(file, aveQual, READ, repeNum)

            # 根据起始位置分段求平均
            qualInfoDF['s'] = (qualInfoDF['s'] / 10).astype('int')
            for i in range(15):
                aveQual = qualInfoDF[qualInfoDF['s'] == i].iloc[:, 2:].mean().tolist()
                SaveAveImg(file, aveQual, READ, repeNum, i)



def SaveAveImg(file, aveQual, READ, repeNum, i=-1):
    fig = plt.figure(figsize=(13, 6.5))
    ax1 = fig.add_subplot(111)
    ax1.plot(aveQual, marker='*')
    # 绘制分割线
    plt.vlines(4, min(aveQual), max(aveQual), colors="r", linestyles="dashed")
    plt.vlines(len(aveQual) - 1 - 5, min(aveQual), max(aveQual), colors="r", linestyles="dashed")
    plt.xlabel('Repeated bps')
    plt.ylabel('Quality')
    if i == -1:
        fig.savefig('../repe/repe_%s/repe%s_ave_img/%sall.jpg' % (repeNum, READ, file.split('.')[0]))
    else:
        fig.savefig('../repe/repe_%s/repe%s_ave_img/%s%s.jpg' % (repeNum, READ, file.split('.')[0], i))

    plt.cla()
    plt.close('all')


def JudgeRead1Or2(flag):
    binNum = bin(flag).replace('0b', '')
    binNumStr = ''.join(binNum)
    length = len(binNumStr)
    # print(binNum)
    # print(length)
    # print(binNumStr[-length])
    if length == 7:
        return 1

    if binNumStr[-7] == '1':
        return 1
    elif binNumStr[-8] == '1':
        return 2


def FindRepeStart(READ, inputFile, repeNum):
    startTime = datetime.datetime.now()
    frontNum = 5
    STARTCUT = 0
    ENDCUT = 0
    bf = pysam.AlignmentFile(inputFile, 'rb')
    cnt = -1

    for r in bf:
        if cnt % 10000 == 0:
            print('Searching Read%s: %s' % (READ, cnt))

        cnt = cnt + 1
        if cnt == 100000:
            break
        rname = r.qname
        rflag = r.flag

        # 去掉read2
        if JudgeRead1Or2(rflag) == READ:
            # 获得bp串
            read = r.query_sequence

            # # 切除开头和结尾的部分碱基
            # readCut = read[STARTCUT:len(read) - ENDCUT]
            # 获得数值质量
            qual = list(r.query_qualities)
            qualSeries = pd.Series(qual)

            # 找到重复区域的index
            repeDict = FindRepe(read, qual, STARTCUT, ENDCUT, READ, repeNum)

            # 保存图像
            # saveRepePlot(qualSeries, repeDict, cnt)

    bf.close()
    endTime = datetime.datetime.now()
    print("Time: ", (endTime - startTime).seconds)
    CalAveQual(READ, repeNum)


if __name__ == "__main__":

    # 使用多进程
    processPool = multiprocessing.Pool(3)

    if not os.path.exists('../repe'):
        os.makedirs('../repe')

    repeNum = 8
    inputFile = '../bams/sampleChr20_sorted.bam'

    if len(sys.argv) == 1:
        repeNum = 10
        inputFile = '../bams/sampleChr20_sorted.bam'
    else:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "hi:n:")
        except getopt.GetoptError:
            print('findRepe.py -i <inputfile>')
            sys.exit(2)
        for opt, arg in opts:
            if opt == '-h':
                print('findRepe.py -i <inputfile>')
                sys.exit()
            elif opt in ("-i"):
                inputFile = arg
            elif opt in ("-n"):
                repeNum = arg

    # 删除并创建结果文件夹
    if os.path.exists('../repe/repe_%s' % repeNum):
        shutil.rmtree('../repe/repe_%s' % repeNum)
    os.makedirs('../repe/repe_%s' % repeNum)

    os.makedirs('../repe/repe_%s/repe1_ave_img' % repeNum)
    os.makedirs('../repe/repe_%s/repe2_ave_img' % repeNum)
    os.makedirs('../repe/repe_%s/repe1_txt' % repeNum)
    os.makedirs('../repe/repe_%s/repe2_txt' % repeNum)

    for READ in [1, 2]:
        processPool.apply_async(FindRepeStart, args=(READ, inputFile, repeNum))

    processPool.close()
    processPool.join()

