import pysam
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import numpy as np
import os

def FindRepe(read, qual, STARTCUT, ENDCUT):
    REPELEN = 3
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
            if len(indexs) == 10:
                WriteToFile(key, indexs, qual)

    return repeDict


def WriteToFile(key, indexs, qual):
    outFile = open("../repe_txt2/%s.txt" % key, 'a')
    startIndex = indexs[0]
    endIndex = indexs[len(indexs) - 1] + len(key) - 1
    print(key, startIndex, endIndex, qual[startIndex:endIndex])
    # print(key, indexs[0], indexs[len(indexs) - 1] + len(key) - 1, qual,  (qual[startIndex:endIndex]))

    outFile.write('(%s,%s)--' % (startIndex, endIndex))
    outFile.write(str(qual[startIndex:endIndex + 1]))
    outFile.write('\n')


def savePlot(qualSeries, repeDict, cnt):

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
    fig.savefig('../repe_img/test%s.jpg' % cnt)
    plt.close()


def CalAveQual():
    rootdir = '../repe_txt'
    files = os.listdir(rootdir)  # 列出文件夹下所有的目录与文件
    for file in files:
        path = os.path.join(rootdir, file)
        if os.path.isfile(path):
            print(path)
            with open(path) as inputFile:

                line0 = inputFile.readline().strip().split('--')

                sumQual = np.array(list(map(int, line0[1][1:-1].split(', '))))

                cnt = 1

                for line in inputFile.readlines():
                    cnt += 1

                    lineSplit = line.strip().split('--')
                    index = lineSplit[0]

                    # str的list转化为int的list
                    qual = list(map(int, lineSplit[1][1:-1].strip().split(', ')))

                    # 逐个相加
                    sumQual = qual + sumQual

                # print(sumQual)

            aveQual = sumQual / cnt
            print(aveQual)

            fig = plt.figure(figsize=(13, 6.5))
            ax1 = fig.add_subplot(111)
            ax1.plot(aveQual)
            fig.savefig('../repe_ave_img/%s.jpg' % file)

def repeStart():
    frontNum = 5
    STARTCUT = 0
    ENDCUT = 0
    bf = pysam.AlignmentFile('../bams/sampleChr20_sorted.bam', 'rb')
    cnt = -1

    for r in bf:
        cnt = cnt + 1

        qname = r.qname

        # 获得bp串
        read = r.query_sequence

        # 切除开头和结尾的部分碱基
        readCut = read[STARTCUT:len(read) - ENDCUT]
        # 获得数值质量
        qual = list(r.query_qualities)
        qualSeries = pd.Series(qual)

        # 找到重复区域的index
        repeDict = FindRepe(read, qual, STARTCUT, ENDCUT)

        # 保存图像
        # savePlot(qualSeries, repeDict, cnt)

        # break
    bf.close()


if __name__ == "__main__":


    # repeStart()
    CalAveQual()


