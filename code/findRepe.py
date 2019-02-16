import pysam
import pandas as pd
import matplotlib.pyplot as plt
import itertools


def FindRepe(queryStr):
    REPELEN = 3
    words = ['A', 'T', 'C', 'G']
    strLen = len(queryStr)
    repeDict = {}

    # 重复的个数
    for i in range(1, REPELEN + 1):
        # 找到words的全排列
        allAlign = itertools.product(words, repeat=i)

        for each in allAlign:

            each = ''.join(list(each))

            startIndex = 0
            index = queryStr.find(each, startIndex)
            indexLists = []

            # 找到所有出现的位置，返回位置的列标
            while startIndex <= strLen - i and index != -1:
                indexLists.append(index)
                startIndex = index + i
                index = queryStr.find(each, startIndex)

            # 存放连续index的列表
            continueLists = []

            # 检查列表，是否有连续的
            if len(indexLists) > 1:
                # print(each, indexLists)

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

    return repeDict


def savePlot(qualSeries, repeDict, cnt):

    # # 绘制初始质量图
    fig = plt.figure(figsize=(13, 6.5))
    ax1 = fig.add_subplot(111)
    ax1.scatter(qualSeries.index, qualSeries.values, s=5)
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
    fig.savefig('result_img/test%s.jpg' % cnt)
    plt.close()

if __name__ == "__main__":

    frontNum = 5
    STARTCUT = 10
    ENDCUT = 30
    bf = pysam.AlignmentFile('bams/sampleChr20_sorted.bam', 'rb')
    cnt = -1


    for r in bf:
        cnt = cnt + 1

        qname = r.qname

        # 获得bp串
        read = r.query_sequence

        # 切除开头和结尾的部分碱基
        readCut = read[STARTCUT:len(read) - ENDCUT]

        # 找到重复区域的index
        repeDict = FindRepe(read)
        # print(repeDict)

        # 获得数值质量
        qual = r.query_qualities
        qualSeries = pd.Series(qual)



        allRepeList = []

        # savePlot(qualSeries, repeDict, cnt)

        for key in repeDict.keys():
            continueLists = repeDict[key]

            for lists in continueLists:
                if len(lists) >= 10:
                    print(cnt, lists)


        # break


    bf.close()

