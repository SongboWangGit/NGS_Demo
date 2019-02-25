
def SaveAsFile(path):

    f = open(path, 'r')

    cnt = 0
    sequences = ''
    outFile = open('operFa.py', 'r')
    # 遍历fa文件，存储到字典中
    for line in f:
        if line.startswith(">"):
            outFile.close()
            name = line.replace('>', '').split()[0]
            outFile = open('../ref/%s_list.txt' % name, 'w')
        else:
            outFile.write(str(line.replace('\n', '').strip()))
            # cnt = cnt + 1


if __name__ == '__main__':

    refFile = '/home/sbwang/user/workspace/genone/REF/GRCh38.d1.vd1.fa'

    SaveAsFile(refFile)