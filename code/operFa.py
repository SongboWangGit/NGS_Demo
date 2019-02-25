

def SaveAsDict(path):

    f = open(path, 'r')

    cnt = 0
    sequences = {}
    # 遍历fa文件，存储到字典中
    for line in f:
        # cnt += 1
        # print(cnt)
        # if cnt == 100:
        #     break
        if line.startswith('>'):
            name = line.replace('>', '').split()[0]
            sequences[name] = ''
        else:
            sequences[name] += line.replace('\n', '').strip()

    # 字典的本地存储
    dict = open('../genome/chr20_dict.txt', 'w')
    dict.write(str(sequences))
    return sequences['chr20']

if __name__ == '__main__':

    refFile = '/home/sbwang/user/workspace/genone/REF/GRCh38.d1.vd1.fa'

    SaveAsDict(refFile)