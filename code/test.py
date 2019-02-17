def SaveAsList(path):

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

if __name__ == '__main__':

    # SaveAsList('../genome/hs_ref_GRCh38.p12_chr20.fa')
    file = open('../genome/chr20.txt')

    for line in file.readlines():
        print(len(line))
        # hhhhhhhh