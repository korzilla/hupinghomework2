import sys,os

cdata = {}
adata = {}
fdata = []
current_key = None
def process_case_line(line):
    global current_key

    d = line.split(",")
    c = d[1]
    current_key = c
    cdata[current_key] = []

def process_vote_line(line):
    d = line.split(",")
    c = d[1]
    cdata[current_key].append(c)

def process_attribute_line(line):
    d = line.split(",")
    adata[d[1]] = [d[3], d[4]]


def prepare_data(filename):
    pass
    #f = open("anonymous-msweb.data", "r")
    f = open(filename, "r")
    if not f:
        print("open data file failed")
        return None
    lines = f.readlines()
    f.close()
    print(f"number of data lines: {len(lines)}")

    for line in lines:
        if len(line) == 0:
            continue
        if line[0] == 'C':
            process_case_line(line)

        elif line[0] == 'V':
            process_vote_line(line)
        elif line[0] == 'A':
            process_attribute_line(line)

    fdata = cdata.values()
    print(f"number of attributes: {len(adata)}")
    print(f"number of cases: {len(cdata)}")
    #print(fdata)
    return fdata

#==========================================================================
# 根据项目，寻找支持度计数
def xunzhao(s, L):
    k = len(s) - 1  # 根据项目长度，决定在频繁几项集中寻找
    if k < 0:  # s为空的情况
        print(s, '不是有效输入')
    # 标志变量，如果s中的项目有和当前比较的不同的，就置0；
    t = 0
    for i in range(len(L[k])):  # 遍历频繁k项集
        t = 1
        for j in range(k + 1):  # 开始比较是不是该项目
            if L[k][i][0][j] != s[j]:
                t = 0
                break
        if t:
            return L[k][i][1]  # 是的话，返回该项目的支持度计数
    print('未找到')  # 遍历结束的话说明未找到
    return -1  # 返回 -1 (无)
 
 
# 得到划分，并顺便求出置信度
def huafen(L, l, conf):  # 划分得到置信度
    X = l[0]  # 需要划分的列表
    lx = len(X)  # 长度
    # 用二进制的性质求真子集
    for i in range(1, 2 ** lx - 1):
        s1 = []
        s2 = []
        for j in range(lx):
            # 二进制末尾是0就进s1,否则就进s2达到划分目的
            if (i >> j) % 2:
                s1.append(X[j])
            else:
                s2.append(X[j])
        conf[str(s1), '->', str(s2)] = l[1] / xunzhao(s1, L)
 
 
# 判断s中的数据在data中的数目
def jishu(s, data):
    c = 0  # 记录出现数目
    # 标志变量，如果s中的数据有在data这一行不存在的，就置0；
    t = 0
    for ii in data:  # 数据每一行
        t = 1
        for jj in s:  # 对于 s 中的每一个项
            if not (jj in ii):
                t = 0
                break
        c += t  # 如果 s 在这一行存在，c++
    return c
 
 
# 输入频繁k-1项集，支持度计数，数据，输出频繁k项集
def Apriori(L, N, data):
    if not L:
        return L
    L_ = []  # 频繁k项集
    L = [i[0] for i in L]  # k-1项集包含哪些
    k = len(L[0]) + 1  # 几项集
    L_len = len(L)
    up = k - 2  # 拼接同项长度
    i = 0
    while i < L_len - 1:  # 只剩最后一个时肯定没法拼
        A = L[i][0:up]  # 拼接前项
        c = i  # 记录走到那个前项了
        i += 1  # i 到下一个
        for j in range(c + 1, L_len):
            if L[j][0:up] != A:  # 前几项不一致就停止
                i = j  # i 快进到发现新键值的地方
                break
            else:  # 前几项一致时
                s = L[c] + L[j][up:]  # 生成预选项
                t = jishu(s, data)  # 得到 s 的支持度计数
                if t >= N:  # 支持度计数大于N
                    L_.append([s, t])  # 添加到频繁项集中
    return L_
 
 
#data = data_fetch()  # 生成数据，存到data
data = prepare_data("anonymous-msweb.data")
maxk = 1
for d in data:
    if len(d) > maxk:
        maxk = len(d)
print(f"maxk: {maxk}")
#N = 2  # 支持度计数
N = 2000
L = {}  # 储存一项集
for i in data:  # 得到所有的一项集
    for j in i:
        if j in L:
            L[j] += 1
        else:
            L[j] = 1
L_1 = []  # 储存频繁一项集
for i in L.keys():  # 取其中所有的频繁集
    if L[i] >= N:
        L_1.append([[i], L[i]])
for  k,v in L_1:
    print(f"{k}:{v}")
L = []  # 频繁1~4项集
L.append(sorted(L_1, key=lambda x: x[0]))  # 按键值排序转为列表
#print('频繁 1 项集', L[0])
for i in range(3):
    L.append(Apriori(L[i], N, data))  # 求频繁2,3,4项集
    #print('频繁', i + 2, '项集', L[i + 1])
print(L[1])
print(L[2])
print(L[3])
# 置信度计算 ===============================================
conf = {}  # 储存置信度
# 对每一个k项集分析
print('置信度：============================================')
for i in L[1:len(L)]:  # 频繁1项集不用看
    for j in i:
        huafen(L, j, conf)  # 划分得到置信度
for key in conf:
    #print(key,conf[key]*100,'%')
    print(f"{key},{conf[key]*100:.2f}%")

def calculate_lift(L, l, lift):
    X = l[0]  # 需要划分的列表
    lx = len(X)  # 长度
    # 用二进制的性质求真子集
    for i in range(1, 2 ** lx - 1):
        s1 = []
        s2 = []
        for j in range(lx):
            # 二进制末尾是0就进s1,否则就进s2达到划分目的
            if (i >> j) % 2:
                s1.append(X[j])
            else:
                s2.append(X[j])
        lift[str(s1), '->', str(s2)] = (l[1] / xunzhao(s1, L)) / (xunzhao(s2, L) / len(data))
 
    
lift = {}
for i in L[1:len(L)]:  # 频繁1项集不用看
    for j in i:
        calculate_lift(L, j, lift)
print('提升度：============================================')
for key in lift:
    #print(key,lift[key]*100,'%')
    print(f"{key},{lift[key]*100:.2f}%")

kulc = {}
def calculate_kulc(L, l, kulc):
    X = l[0]  # 需要划分的列表
    lx = len(X)  # 长度
    # 用二进制的性质求真子集
    for i in range(1, (2 ** (lx - 1))):
        s1 = []
        s2 = []
        for j in range(lx):
            # 二进制末尾是0就进s1,否则就进s2达到划分目的
            if (i >> j) % 2:
                s1.append(X[j])
            else:
                s2.append(X[j])
        a1 = l[1] / xunzhao(s1, L)
        a2 = l[1] / xunzhao(s2, L)

        print(f"{str(s1)} -> {str(s2)}, kulc= {100*(a1+a2)/2:.2f}%")
        #lift[str(s1), '->', str(s2)] = (l[1] / xunzhao(s1, L)) / (xunzhao(s2, L) / len(data))
        #conf[str(s1), '->', str(s2)] = l[1] / xunzhao(s1, L)
print('Kulc：============================================')
#a = [[1, 2], 15]
#calculate_kulc(L, a, kulc)
#sys.exit(0)
for i in L[1:len(L)]:  # 频繁1项集不用看
    for j in i:
        calculate_kulc(L, j, kulc)
#==========================================================================
#prepare_data()
