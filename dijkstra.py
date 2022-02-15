import json
import math
from heapq import heappush, heappop, heapify
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

rn = 0.4 #筐体の半径
acc_round = 5 #角の演算精度(直角に曲がってる点を分割するため)

# 衝突検知として、外積より線分の交差を検知
# https://qiita.com/zu_rin/items/e04fdec4e3dec6072104
def closses(va,vb,vc,vd):
    # (AB*AC)
    s=(vb[0]-va[0])*(vc[1]-va[1]) - (vc[0]-va[0])*(vb[1]-va[1])
    # (AB*AD)
    t=(vb[0]-va[0])*(vd[1]-va[1]) - (vd[0]-va[0])*(vb[1]-va[1])
    if s*t <= 0:
        # (CD*CA)
        s=(vd[0]-vc[0])*(va[1]-vc[1]) - (va[0]-vc[0])*(vd[1]-vc[1])
        # (CD*CB)
        t=(vd[0]-vc[0])*(vb[1]-vc[1]) - (vb[0]-vc[0])*(vd[1]-vc[1])
        if s*t <= 0:
            return True
    return False

def nlen(n1,n2=(0.0,0.0)):#ベクトルの長さを返す
    return math.sqrt((n1[0]-n2[0])**2+(n1[1]-n2[1])**2)

def nnormalize(n,l=1.0):#ベクトルを長さlに正規化する
    len = nlen(n)/l
    return (n[0]/len,n[1]/len)

def nangle(n1,n2=(1.0,0.0)):#ベクトルの角度を返す
    abs1 = nlen(n1)
    abs2 = nlen(n2)
    ext = n1[0]*n2[0]+n1[1]*n2[1]#内積
    angle = math.acos(ext/(abs1*abs2))

    if n1[0]*n2[1]-n1[1]*n2[0]<0: #外積のz成分がマイナス ->角度がマイナス
        return -angle
    return angle

#https://ikatakos.com/pot/programming_algorithm/route_search/dijkstra
def trace(s, t, ancestors):
    route = [t]
    c = t
    while True:
        a = ancestors[c]
        assert a is not None, 'Failed to trace'
        route.append(a)
        if a == s:
            break
        c = ancestors[c]
    route.reverse()
    return route
def dijkstra(nodes={(2, 0): []}, start=(0,0),goal=(0,0)):
    heap = [(nodes[start][l],l, start) for l in nodes[start]]

    heapify(heap)
    visited = set()
    ancestors = {}
    while heap:
        cost, node, ancestor = heappop(heap)
        if node in visited:
            continue
        visited.add(node)
        ancestors[node] = ancestor
        if node == goal:
            return cost, trace(start, goal, ancestors)
        for node2 in nodes[node]:
            if node2 not in visited:
                heappush(heap, (cost + nodes[node][node2], node2, node))
    return float('inf'), None

#構造物のリストから通る可能性のあるパスを生成する
def genpath(structs,start,goal):
    nodes = {}
    nodes[start] = {} #Start
    nodes[goal] = {} #Goal
    lines = []

    for i in range(len(structs)):
        struct = structs[i]
        edges = []
        for line in struct: #構造物の辺
            if not line['s'] in edges:
                edges.append(line['s'])
            if not line['e'] in edges:
                edges.append(line['e'])
            lines.append(line)#ラインのリストに追加
        center = (0,0)
        for eg in edges: #構造物の頂点
            center = (center[0]+eg[0],center[1]+eg[1])
        center = (center[0]/len(edges),center[1]/len(edges))
        for line in range(len(struct)): #構造物の辺
            s0 = struct[(line-1)%len(struct)]['s']
            s1 = struct[line]['s']
            e1 = struct[line]['e']
            if "i" in struct[line]:continue

            if not "c" in struct[line]:#直線
                # norm = nnormalize((e1[0]-s1[0],e1[1]-s1[1]),rn)
                nvect = (0,0)#法線ベクトル
                fromcenter = (s1[0]-center[0],s1[1]-center[1])

                if abs(e1[0]-s1[0])<abs(e1[1]-s1[1]):#y軸方向(0,1)or(0,-1)
                    nvect = nnormalize((fromcenter[0],0),rn)
                else:#x軸方向
                    nvect = nnormalize((0,fromcenter[1]),rn)

                if not "i" in struct[line]:
                    nodes[(s1[0]+nvect[0],s1[1]+nvect[1])] = {}
                if not "i" in struct[(line+1)%len(struct)]:
                    nodes[(e1[0]+nvect[0],e1[1]+nvect[1])] = {}
                
                n1 = nnormalize((s1[0]-s0[0],s1[1]-s0[1]))#s0 to s1
                n2 = nnormalize((e1[0]-s1[0],e1[1]-s1[1]))#s1 to e1
                if abs(nangle(n1,n2)) > math.pi/3:
                    theta = nangle(n2,n1)
                    theta0 = nangle((0,1),n1)
                    for j in range(acc_round+1):
                        v = (rn*math.cos(theta0-j*theta/acc_round),rn*math.sin(theta0-j*theta/acc_round))
                        nodes[(s1[0]+v[0],s1[1]+v[1])] = {}
                    # way = nnormalize((n1[0]-n2[0],n1[1]-n2[1]),rn)
                    # nodes[(s1[0]+way[0],s1[1]+way[1])] = {}

                    
            else:#扇形
                c1 = struct[line]['c']
                r1 = nlen((s1[0]-c1[0],s1[1]-c1[1]))+rn
                sr = (s1[0]-c1[0],s1[1]-c1[1])#中心をc1としたs1
                er = (e1[0]-c1[0],e1[1]-c1[1])#中心をc1としたe1

                theta = nangle(er,sr)
                theta0 = nangle((1,0),sr)

                for j in range(acc_round+1):
                    v = (r1*math.cos(theta0-j*theta/acc_round),r1*math.sin(theta0-j*theta/acc_round))
                    nodes[(c1[0]+v[0],c1[1]+v[1])] = {}
    
    for n1 in nodes.keys(): #適当な点1
        for n2 in nodes.keys(): #適当な点2
            if n2[0] == 0 or n2[1] == 0:#もしx･y軸を通るならやめる
                continue
            cl = False #辺と交差してるか
            if n1 == n2:
                continue
            for l in lines:
                if closses(n1,n2,l['s'],l['e']):
                    cl = True
            if cl == True:
                continue

            cls = False #経路がどこかと接触しているか
            for i in range(len(structs)):
                struct = structs[i]
                for line in range(len(struct)): #構造物の辺
                    s1 = struct[line]['s']
                    r = rn-0.02
                    center = struct[line]['s']
                    if "c" in struct[line]:
                        center = struct[line]['c']
                        r = nlen((s1[0]-center[0],s1[1]-center[1]))+rn-0.02
                    
                    #https://yttm-work.jp/collision/collision_0006.html
                    A = (n2[0]-n1[0],n2[1]-n1[1]) #start to end
                    norm_A = nnormalize(A)
                    B = (center[0]-n1[0],center[1]-n1[1]) #start to center
                    D = (center[0]-n2[0],center[1]-n2[1]) #end to center
                    if abs(B[0] * norm_A[1] - norm_A[0] * B[1]) >= r:#あたってる可能性なし 外積
                        continue
                    dot01 = B[0] * A[0] + B[1] * A[1]
                    dot02 = D[0] * A[0] + D[1] * A[1]
                    if dot01 * dot02 <= 0:#接触
                        cls = True
                        break
                    if nlen(B) <= r or nlen(D) <= r:#接触
                        cls = True
                        break
            if cls == True:
                continue
            if not n2 in nodes[n1]:
                nodes[n1][n2] = nlen(n1,n2)#コストは距離とする
    return nodes,lines

if __name__ == '__main__':
    file = open("./structs.json" , "r")
    jsondata = json.load(file)
    structs = jsondata["structs"]
    start = tuple(jsondata["start"])
    goal = tuple(jsondata["goal"])

    nodes, lines = genpath(jsondata["structs"],start,goal)#パス生成
    opt_cost, opt_root = dijkstra(nodes,start,goal)#最短経路探索
    print(opt_root)

    fig , ax = plt.subplots()
    ax.set_xlim([0,27])
    ax.set_ylim([0,14])

    ax.text(start[0],start[1],"Start",size=12,verticalalignment="top")
    ax.plot(start[0],start[1],marker='.',color="orange",markersize=30)
    ax.text(goal[0],goal[1],"Goal",size=12,verticalalignment="top")
    ax.plot(goal[0],goal[1],marker='.',color="blue",markersize=30)

    for ln in lines:
        if "c" in ln:
            c1 = ln['c']
            s1 = ln['s']
            e1 = ln['e']
            sr = (s1[0]-c1[0],s1[1]-c1[1])#中心をc1としたs1
            er = (e1[0]-c1[0],e1[1]-c1[1])#中心をc1としたe1
            thetae = math.degrees(nangle((1,0),er))
            thetas = math.degrees(nangle((1,0),sr))
            height = nlen(ln['e'],ln['c'])*2.0
            if nangle(sr,er) < 0: #swap
                tmp = thetas
                thetas = thetae
                thetae = tmp
            c1 = mpatches.Arc( xy=ln['c'] , height=height, theta1=thetas, theta2=thetae,width=height, color="green",lw=1)
            ax.add_patch(c1)
        else:  
            ax.plot((ln['s'][0],ln['e'][0]), (ln['s'][1],ln['e'][1]),color="green",lw=1)

    for i in range(len(opt_root)-1):
        fr = opt_root[i]
        to = opt_root[i+1]
        ax.plot((fr[0],to[0]), (fr[1],to[1]),color="red")

    # for node in nodes.keys(): #path
        # ax.plot(node[0], node[1],marker='.',color="black")
        # for to in nodes[node]:
        #     ax.plot((node[0],to[0]), (node[1],to[1]),marker='.',color="black")

    plt.show()

    # cost = 0
    # for i in range(len(opt_root)-1):
    #     fr = opt_root[i]
    #     to = opt_root[i+1]
    #     cost += nlen(to,fr)
    #     print(cost)


    # print("cost")

    # for i in range(len(opt_root)-2):
    #     fr = opt_root[i]
    #     to = opt_root[i+1]
    #     to2 = opt_root[i+2]
    #     bef = (to[0]-fr[0],to[1]-fr[1])
    #     aft = (to2[0]-to[0],to2[1]-to[1])
    #     ext = aft[0]*bef[1] - aft[1]*bef[0]
    #     if ext < 0:
    #         print(-1)
    #     else:
    #         if ext > 0:
    #             print(1)
    #         else:
    #             print(0)