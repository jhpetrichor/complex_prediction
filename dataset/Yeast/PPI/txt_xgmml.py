import networkx as nx

# 读取.txt文件中的边列表
with open('krogan_extented.txt', 'r') as file:
    edges = [line.strip().split() for line in file]

# 创建一个空的无向图
graph = nx.Graph()

# 添加边到图中
for edge in edges:
    source = edge[0]
    target = edge[1]
    graph.add_edge(source, target)

# 将图保存为.xgmml文件
nx.write_gml(graph, './xgmml/krogan_extented.xgmml')