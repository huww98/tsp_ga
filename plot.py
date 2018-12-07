import matplotlib.pyplot as plt

results = []

for x in range(100):
    xList = []
    yList = []

    curveFilePath = './results/tspgacurve_{}.txt'.format(x)
    pathFilePath = './results/tsppath_{}.txt'.format(x)
    with open(curveFilePath, 'r') as f:
        content = f.read()
        lines = content.splitlines()
        for line in lines:
            if line == '':
                break
            cell = line.split('\t')
            xList.append(int(cell[0]))
            yList.append(int(cell[1]))
    px = []
    py = []
    with open(pathFilePath, 'r') as f:
        content = f.read()
        lines = content.splitlines()
        for line in lines:
            if line == '':
                break
            cell = line.split()
            px.append(int(cell[1]))
            py.append(int(cell[2]))
        px.append(px[0])
        py.append(py[0])
    results.append((xList, yList, px, py))

def resultSortKey(r):
    return r[1][-1]

bestX, bestY, bestPx, bestPy = min(results, key=resultSortKey)
print(bestY[-1])
worstX, worstY, worstPx, worstPy = max(results, key=resultSortKey)
print(worstY[-1])
plt.plot(bestX, bestY, label="best")
plt.plot(worstX, worstY, label="worst")
plt.legend()
plt.savefig("./results/curve.png")
plt.close()

plt.plot(bestPx, bestPy, marker='x')
plt.savefig("./results/best.png")
plt.close()

plt.plot(worstPx, worstPy, marker='x')
plt.savefig("./results/worst.png")
plt.close()

finals = [r[1][-1] for r in results]
plt.hist(finals)
plt.savefig("./results/hist.png")
plt.close()

print(sum([1 if f == bestY[-1] else 0 for f in finals]))
