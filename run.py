from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess

maxProcess = 8
runCount = 100
exePath = R"C:\Users\huww\source\homeworks\tsp_ga\Release\tsp_ga.exe"
outPutDir = R".\results"
def run(seed):
    with open(outPutDir + "\\output_{}.txt".format(seed), "w") as f:
        subprocess.check_call([exePath, str(seed)], stdout=f)

with ThreadPoolExecutor(max_workers=maxProcess) as executor:
    futures = {executor.submit(run, seed): seed for seed in range(runCount)}
    for future in as_completed(futures):
        seed = futures[future]
        try:
            data = future.result()
            print('completed {}'.format(seed))
        except Exception as exc:
            print('error {}'.format(seed))
