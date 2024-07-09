import json
import sys
import os
import time

min_size = int(sys.argv[-3])
max_size = int(sys.argv[-2])
step = int(sys.argv[-1])

try:
    with open('bench_data_grid.json','r') as f:
        bench_data = json.load(f)
except:
    bench_data = {}

for lat_size in range(min_size, max_size, step):
    if str(lat_size)+'x'+str(lat_size)+'.gr' not in os.listdir():
        os.system('python3 produce_grid_instances.py '+str(lat_size))
    
    if lat_size not in bench_data.keys():
        t0 = time.time()
        os.system('java -Xmx10g -Xms10g -Xss10m tw.exact.MainDecomposer < '+str(lat_size)+'x'+str(lat_size)+'.gr > '+str(lat_size)+'x'+str(lat_size)+'.td')
        bench_data[lat_size] = time.time() - t0
    print("lat_size: ", lat_size, "rt", bench_data[lat_size])

for key, val in bench_data.items():
    print(key, ':', val)
