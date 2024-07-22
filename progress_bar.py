from tqdm import tqdm
import numpy as np
import time

def clamp(minimum, x, maximum):
    return max(minimum, min(x, maximum))

class ProgressBar:
    def __init__(self, total = 100, desc = "Progress", bar = True):
        self.bar = bar
        self.desc = desc
        self.total = total
        self.count = 0
        self.last_time = time.perf_counter_ns()
        if(self.bar):
            self.pbar = tqdm(total = self.total, desc = self.desc)
    def update(self, res, cv):
        if(self.bar):
            self.pbar.n =  round(clamp(0, cv, 1)*100.0, 2)
            self.pbar.refresh()
        else:
            print(f"{round(clamp(0.0, cv, 1.0) * 100.0, 2)}% - Current best solution -> {np.round(res, 2)}")
    def update_xk(self, x):
        self.count+=1
        #every 30 calls, calculate residual and display progress
        if(self.bar):
            self.pbar.n = self.count
            self.pbar.refresh()
        elif(time.perf_counter_ns() - self.last_time > 2e9):
            print(f"Current best solution -> {np.round(x, 2)}")
            self.last_time = time.perf_counter_ns()
    def update_da(self, x, f, context):
        self.count+=1
        #every 30 calls, calculate residual and display progress
        if(self.bar):
            self.pbar.n = self.count
            self.pbar.refresh()
        elif(not self.bar and time.perf_counter_ns() - self.last_time > 2e9):
            print(f"New minima's residual: {round(f, 2)} - New minima: {np.round(x, 2)}")
            self.last_time = time.perf_counter_ns()
        return False
    def close(self):
        if(self.bar):
            self.pbar.n = self.total
            self.pbar.refresh()
            self.pbar.close()