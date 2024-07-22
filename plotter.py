import matplotlib.pyplot as plt


class Plotter:
    def __init__(self):
        plt.figure(figsize=(4, 3))
    def plot(self, x_data, y_data, alpha, label):
        plt.plot(x_data, y_data, alpha=alpha, label = label)
    def fill_between(self, x, y1, y2=0, alpha=1, label = None):
        plt.fill_between(x, y1, y2, alpha=alpha, label = label)