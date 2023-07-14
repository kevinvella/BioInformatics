import argparse
import matplotlib.pyplot as plt
import pyperf
import pylab
import scipy.stats as stats

def display_histogram_scipy(bench, mean, bins):
    values = bench.get_values()
    values = sorted(values)
    if mean:
        fit = stats.norm.pdf(values, bench.mean(), bench.stdev())
        pylab.plot(values, fit, '-o', label='mean-stdev')
    else:
        fit = stats.norm.pdf(values, bench.mean(), bench.stdev())
        pylab.plot(values, fit, '-o', label='mean-stdev')

    plt.legend(loc='upper right', shadow=True, fontsize='x-large')
    pylab.hist(values, bins=bins)
    pylab.show()

def main():
    
    
    suite = pyperf.BenchmarkSuite.load("benchmark_pigeonHoleTextSearch.json")
    bench = suite.get_benchmark("benchmark_pigeonHoleTextSearch")
    
    display_histogram_scipy(bench, '', 25)
    xxxx = ""

if __name__ == "__main__":
    main()