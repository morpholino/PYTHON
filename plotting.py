import pandas as pd
brooklyn_one_bed = pd.read_csv('brooklyn-one-bed.csv') #reading tables into pd.array
june = london_data.loc[london_data["month"] == 6]["TemperatureC"] #get temperature for filtered data (june)
#transactions is a pd.dataframe
transactions = transactions.drop(["Unnamed: 0"], axis = 1) #drop a column from data
times = transactions["Transaction Time"].values # data to a separate numpy array
times = transactions["Transaction Time"] #keep column and row info
times = transactions["Transaction Time"].unique() #list of unique items in column
times = transactions[transactions["Transaction Time"] <= whatever] #filter dataframes
times_hist = np.histogram(times, range = (0, 24), bins = 4) #bin data as a histogram
#histogram produces two arrays, first array contains counts per bin, second contains the min (max) values per bin

import numpy as np
#common characteristics of data:
#centrality, modality, spread (range), skew, outliers
average_age = np.average(author_ages) #must be an array
median_age = np.median(author_ages) #must be an array; if the array has even number of items, 
									#it is two values/average of two mid values
#right-skew if average is greater than median, left-skew if opposite, or symmetric distro
example_sum = np.sum(example_array)
example_variance = np.var(example_array) #diff from average squared
example_std = np.std(example_array) #sqrt of variance, i.e. ** 0.5
example_min = np.amin(times) #amax for array max

from scipy import stats
example_mode = stats.mode(example_array) #most frequent value in data
#unimodal has one peak, bimodal/multimodal more peaks; no peaks - uniform
str(example_mode[0][0]) #the value
str(example_mode[1][0]) #its count

data_q3 = np.quantile(example_array, 0.75) #upper quartile value
data_quantiles = np.quantile(example_array, [0.2, 0.4, 0.6, 0.8]) #quintile values; 1.0 not needed
#interquartile range (IQR) ignores the outliers, so you know the range around which your data is centered
from scipy.stats import iqr
interquartile_range = iqr(dataset) #or stats.iqr(dataset)

from matplotlib import pyplot as plt
# Use plt.hist() below
plt.hist(author_ages, range=(10, 80), alpha = 0.75, bins=14,  edgecolor='black') #tzn 10-15,15-20,20-25,...
#alpha >> transparency
#if more histograms need be plotted, just repeat plt.hist with another data
plt.title("Age of Top 100 Authors at Publication")
plt.axvline(average_age, color='r', linestyle='solid', linewidth=2, label="Mean") #add vertical line
#dotted and dashed lines ale alternatives
plt.xlabel("Age")
plt.ylabel("Count")
#if you have variability in data, use ribbon plot:
plt.fill_between(data, y_lower, y_upper, alpha=0.2)

plt.show()

#if two subplots are needed, then start with 
plt.figure(1) #might not be necessary
plt.figure(figsize=(20,6)) #if you need to define the size of the result
plt.subplot(211)
whatever plot

plt.subplot(212)
whatever plot below that

plt.tight_layout() #for plots less squished

#boxplot is cool because it plots automatically several datasets side-by-side
plt.boxplot([dataset_one, dataset_two], labels = ["Label 1", "Label 2"])
#!!whiskers extend from the box 1.5 times the size of the IQR, outliers are shown as dots
#to quickly fill boxplot datasets from pandas:
states = chest_pain["Provider State"].unique()
datasets = []
for state in states:
  datasets.append(chest_pain[chest_pain['Provider State'] == state][' Average Covered Charges '].values)
plt.boxplot(datasets, labels=states)
plt.show()