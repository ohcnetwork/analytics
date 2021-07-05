# What is the effective reproduction rate of a pandemic, R<sub>t</sub> ?

The effective reproduction rate of a pandemic, R<sub>0</sub>, is defined as is the average number of secondary cases per infectious case in a population made up of both susceptible and non-susceptible hosts. However, in COVID-19 pandemic, a static R0 does not adequately reflect the reality in time and space due to changes in social behaviors and restrictions. 

Computation of Rt for COVID-19 enables understanding how effectively a local or state government handles the pandemic and gives the authority helpful information in decision to loosen and tighten measures of social restrictions. As the pandemic spreads with great acceleration, R<sub>t</sub> is much more larger than 1. On the contrary, when the pandemic slows down and dies out, Rt is smaller than 1 and approaches 0. This project focuses on computation of Rt for every state on the Kerala based on the number of new cases k reported daily by the state's Department of Health. The value of R<sub>t</sub> is related to that of a day before R<sub>t-1</sub>, and every previous value of n days before, R<sub>t-n</sub>.

## Algorithms used in this project:

### Bettencourt & Ribeiro
We will use the same approach outlined by Kevin in his article, which uses the paper [1] as its underlying basis.

### Modeling arrivals
The first step is to model the 'arrival' process of infections. A popular choice for the distribution of arrivals amongst statisticians is the Poisson Distribution. Accordingly, if we let λ represent the average rate of infections per day, then the probability that we are likely to see k new cases on a day, is given by
   
P(k|λ) = λ<sup>k</sup> * e<sup>−λ</sup> / k!

Given this setup, we can construct the probability distribution of new cases for a set of λs.

### Poisson likelihood
Modeling the arrival process as a Poisson distribution allows us to predict the distribution of new cases in a day as a function of the arrival rate λ. However, in reality, we only observe the number of arrivals. So the key question now is how do we go from the observed number of arrivals to get a sense of the distribution of λ. Thankfully, the answer to this question is simple.

L(λ|k) = λ<sup>k</sup> * e<sup>−λ</sup> / k!

The distribution of λ over k is called the likelihood function, and it has the same mathematical expression as the probability mass function we used earlier. We can visualize the likelihood function by fixing the number of new cases observed (k), and computing the likelihood function over a range of values of λ.

### Relating λ and Rt
According to this paper by Bettencourt & Ribeiro, the relationship between arrival rate λ and effective reproduction rate is defined as follows:

λ = k<sub>t−1</sub> * e<sup>γ(Rt−1)</sup>

Note that γ here is the reciprocal of the serial interval (about 4 days for COVID19), and k<sub>t−1</sub> is the number of new cases observed in the time interval t−1.

We can use this expression for λ and reformulate the likelihood function in terms of R<sub>t</sub>.

L(R<sub>t</sub>|k) = λ<sup>k</sup> * e <sup>−λ</sup> / k!

# Dataset
For this project, we use the data from https://keralastats.coronasafe.live/histories.json page where daily counts of new COVID-19 cases are reported for every district in the state of Kerala. The data is then cleaned and wrangled in a proper dataframe containing the daily count of each district. We select Ernakulam to compute the district's effective reproduction rate of the COVID-19 pandemic, R<sub>t</sub>. Every district's R<sub>t</sub> can be computed at users' choice by modifying the vector districts in the Analytics markdown file.

# Computation Steps

The process to compute Rt can be briefly described as follows:

1. Import the all districts' daily counts to a dataframe
2. Initialize a value of γ and a set of discrete values of R<sub>t</sub>
3. Select one or more states of interest
4. Smooth out the daily counts to flatten the choppy data points
5. Compute the log-likelihood distribution P(k|R<sub>t</sub>)
6. Compute the posterior P(R<sub>t</sub>|k<sub>t</sub>)
7. Compute the estimate of R<sub>t</sub>, including the most-likely, the max, and the min values of R<sub>t</sub>
8. Add Lockdown data to understand its impact


# References
[1] https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0002185
[2] https://www.datacamp.com/community/tutorials/replicating-in-r-covid19
[3] https://github.com/calldrj/COVID19.Effective.Reproduction.Rate