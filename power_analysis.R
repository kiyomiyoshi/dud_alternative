library(pwr)
pwr.t.test(d = 0.7867221, n = 24, sig.level = 0.05, 
           type =  "one.sample", alternative = "two.sided")

pwr.t.test(d = 0.7867221, n = 19, sig.level = 0.05, 
           type =  "one.sample", alternative = "greater")