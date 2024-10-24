# Danielle De La Pascua
# 10-24-24
# Monte Carlo R translation

# Compute grand mean and variance of family mean resistance for control treatment


# Compute mean and variance across all families of the difference in family mean resistance in damage vs control plants


# Draw population level mean and variance for constitutive resistance from respective sampling distributions


# Draw a constitutive value for each family equal to the number of populations in the experiment from a log normal dist with the population-level mean and variance 


# Calculate Cmin (smalled value)


# Draw population-level mean (î) and variance for induced resistance from each respective sampling distributions


# draw an induced resistance for each family
# formula I = L { î + Cmin , V } - Cmin 


# Add the induced resistance (which may be negative) to the constitutive resistance for each family to get the expected resistance in the damage treatment


# Draw resistance for the replicate plants in a family from two log normal distributions, one for the control and one for the damage treatment.
# The coefficients of variation for those distributions is estimated by the variability among replicate plants within treatments within the experimental data 


# Finally, compute for all families the mean resistance from the control and damage sample, the difference measure of induced resistance, and the correlation between constitutive and induced resistance. 


# Repeat the entire process many times and identify the lower 5th percentile of the resulting distribution of correlation coefficients as the most negative correlation likely to arise by chance when there is no trade-off.