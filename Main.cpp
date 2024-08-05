#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>

// Function to compute the mean of a vector
double mean(const std::vector<double>& data) {
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    return sum / data.size();
}

// Function to compute the variance of a vector
double variance(const std::vector<double>& data, double mean) {
    double var = 0.0;
    for (double value : data) {
        var += (value - mean) * (value - mean);
    }
    return var / (data.size() - 1);
}

// Function to compute the skewness of a vector
double skewness(const std::vector<double>& data, double mean, double variance) {
    double skew = 0.0;
    double var_pow_1_5 = std::pow(variance, 1.5);

    for (double value : data) {
        skew += std::pow(value - mean, 3);
    }

    skew /= (data.size() * var_pow_1_5);
    return skew;
}

// Function to calculate the pairwise correlation coefficients
std::vector<double> calculatePairwiseCorrelations(const std::vector<std::vector<double>>& returns) {
    size_t n = returns[0].size();
    size_t m = returns.size();

    std::vector<double> correlations;

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = i + 1; j < m; ++j) {
            // Calculate correlation between returns[i] and returns[j]
            double mean_i = mean(returns[i]);
            double mean_j = mean(returns[j]);
            double var_i = variance(returns[i], mean_i);
            double var_j = variance(returns[j], mean_j);

            double cov = 0.0;
            for (size_t k = 0; k < n; ++k) {
                cov += (returns[i][k] - mean_i) * (returns[j][k] - mean_j);
            }
            cov /= (n - 1);

            double correlation = cov / (std::sqrt(var_i) * std::sqrt(var_j));
            correlations.push_back(correlation);
        }
    }

    return correlations;
}

int main() {
    // Example returns data (rows: assets, columns: time periods)
    std::vector<std::vector<double>> returns = {
        {0.01, 0.02, -0.01, 0.03, 0.02},
        {0.02, 0.01, 0.00, 0.04, 0.01},
        {-0.01, 0.00, 0.01, 0.02, -0.01}
    };

    std::vector<double> correlations = calculatePairwiseCorrelations(returns);

    double mean_corr = mean(correlations);
    double var_corr = variance(correlations, mean_corr);
    double skew = skewness(correlations, mean_corr, var_corr);

    std::cout << "Mean Correlation: " << mean_corr << std::endl;
    std::cout << "Variance of Correlations: " << var_corr << std::endl;
    std::cout << "Skewness of Correlations: " << skew << std::endl;

    return 0;
}
