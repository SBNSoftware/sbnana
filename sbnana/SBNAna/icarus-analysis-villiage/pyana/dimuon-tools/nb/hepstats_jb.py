import numpy as np
from scipy.stats import poisson,norm
from scipy.optimize import minimize
from scipy.special import loggamma
from scipy.interpolate import interp1d

# Define this by hand as the scipy version doesn't play well with real-valued
# arguments (as opposed to integer)
def poisson_logpdf(n,mu):
    return -mu + n * np.log(mu) - loggamma(n + 1.0)

# Defines a simple 1 bin model with a benchmark (MC) signal s, CV background b
# and relative uncertainties on the signal and background sig_s and sig_b
# Uncertainties are taken to be Gaussian
class Model:
    def __init__(self,s,b,sig_s,sig_b):
        # Record the model parameters
        self.s = s
        self.b = b
        self.sig_s = sig_s
        self.sig_b = sig_b
        # Cache the time-intensive results for the model to save time
        self.max_L_cache = {}
        self.profile_L_cache = {}
        self.CLs_cache = {}
        self.CLsb_cache = {}
        self.CLb_cache = {}
        # Store samples from the likelihood at various signal strengths
        self.samples = {}

    # Log likelihood given data
    # The likelihood in the bin is taken to be Poisson
    # The normalization uncertainthy is from a Gaussian (truncated at 0)
    def ll(self,n,pars):
        mu = pars[0]
        gamma_s = pars[1]
        gamma_b = pars[2]
        res = poisson_logpdf(n,mu * gamma_s * self.s + gamma_b * self.b)
        res = res + norm.logpdf(gamma_s,1.0,self.sig_s)
        res = res + norm.logpdf(gamma_b,1.0,self.sig_b)
        return res

    # Twice the negative log likelihood, used in constructing test stat
    def twice_nll(self,n,pars):
        return -2.0 * self.ll(n,pars)

    # Twice negative log of maximum likelihood
    def max_L_fit(self,n):
        if n in self.max_L_cache:
            return self.max_L_cache[n]
        # If best fit signal will be negative, take best LL to be at signal = 0
        if n - self.b < 0:
            res = self.profile_L_fit(0.0,n)
            return np.insert(res[0],0,0.0),res[1]
        else:
            return [(float(n) - self.b) / self.s,1.0,1.0], self.twice_nll(n,[(float(n) - self.b) / self.s,1.0,1.0])

    # Twice negative log of profile likelihood
    def profile_L_fit(self, mu_val, n):
        if (mu_val,n) in self.profile_L_cache:
            return self.profile_L_cache[(mu_val,n)]
        # Minimize the negative log likelihood
        # Can be done analytically, but it's a mess, so let's do it numerically
        init_vals = np.array([1.0,1.0])
        bound_vals = (1.0e-10,10.0),(1.0e-10,10.0)
        minres = minimize(lambda pars : self.twice_nll(n,[mu_val,pars[0],pars[1]]),init_vals,bounds = bound_vals)
        self.profile_L_cache[(mu_val,n)] = minres.x,minres.fun
        # Beware if the resulting signal or background scales are negative
        if minres.x[0] < 0 or minres.x[1] < 0:
            print('Warning: expected signal or background less than zero')
        return self.profile_L_cache[(mu_val,n)]

    # The test statistic qmu_tilde
    # If the signal is less than the best fit, give 0 (better for upper limits)
    # If the signal is greater than best fit, it's the twice the negative
    # log of the profile likelihood ratio
    # Note: for a normal distribution, this would just be chi^2
    def test_stat(self,mu_val,n):
        gfit, prof_L = self.profile_L_fit(mu_val,n)
        mugfit, max_L = self.max_L_fit(n)
        #return prof_L - max_L
        if mugfit[0] <= mu_val:
            return prof_L - max_L
        else:
            return 0.0

    # Generate and store smaples from the likelihood
    def generate_samples(self,mu_val):
        self.samples[mu_val] = np.empty((10000,),dtype=np.int64)
        for i in range(10000):
            # The distributions for signal and background depend on only 1 variable
            # Sample these first from a normal distribution
            gs = norm.rvs(1.0,self.sig_s)
            # Delete negative cases: events should be positive
            while gs < 0.0:
                gs = norm.rvs(1.0,self.sig_s)
            gb = norm.rvs(1.0,self.sig_b)
            while gb < 0.0:
                gb = norm.rvs(1.0,self.sig_b)
            # Finally, sample the signal from the Poisson distribution
            n_samp = poisson.rvs(mu_val * gs * self.s + gb * self.b)
            self.samples[mu_val][i] = n_samp

    # Get the p-value
    def pvalue(self, mu_val, n):
        tot_pass = 0
        test_stat_obs = self.test_stat(mu_val,n)
        # Reset samples, just in case
        if mu_val not in self.samples:
            self.generate_samples(mu_val)
        # Take our samples and determine what fraction have test statistic
        # larger than observed
        # That's the p-value!
        for n_samp in self.samples[mu_val]:
            if self.test_stat(mu_val,n_samp) >= test_stat_obs:
                tot_pass = tot_pass + 1
        return float(tot_pass) / 10000.0

    # Get the p-value (bkg)
    def pvalue_bkg(self, mu_val, n):
        tot_pass = 0
        # For background, the test stat should be compared with that for no signal
        test_stat_obs = self.test_stat(0.0,n)
        # Reset samples, just in case
        if 0.0 not in self.samples:
            self.generate_samples(0.0)
        for n_samp in self.samples[0.0]:
            if self.test_stat(mu_val,n_samp) >= test_stat_obs:
                tot_pass = tot_pass + 1
        return float(tot_pass) / 10000.0

    # Get CLb: the p-value for the background model
    def CLb(self,mu_val,n):
        if n not in self.CLb_cache:
            self.CLb_cache[(mu_val,n)] = self.pvalue_bkg(mu_val,n)
        return self.CLb_cache[(mu_val,n)]

    # Get CLsb: the p-value for signal + background model
    def CLsb(self,mu_val,n):
        if (mu_val,n) not in self.CLsb_cache:
            self.CLsb_cache[(mu_val,n)] = self.pvalue(mu_val,n)
        return self.CLsb_cache[(mu_val,n)]

    # Get CLs: the ratio CLsb / CLb
    def CLs(self,mu_val,n):
        if (mu_val,n) not in self.CLs_cache:
            self.CLs_cache[(mu_val,n)] = self.CLsb(mu_val,n) / self.CLb(mu_val,n)
        return self.CLs_cache[(mu_val,n)]

    # Expected CLs
    def CLs_exp(self,mu_val):
        # Evaluate expectation at 0 signal
        if 0.0 not in self.samples:
            self.generate_samples(0.0)
        # Get the CLs
        CLs_samp = np.asarray([self.CLs(mu_val,n_samp) for n_samp in self.samples[0.0]])
        # We use the normal distribution percentiles to determine the 1 sigma and 2 sigma bands
        # This is probably not strictly correct, but it's what "everyone" does
        normal_percentiles = np.asarray([2.27501319, 15.86552539, 50.0, 84.13447461, 97.72498681])
        return [np.percentile(CLs_samp,perc) for perc in normal_percentiles]
