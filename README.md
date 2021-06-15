# Economic-Scenario-Generator
With the current low yield environment, life insurance company are required more and more transparency of their products as well as continuously presenting new investment opportunity. Insurers need to have the pre-contractual risk indicators in order to reflex the possible performance and risks of the products which the subscribers wish to know and could also assist the launch of potential products. The determination of these risk indicators for these products involves stochastic modelling and simulation of different economic scenarios in a market consistent or risk-neutral world. This process is often time-consuming and require very complex mathematics. Under numerous regulations, Insurers and asset managers need to have a standard method of calculation process with the help of computer software to ease the computational heavy. An economic scenario generator (or ESG) which is a computer software can be a key element tool for this process.

As a market consistent ESG, the ESG models have to be calibrated to market prices or volatilities, quoted on the evaluation date, of financial instrument that can reflect the nature and term of the insurer’s liabilities. Those are the two principal constraints imposed by the regulator for constructing a market consistent ESG:
	The market consistency criterion: the aim of market consistent scenarios is to reproduce market prices and implied volatilities.
	The martingale criterion: if there is no arbitrage opportunity, the discounted prices of financial assets are martingales under the risk-neutral probability [13].

	Principal component analysis
The principal component analysis (or PCA) is the most popular method when dealing with yield curve. It allows us to take a set of yield curves, using linear mathematical methods, and then define a reduced form model for the yield curve. This reduced form model retains only a small number of principal but preserved most variance (or information) of the set of yield curve and therefore can reproduce the majority of yield curves. In other words, using PCA, yield curves can be approximated by linear combination of a small number of principal components, with relative small errors. The principal components and the output of the reduced model represent the same quantity at the yield curves set analyzed [6].
The main idea behind PCA mechanism is that a high dimensional system can be approximated to a reasonable degree of accuracy by a system with a smaller number of dimensions by exploiting correlations between the system variables. The PCA provides a natural way of simplifying the statistical model whilst maintaining its ability to reproduce the majority of yield curves that can be produced by the structural model, by simply reducing the number of principal components that are used [12]. This reduced model includes only the quantity of principal components that are deemed important, as indicated by the corresponding eigenvalues. The PCA method is required whenever an interest rate or interest rates are observed at multiple times in the future. Another purpose of the PCA when applied to yield curve is to capture the correlation between interest rate movements at different points of time in the future. Refer to appendix 1 for the mathematics behind PCA.

	Bootstrapping 
Efron (1979) [10] introduced the bootstrap method as a means of estimating the sampling distribution of statistics based on independent observations. 
Let be i.i.d. random variables X with the cumulative distribution function F. Let X_n= (X1,X2,…,Xn). Given a specified random variable R(X_n, F), possibly depending on both X_n and unknown distribution F, we want to estimate the sampling distribution of R based on X_n itself. 
	Construct the empirical distribution function F_n(x) = n^(-1) ∑_(i=1)^n▒I (Xi < x) of X_n.
	With F_n fixed, draw a random sample of size m = m_n with replacement from F_n  
X_i^*~ F_n , i = l, 2,.., m. Call this the bootstrap sample, and let X_n^*= (X_1^*,X_2^*,…,X_m^*)
	Approximate the sampling distribution of R(X_n, F) by the conditional distribution of R(X_n^*,F_n) given X_n (called the bootstrap distribution of R (X_n, F)), where n and m are called the sample size and the resample size respectively.
 
In comparison with Monte Carlo simulation, which is a way of simulating prices over time is by choosing a priori assumption regarding the distribution and trying to adapt it to the stock to be simulated, a method to simulate a stock price without doing any a priori assumptions about the distribution of the log return is non-parametric bootstrap. As not always the assumption about the distribution is true, for instance the assumption of normally-distributed log return in the Black and Schole framework in the case of exotic option pricing with consideration of the volatility smile. The idea of bootstrapping is that given a sample, which is a good approximation of the population, randomly drawing observations from it creating a large number of resamples on which some statistic is calculated. The relative frequency histogram of the statistic reveals a good estimate of its distribution. 

	Interest rate models 
There has been numerous researches and developments of interest rate modeling. Interest rate models can be broadly categorized into Short rate models, forward rate models, Libor market models and Swap market models. Each type of these classes of models targets the calibration of different observed/unobserved interest rate instrument in the market. For examples, short rate models describe the dynamics of the instantaneous spot rate while forward rate models choose the instantaneous forward rate to calibrate. On the other hand, Libor and Swap market models try to model the evolution of rates which are directly observed in the market: libor rates and swap rates. The short rate class can be divided further into 2 types of model: equilibrium models and no-arbitrage models. The term “endogenous term structure” usually implies equilibrium models as the term structure of interest rate in an output of these model. The most common used this class of model is Vasicek model (1977) and Cox-Ingersoll-Ross model (1985). The idea of no-arbitrage models is to have the exact match of current term structure of interest rate. The Hull-White model (1990) and the Black-Karasinski model (1991) are two representatives for this class of models [5].
The first criteria to select an interest rate model for the market consistent economic scenario generator is to take the current term structure of interest rate as an input, this rules out the equilibrium models or endogenous term structure models [6]. The second priority to consider is model calibration to market prices should be relatively simple and not much computational heavy. In this paper, we choose the Hull-White one-factor short rate model as it satisfies both conditions.

2.4.1 Vasicek model
Introduced in 1977 by Oldřich Vašíček, the Vasicek model was the first one to capture mean reversion, an essential characteristic of the interest rate that sets it apart from other financial prices [4]. As it is an equilibrium model, it does not take the current term structure as the input, rather it uses a period of historical data to reflects the endogenous term structure of interest rate.
This does not satisfy the market consistent or risk neutral property. We included this model in our calculation to show the different with Hull-White model.
In Vasicek’s model, the risk-neutral process for r is
dr(t)=a(b-r(t))dt+ σdW(t),r(0)=r_0
where:
	r_t is the instantaneous short rate in risk-neutral world where, in a very short time period between t and t+∆t, investors earn on average r_t+∆t. The process for r in the real world is irrelevant.
	r_0 is a positive constant.
	a is the reversion rate of the process while b is the long term mean of instantaneous rate.
	σ(t) is the volatility of the instantaneous short rate.
	W(t) is a Wiener process. 
a, b, and σ are constants . This model incorporates mean reversion. The short rate is pulled to a level b at rate a. Vasicek equation can be used to obtain the following expression for the price at time t of a zero-coupon bond that pays $1 at time T:
P(t,T)= A(t,T) e^(-B(t,T)r(t))
In this equation r(t) is the value of instantaneous rate r at time t,
B(t,T) =  (1- e^(-a(T-t)))/a
A(t,T)=exp⁡[((B(t,T)-T+t)(a^2 b-σ^2/2))/a^2 -  σ^2/4a 〖B(t,T)〗^2]
The yield of zero coupon bond maturity T at time t is calculated by:
R(t,T) =  (-ln⁡(P(t,T))/(T-t)

2.4.2 Hull-White one factor model
In a paper published in 1990, Hull and White explored extensions of the Vasicek model that provide an exact fit to the initial term structure [3]. The Hull-White one-factor model assumes that the instantaneous short-rate process r is normally distributed and has dynamics given by
dr(t)=(θ(t)-a(t)r(t))dt+ σ(t)dW(t),r(0)=r_0
where: 
	r_t is the instantaneous short rate in risk-neutral world where, in a very short time period between t and t+∆t, investors earn on average r_t+∆t. The process for r in the real world is irrelevant.
	r_0 is the spot rate at very short maturity (1month), represent the instantaneous short rate. 
	a(t) r(t) θ(t) is deterministic functions of time, a(t) is rate at which r(t) reverts towards its expected value.
	σ(t) is the volatility of the instantaneous short rate.
	W(t) is a Wiener process. 
As Hull and White also developed a two factor model, which is not under the scope of this paper, let us denote in this paper the Hull-White one factor model as Hull-White model. Since some or all parameters of Hull-White are allowed to vary with time, this ensures that the model can match the current term structure of interest rates exactly, which the Vasicek model cannot do. The advantage of Hull-White model is having explicit close-form formulas for zero-coupon bond. Therefore the calibration process is straightforward. Also the normally distributed short rate process can produce negative interest rates, which can reflect current market situation. This not possible for Black-Karasinski model [7]. 
Calibration of Hull-White model means determining the functions a(t) and σ(t) to fit prices of interest rate derivatives on the market, and the function θ(t)  to match the exact current yield curve. Hull and White (1994) introduced a method to assume that both a(t) and σ(t) are time-independent. By fixing a and σ, the model behaves like Vasicek model to some extend, thus the name Hull-White model [2]. However, the disadvantage of this approach is that the future volatility shall behave differently with today volatility.
Setting a(t) = a and σ(t) = and σ as 2 positive constants and θ(t) is chosen so that the model exactly matches the current term structure of interest rates, we shall have:
dr(t)=(θ(t)-ar(t))dt+ σdW(t),r(0)=r_0
θ(t)=  (∂f^M (0,t))/∂t+ 〖af〗^M (0,t)+  σ^2/2a(1- e^(-2at))
where (∂f^M (0,T))/∂T  is the partial derivative of f^M with respect to its second argument, and f^M (0,T) is the market instantaneous forward rate at time 0 for maturity T.
f^M (0,T)= -(∂lnP^M (0,T))/∂T
where P^M (0,t)  is the zero-coupon price for maturity T.
By integrating equation, we get: 
r(t)=r(s) e^(-a(t-s) )+ ∫_s^t▒〖e^(-a(t-u) ) θ(u)du〗+σ∫_s^t▒〖e^(-a(t-u) ) dW(u)〗
=r(s) e^(-a(t-s) )+α(t)-α(s) e^(-a(t-s) )+σ∫_s^t▒〖e^(-a(t-u) ) dW(u)〗
where,
α(t)= f^M (0,T)+  σ^2/(2a^2 ) 〖(1-e^(-at))〗^2
Bond prices at time t in the Hull–White model are given by
P(t,T)= A(t,T) e^(-B(t,T)r(t))
where P(t, T) is the price at time t of a zero-coupon bond with maturity T
B(t,T) =  (1- e^(-a(T-t)))/a
ln⁡A(t,T)=ln (P(0,T))/(P(0,t))+B(t,T) f^M (0,t)-  1/(4a^3 ) σ^2 (e^2at-1)(e^(-aT)-e^(-at) )^2
or: 
A(t,T)=(P(0,T))/(P(0,t))*exp⁡[B(t,T) f^M (0,t)-  σ^2/4a (1-e^(-2at) ) 〖B(t,T)〗^2]
Finally, it is straight forward to deduce the yield of a zero coupon bond maturity T at time t:
R(t,T) =  (-ln⁡(P(t,T))/(T-t)
