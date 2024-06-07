import math

def norm_cdf(x):
    """
    Calculate the cumulative distribution function (CDF) of the standard normal distribution.
    
    :param x: Value at which to calculate the CDF
    :return: CDF value
    """
    return (1.0 + math.erf(x / math.sqrt(2.0))) / 2.0

def black_scholes_barrier_option(S, K, T, r, sigma, H, option_type, barrier_type):
    """
    Calculate the price of a European barrier option using the Black-Scholes formula.
    
    :param S: Current price of the underlying asset
    :param K: Strike price of the option
    :param T: Time to expiration (in years)
    :param r: Risk-free interest rate (annual)
    :param sigma: Volatility of the underlying asset (annual)
    :param H: Barrier level
    :param option_type: 'call' for call option, 'put' for put option
    :param barrier_type: 'in' for knock-in, 'out' for knock-out
    :return: Option price
    """
    d1 = (math.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * math.sqrt(T))
    d2 = d1 - sigma * math.sqrt(T)
    lambda_ = (r + 0.5 * sigma**2) / sigma**2
    x1 = math.log(S / H) / (sigma * math.sqrt(T)) + lambda_ * sigma * math.sqrt(T)
    y1 = math.log(H / S) / (sigma * math.sqrt(T)) + lambda_ * sigma * math.sqrt(T)
    y = (math.log(H**2 / (S * K)) + (r + 0.5 * sigma**2) * T) / (sigma * math.sqrt(T))
    
    if option_type == 'call' and barrier_type == 'out':
        call_price = (S * norm_cdf(d1) - K * math.exp(-r * T) * norm_cdf(d2)) - \
                     (H * norm_cdf(y1 - sigma * math.sqrt(T)) - H * norm_cdf(-y1 + sigma * math.sqrt(T)) - 
                      S * norm_cdf(-d1 + 2 * lambda_ * sigma * math.sqrt(T)) + K * math.exp(-r * T) * norm_cdf(-d2 + 2 * lambda_ * sigma * math.sqrt(T)))
        return call_price
    elif option_type == 'call' and barrier_type == 'in':
        call_price = H * norm_cdf(y1 - sigma * math.sqrt(T)) - H * norm_cdf(-y1 + sigma * math.sqrt(T)) - \
                     S * norm_cdf(-d1 + 2 * lambda_ * sigma * math.sqrt(T)) + K * math.exp(-r * T) * norm_cdf(-d2 + 2 * lambda_ * sigma * math.sqrt(T))
        return call_price
    elif option_type == 'put' and barrier_type == 'out':
        put_price = (K * math.exp(-r * T) * norm_cdf(-d2) - S * norm_cdf(-d1)) - \
                    (H * norm_cdf(-y1 + sigma * math.sqrt(T)) - H * norm_cdf(y1 - sigma * math.sqrt(T)) - 
                     S * norm_cdf(d1 - 2 * lambda_ * sigma * math.sqrt(T)) + K * math.exp(-r * T) * norm_cdf(d2 - 2 * lambda_ * sigma * math.sqrt(T)))
        return put_price
    elif option_type == 'put' and barrier_type == 'in':
        put_price = H * norm_cdf(-y1 + sigma * math.sqrt(T)) - H * norm_cdf(y1 - sigma * math.sqrt(T)) - \
                    S * norm_cdf(d1 - 2 * lambda_ * sigma * math.sqrt(T)) + K * math.exp(-r * T) * norm_cdf(d2 - 2 * lambda_ * sigma * math.sqrt(T))
        return put_price
    else:
        raise ValueError("Invalid option or barrier type. Use 'call' or 'put' for option_type and 'in' or 'out' for barrier_type.")

# Exemplo de uso:
S = 12.74  # Preço atual do ativo subjacente
K = 14.1414  # Preço de exercício da opção
T = 282/252  # Tempo até a expiração (em anos)
r = 6.101856  # Taxa de juros livre de risco (anual)
sigma = 0.2745  # Volatilidade do ativo subjacente (anual)
H = 15.28  # Nível de barreira

# Calculando diferentes tipos de opções de barreira
call_in_price = black_scholes_barrier_option(S, K, T, r, sigma, H, 'call', 'in')
call_out_price = black_scholes_barrier_option(S, K, T, r, sigma, H, 'call', 'out')
put_in_price = black_scholes_barrier_option(S, K, T, r, sigma, H, 'put', 'in')
put_out_price = black_scholes_barrier_option(S, K, T, r, sigma, H, 'put', 'out')

print(f"Call In Option Price: {call_in_price:.2f}")
print(f"Call Out Option Price: {call_out_price:.2f}")
print(f"Put In Option Price: {put_in_price:.2f}")
print(f"Put Out Option Price: {put_out_price:.2f}")
