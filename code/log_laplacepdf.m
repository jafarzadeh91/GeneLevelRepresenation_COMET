function density = log_laplacepdf(x, mu, sudo_cov)
    density = (-1*abs(x-mu)/sudo_cov)-log(2*sudo_cov);

end