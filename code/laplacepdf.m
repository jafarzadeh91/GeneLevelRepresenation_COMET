function density = laplacepdf(x, mu, sudo_cov)
    density = exp(-1*abs(x-mu)/sudo_cov)/(2*sudo_cov);

end