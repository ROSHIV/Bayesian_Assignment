function b = Bernoulli_f(theta,gama)
b = (theta.^gama)*(1-theta).^(1-gama);
end
