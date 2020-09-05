function noise_level = gnoise(sigma)
% try using 0.01 times the RMS intensity to be the sigma value 
  noise_level = ((rand + rand + rand + rand + rand + rand + rand + rand +rand + rand + rand + rand) - 6)*sigma; 
end