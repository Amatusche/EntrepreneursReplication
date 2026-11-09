function lab_fun(lab,tau,w,z,rb,b,T,varphi,sigma,gamma,min_c)
    diff = max((1-tau)*w*z*lab+rb*b+T,min_c)^(-sigma*gamma)*((1-tau)*w*z/varphi)^(gamma) - lab;

end
