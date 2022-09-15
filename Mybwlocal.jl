

function bwlocallinear1(xdata::RealVector, ydata::RealVector, kernel::Function=gaussiankernel)
    n=length(xdata)
    length(ydata)==n || error("length(ydata) != length(xdata)")
    w = ones(n)
    if kernel == gaussiankernel
        h0= bwnormal(xdata)
        hlb = 0.1*h0
        hub = 10*h0
    elseif kernel == betakernel
        h0 = midrange(xdata)
        hlb = h0/n
        hub = 0.25
    elseif kernel == gammakernel
        h0 = midrange(xdata)
        hlb = h0/n
        hub = h0
    elseif kernel == epkernel
      h0 = midrange(xdata)
      hlb = 0.1*h0
      hub = 10*h0
    end
    Optim.minimizer(Optim.optimize(h->AIClocallinear(xdata, ydata, kernel,h, w, n), hlb, hub))
end
