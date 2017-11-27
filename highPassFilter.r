# A FFT-based high-pass filter for a vector $x.
# $f is the threshold frequency in % of maximum frequency.
# This implementation stuffs both ends of $x to avoid overshooting due to the discrete FT

mgHighPassFilter <- function(x, f)
{
    n = length(x)

    # allow no negative values to prevent wrapping during inverse FFT
    x_min = min(x)
    x = x - x_min

    # values of end points used for stuffing,
    # averaging over some more endpoints could improve robustness
    x0 = x[1]
    x1 = x[n]

    # stuff spectrum to prevent swings at spectral ends
    x_stuffed <- c(rep(x0, n), x, rep(x1, n));

    # go to frequency domain
    xf = fft(x_stuffed)

    # cut-off: given f and f_nyquist
    f0 = 3 * n / 2 / 100 * f
    f1 = 3 * n / 2
    
    # cut out the high frequencies harshly. try a softer transition?
    xf[f0:f1] = 0
    xf[(3 * n - f1):(3 * n - f0)] = 0
    
    # back to time domain
    x_hp = Mod(fft(xf, inverse = TRUE))
     
    # remove stuffed values   
    x_hp <- x_hp[(n + 1):(2 * n)]

    # correct scaling
    x_hp <- x_hp / n / 3

    # add original offset
    x_hp = x_hp + x_min
    
    return(x_hp)
}
