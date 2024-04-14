#-------------------------------------------------
# Transcode from Python to Julia.
#-------------------------------------------------
module FFT
#-------------------------------------------------
export dif_fft,fft,ifft,fft_shift,ifft_shift
export fft_checking,test_FFT
#-------------------------------------------------
function dif_fft(flg::Float64, x::Vector{Complex{T}}) where T
    n = length(x)
    nu = floor(Int, log2(n))
    # butterfly
    for nm in 1:nu
        le = 2 ^ (nu - nm +1 )
        le1 = le ÷ 2
        println("le=",le,", le1=",le1)
        u = 1 + 1im * 0
        w = cos(pi / le1) + 1im * sin(flg * pi / le1)
        for j in 1:le1
            for i in j:le:n-1
                ip = i + le1
                xp = x[i] + x[ip]
                x[ip] = (x[i] - x[ip]) * u
                x[i] = xp
            end
            u *= w
        end
    end
    # bit reverse
    nv2 = n ÷ 2
    j = 0
    for i in 1:n-1
        if i < j
            x[j], x[i] = x[i], x[j]
        end
        k = nv2
        while k <= j
            j -= k
            k ÷= 2
        end
        j += k
    end
    # normalize for ifft
    if flg >= 0
        for i in 1:n
            x[i] /= n
        end
    end
    return x
end
#-------------------------------------------------
# fft function
function fft(x::Vector{Complex{T}}) where T
    x = dif_fft(-1.0, x)
    return x
end
#-------------------------------------------------
# ifft function
function ifft(x::Vector{Complex{T}}) where T
    x = dif_fft(1.0, x)
    return x
end
#-------------------------------------------------
# fft_shift function
function fft_shift(x::Vector{Complex{T}}) where T
    nn2 = length(x)
    c = div(nn2, 2)
    if c % 2 == 0
        for k in 1:c
            x[k + c], x[k] = x[k], x[k + c]
        end
    else
        tmp = x[1]
        for k in 1:c
            x[k] = x[c + k + 1]
            x[c + k + 1] = x[k + 1]
        end
        x[c] = tmp
    end
    return x
end
#-------------------------------------------------
# ifft_shift function
function ifft_shift(x::Vector{Complex{T}}) where T
    nn2 = length(x)
    nn = div(nn2, 2)
    if nn2 % 2 == 0
        for k in 1:nn
            x[k], x[k + nn] = x[k + nn], x[k]
            x[k + nn], x[k] = x[k], x[k + nn]
        end
    else
        tmp = x[nn2]
        for k in nn-1:-1:1
            x[nn + k + 1] = x[k]
            x[k] = x[nn + k]
        end
        x[nn] = tmp
    end
    return x
end
#-------------------------------------------------
norm(x)=sqrt(sum(abs2, x))
#-------------------------------------------------
function test_FFT()
    nn2=4
    T=Float64
    x = zeros(Complex{T}, nn2)
    y = zeros(Complex{T}, nn2)
    z = zeros(Complex{T}, nn2)
    x[1] = 1 + 1im
    x0 = copy(x)
    y = FFT.dif_fft(-1.0, x)
    z = FFT.dif_fft(1.0, y)
    println("\n\nup/down:loopback_error=", norm(z - x0))
    l=  norm(z - x0) < 1e-4
    println("test FFT.test_FFT = ", l)
    return l
end
#-------------------------------------------------
end # module FFT
#-------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    FFT.test_FFT()
end
