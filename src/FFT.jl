#----------------------------------------------------
# FFTパッケージ構造体
module FFTC
    # dif_fft関数
    function dif_fft(flg::Float64, x::Vector{Complex{T}}) where T
        n = length(x)
        nu = floor(Int, log2(n))
        # butterfly
        for nm in 0:nu-1
            le = 2^(nu - nm)
            le1 = le ÷ 2
            u = 1 + 1im * 0
            w = cos(pi / le1) + 1im * sin(flg * pi / le1)
            for j in 0:le1-1
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
        for i in 0:n-2
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
    #------------------------------------------------
    # fft関数
    function fft(x::Vector{Complex{T}}) where T
        x = dif_fft(-1.0, x)
        return x
    end
    #------------------------------------------------
    # ifft関数
    function ifft(x::Vector{Complex{T}}) where T
        x = dif_fft(1.0, x)
        return x
    end
    #------------------------------------------------
    # fft_shift関数
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
    #------------------------------------------------
    # ifft_shift関数
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
    #------------------------------------------------
    # fft_checking関数
    function fft_checking(nn2::Int)
        x = zeros(Complex{T}, nn2)
        y = zeros(Complex{T}, nn2)
        z = zeros(Complex{T}, nn2)
        x[2] = 1 + 1im
        x0 = copy(x)
        y = dif_fft(-1.0, x)
        z = dif_fft(1.0, y)
        println("\n\nup/down:loopback_error=", norm(z - x0))
    end
    #------------------------------------------------
end # module FFT
