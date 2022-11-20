=begin
# MIT LICENCE

Copyright © 2022 hirooarata

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=end

# A program of the FFT in Ruby with complex/matrix gem
# version 5.0 2022/11/06 19:30 Python, BUG
# version 5.1 2022/11/06 19:30 Python, OK
# version 5.2 2022/11/20 21:30 Ruby, OK

include Math
require 'complex'
require 'matrix'

def dif_fft_ruby(flg, x)
    #      ver 1.0
    #       dif_fft_ruby (frequeny-inplace type)
    #         The algorithm of FFT is referred from  "Discrete-Time Signal Processing"
    #         by A.V.Oppenheim, R.W.Schafer, p.599
    #         prentice-hall,Englewood cliffs, New Jersey 07632
    #         ISBN 013216292X
    #       --- flg= -1.0 .. regular FFT transform
    #       ---    = +1.0 .. inverse FFT transform
    #       --- x[n] .. data (complex)
    n = x.size
    nu = Math.log2(n).to_i
    # p "dif_fft:n=", n, "nu=log2(n)=", nu,"\n"

    # butterfly
    nu.times do | ll | # for ll in range(nu):
        le = 2 ** (nu - ll)
        le1 = le /  2
        u= Complex(1, 0)
        w= Complex( cos(PI / le1) ,  sin(flg * PI / le1) )
        le1.times do | j |   # for j in range(le1):
            j.step( n - 1, le ) do | i |   # for i in range(j, n, le):
                ip = i + le1
                xp = x[i] + x[ip]
                x[ip] = (x[i] - x[ip]) * u
                x[i] = xp
            end
            u = u * w
        end
    end

    # bit reverse
    nv2 = n /  2
    j = 0
    (n-1).times do | i |  # for i in range(n - 1):
        if i < j then
            x[j], x[i] = x[i], x[j]
        end
        k = nv2
        while k <= j do
            j -= k
            k /= 2
        end
        j += k
    end

    # normalize for ifft
    if flg >= 0 then
        n.times do | i |  # for i in range(n):
            x[i] /= n
        end
    end
    return x
end

#
def fft(x)
    dif_fft_ruby(-1, x)
    return x
end

def ifft(x)
    dif_fft_ruby(1, x)
    return x
end


def fft_shift(x)
    nn2 = x.size
    nn = nn2 / 2
    nn.times do | k |  # for k in range(nn):
        t=x[k]
        x[k]= x[k+nn]
        x[k+ nn]=t
    end
    return x
end

def fft_checking(nn2)
    x = Vector.zero(nn2)
    x *= 0+0i
    x[1] = 1 + 1i
    x[nn2-1] = 1 + 1i
    puts "x=#{x}"
    x0 = x
    y = dif_fft_ruby(1.0, x) #  From frquency domain to time domain.
    q = fft_shift(y)
    y = fft_shift(q)
    puts "y=#{y}"
    z = dif_fft_ruby(-1.0, y) #  From time domain to frquency domain.
    puts "z=#{z}"
    puts "z-x=#{z-x0}"
    puts "up/down:loopback_error=#{(z-x0).norm}"
end

fft_checking(8)

=begin
~ % ruby fft.rb
x=Vector[0+0i, 1+1i, 0+0i, 0+0i, 0+0i, 0+0i, 0+0i, 1+1i]
y=Vector[1/4+1/4i, 0.17677669529663687+0.1767766952966368i, 0.0+0.0i, -0.17677669529663692-0.17677669529663687i, -1/4-1/4i, -0.17677669529663687-0.1767766952966368i, 0.0+0.0i, 0.17677669529663692+0.17677669529663687i]
z=Vector[0.0+0.0i, 0.9999999999999998+1.0i, 0.0+0.0i, 5.551115123125783e-17+0.0i, 0.0+0.0i, 1.6653345369377348e-16+0.0i, 0.0+0.0i, 1.0+1.0i]
z-x=Vector[0.0+0.0i, 0.0+0.0i, 0.0+0.0i, 0.0+0.0i, 0.0+0.0i, 0.0+0.0i, 0.0+0.0i, 0.0+0.0i]
up/down:loopback_error=0.0
=end
