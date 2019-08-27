using BenchmarkTools

function func1(a::Array{Float64,4})

   #div = foo(a)
   @time div = zeros(size(a))

   @time div .= foo(a)

   a .= a + div

end

function func2(a::Array{Float64,4})

   let div = ones(10,10,10,8)
      #div = foo(a)
      a .= a + div
      #div[1] += 1
      #@show div[1]
   end

end

function foo(a::Array{Float64,4})
   div = a .+ 1
end

function foo(a::Array{Float64,1})
   div = a .+ 1
end

const func3 = let y = zeros(100)
    function inner_foo(x)
        for i in 1 : 100
            y[i] = 2*x[i] - 1
        end
        return sum(y)
    end
end

function test()
   a = zeros(10,10,10,8)
   #=
   @time for i=1:10
      func1(a)
   end
   =#
   #=
   @time for i=1:10
      func2(a)
   end
   =#
   b = zeros(100)
   @time for i=1:10
      func3(b)
   end
end
