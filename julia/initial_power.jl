module InitialPower

export tensor_param_indeptilt, tensor_param_rpivot, tensor_param_AT,
       InitialPowerLaw, scalar_power, tensor_power, effective_ns

using ..Constants: const_pi

const tensor_param_indeptilt = 1
const tensor_param_rpivot = 2
const tensor_param_AT = 3

struct InitialPowerLaw
    ns::Float64
    nrun::Float64
    nrunrun::Float64
    nt::Float64
    ntrun::Float64
    r::Float64
    pivot_scalar::Float64
    pivot_tensor::Float64
    As::Float64
    At::Float64
    tensor_parameterization::Int
    curv::Float64
end

function InitialPowerLaw(; ns=1.0, nrun=0.0, nrunrun=0.0, nt=0.0, ntrun=0.0,
                           r=0.0, pivot_scalar=0.05, pivot_tensor=0.05,
                           As=1.0, At=1.0,
                           tensor_parameterization=tensor_param_indeptilt,
                           curv=0.0)
    return InitialPowerLaw(ns, nrun, nrunrun, nt, ntrun, r,
                           pivot_scalar, pivot_tensor, As, At,
                           tensor_parameterization, curv)
end

function scalar_power(ip::InitialPowerLaw, k::Float64)
    lnrat = log(k / ip.pivot_scalar)
    return ip.As * exp(lnrat * (ip.ns - 1 + lnrat * (ip.nrun / 2 + ip.nrunrun / 6 * lnrat)))
end

function tensor_power(ip::InitialPowerLaw, k::Float64)
    lnrat = log(k / ip.pivot_tensor)
    k_dep = exp(lnrat * (ip.nt + ip.ntrun / 2 * lnrat))
    if ip.tensor_parameterization == tensor_param_indeptilt
        val = ip.r * ip.As * k_dep
    elseif ip.tensor_parameterization == tensor_param_rpivot
        val = ip.r * scalar_power(ip, ip.pivot_tensor) * k_dep
    else
        val = ip.At * k_dep
    end
    if ip.curv < 0
        val *= tanh(const_pi / 2 * sqrt(-k^2 / ip.curv - 3))
    end
    return val
end

effective_ns(ip::InitialPowerLaw) = ip.ns

end # module
