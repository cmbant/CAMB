module Config

export version, FeedbackLevel, output_file_headers, DebugMsgs, DebugEvolution,
       DebugParam, ThreadNum, do_bispectrum, hard_bispectrum,
       full_bessel_integration, OutputDenominator, C_Temp, C_E, C_Cross,
       C_Phi, C_PhiTemp, C_PhiE, CT_Temp, CT_E, CT_B, CT_Cross,
       C_name_tags, CT_name_tags, lens_pot_name_tags, OmegaKFlat, base_tol,
       highL_unlensed_cl_template, lmax_extrap_highl, highL_CL_template,
       lensed_convolution_margin, global_error_flag, global_error_message,
       error_reionization, error_recombination, error_inital_power,
       error_evolution, error_unsupported_params, error_darkenergy, error_ini,
       error_nonlinear, error_allocation, global_error

using ..Constants: const_twopi

version = "1.6.1"

FeedbackLevel = 0
output_file_headers = true
DebugMsgs = false
const DebugEvolution = false
DebugParam = 0.0
ThreadNum = 0

do_bispectrum = false
const hard_bispectrum = false
const full_bessel_integration = false
const OutputDenominator = const_twopi

const C_Temp = 1
const C_E = 2
const C_Cross = 3
const C_Phi = 4
const C_PhiTemp = 5
const C_PhiE = 6

const CT_Temp = 1
const CT_E = 2
const CT_B = 3
const CT_Cross = 4

const C_name_tags = ["TT","EE","TE","PP","TP","EP"]
const CT_name_tags = ["TT","EE","BB","TE"]
const lens_pot_name_tags = ["TT","EE","BB","TE","PP","TP","EP"]

const OmegaKFlat = 5e-7
const base_tol = 1.0e-4
highL_unlensed_cl_template = "HighLExtrapTemplate_lenspotentialCls.dat"
const lmax_extrap_highl = 8000
highL_CL_template = Array{Float64}(undef, 0, 0)
const lensed_convolution_margin = 100

global_error_flag = 0
global_error_message = ""

const error_reionization = 1
const error_recombination = 2
const error_inital_power = 3
const error_evolution = 4
const error_unsupported_params = 5
const error_darkenergy = 6
const error_ini = 7
const error_nonlinear = 8
const error_allocation = 9

function global_error(; message="", id=nothing)
    global global_error_message
    global global_error_flag
    global_error_message = message
    if id !== nothing
        if id == 0
            error("Error id must be non-zero")
        end
        global_error_flag = id
    else
        global_error_flag = -1
    end
end

end # module
