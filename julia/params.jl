module Params

export TransferParams, AccuracyParams, CAMBparams, default_params

struct TransferParams
    high_precision::Bool
    accurate_massive_neutrinos::Bool
    kmax::Float64
    k_per_logint::Int
    PK_num_redshifts::Int
    PK_redshifts::Vector{Float64}
end

TransferParams(; high_precision=false, accurate_massive_neutrinos=false,
                kmax=0.9, k_per_logint=0, PK_num_redshifts=1,
                PK_redshifts=Float64[]) =
    TransferParams(high_precision, accurate_massive_neutrinos,
                   kmax, k_per_logint, PK_num_redshifts, PK_redshifts)

struct AccuracyParams
    AccuracyBoost::Float64
    lSampleBoost::Float64
    lAccuracyBoost::Float64
    AccuratePolarization::Bool
    AccurateBB::Bool
    AccurateReionization::Bool
    TimeStepBoost::Float64
    BackgroundTimeStepBoost::Float64
    IntTolBoost::Float64
    SourcekAccuracyBoost::Float64
    IntkAccuracyBoost::Float64
    TransferkBoost::Float64
    NonFlatIntAccuracyBoost::Float64
    BessIntBoost::Float64
    LensingBoost::Float64
    NonlinSourceBoost::Float64
    BesselBoost::Float64
    LimberBoost::Float64
    SourceLimberBoost::Float64
    KmaxBoost::Float64
    neutrino_q_boost::Float64
end

AccuracyParams(; AccuracyBoost=1.0, lSampleBoost=1.0, lAccuracyBoost=1.0,
                AccuratePolarization=true, AccurateBB=false,
                AccurateReionization=true, TimeStepBoost=1.0,
                BackgroundTimeStepBoost=1.0, IntTolBoost=1.0,
                SourcekAccuracyBoost=1.0, IntkAccuracyBoost=1.0,
                TransferkBoost=1.0, NonFlatIntAccuracyBoost=1.0,
                BessIntBoost=1.0, LensingBoost=1.0, NonlinSourceBoost=1.0,
                BesselBoost=1.0, LimberBoost=1.0, SourceLimberBoost=1.0,
                KmaxBoost=1.0, neutrino_q_boost=1.0) =
    AccuracyParams(AccuracyBoost, lSampleBoost, lAccuracyBoost,
                   AccuratePolarization, AccurateBB, AccurateReionization,
                   TimeStepBoost, BackgroundTimeStepBoost, IntTolBoost,
                   SourcekAccuracyBoost, IntkAccuracyBoost, TransferkBoost,
                   NonFlatIntAccuracyBoost, BessIntBoost, LensingBoost,
                   NonlinSourceBoost, BesselBoost, LimberBoost,
                   SourceLimberBoost, KmaxBoost, neutrino_q_boost)

struct CAMBparams
    WantCls::Bool
    WantTransfer::Bool
    WantScalars::Bool
    WantTensors::Bool
    WantVectors::Bool
    WantDerivedParameters::Bool
    Want_cl_2D_Array::Bool
    Want_CMB::Bool
    Want_CMB_lensing::Bool
    DoLensing::Bool
    NonLinear::Int
    Transfer::TransferParams
    want_zstar::Bool
    want_zdrag::Bool
    Min_l::Int
    Max_l::Int
    Max_l_tensor::Int
    Max_eta_k::Float64
    Max_eta_k_tensor::Float64
    ombh2::Float64
    omch2::Float64
    omk::Float64
    omnuh2::Float64
    H0::Float64
    TCMB::Float64
    Yhe::Float64
    Num_Nu_massless::Float64
    Num_Nu_massive::Int
    Nu_mass_eigenstates::Int
    share_delta_neff::Bool
    Nu_mass_degeneracies::Vector{Float64}
    Nu_mass_fractions::Vector{Float64}
    Nu_mass_numbers::Vector{Int}
    Accuracy::AccuracyParams
end

function default_params()
    transfer = TransferParams()
    acc = AccuracyParams()
    return CAMBparams(true, false, true, false, false, true, true,
                      true, true, true, 0, transfer, false, false, 2,
                      2500, 600, 5000.0, 1200.0, 0.0, 0.0, 0.0, 0.0, 67.0,
                      2.7255, 0.24, 3.044, 0, 0, false,
                      zeros(5), zeros(5), zeros(Int,5), acc)
end

end # module
