using JuMP
using Cbc
using Printf

const a_coef = [0.257986,  0.0389016, -0.00021708, 0.0468684, -0.00094284, -0.00034344];
const b_coef = [0.933884, -0.058212,   0.00450036, 0.00243,    0.000486,   -0.001215];
const p_coef = [4.0033, -3.5162, -5.4302,  5.6807, -1.1989, -0.1963, 0.8593];
const r_coef = [0.9804, -2.6207,  3.6708, -2.6975,  0.0446,  1.2533, 0.2494];
const c_coef = [0.222903, 0.313387, 0.46371];
const Q_norm = 95; # kW
const Q_stor = 500; # kWh
const COP = 3.5;
const Cf = 3.915; # kJ/kg-K
const m_ice = 5.24; # kg/s
const T_cw_ch = -5; # degrees C
const T_cw_norm = 4.4; # degrees C
const T_fr = 0; # degrees C

poly(c::Vector{Float64}, x) = c[1] + sum(c[2:end] .* (x .^ Array(1.0:length(c)-1)));

Ψ_12(T_cw, T_out, coef::Vector{Float64}) = coef[1] + coef[2] * T_cw + coef[3] * T_cw^2 + coef[4] * T_out + coef[5] * T_out^2 + coef[6] * T_cw * T_out;
Ψ_3(PLR) = poly(c_coef, PLR); # PLR: partial load ratio
Q_avail(T_cw, T_out) = Q_norm * Ψ_12(T_cw, T_out, a_coef);
P_chiller(q_chiller, T_cw, T_out) = Q_avail(T_cw, T_out) * Ψ_12(T_cw, T_out, b_coef) * Ψ_3(q_chiller / Q_avail(T_cw, T_out));

u_UB(l) = poly(p_coef, l) * m_ice * Cf * (T_fr - T_cw_ch);
u_LB(l) = poly(r_coef, l) * m_ice * Cf * (T_cw_norm - T_fr);

function TESS_scheduling_PL_Bin(energyPrice::Vector{Float64}, reservePrice::Vector{Float64}, T_out::Vector{Float64}, Q_cool::Vector{Float64}, initSoC::Float64, finalSoC::Float64)
    numHours = size(energyPrice, 1);
    if numHours ≠ size(reservePrice, 1) || numHours ≠ size(T_out, 1) || numHours ≠ size(Q_cool, 1)
        throw(DimensionMismatch("The lengths of energyPrice, reservePrice, and Q_cool should match."))
    end

    iobuffer = IOBuffer();

    m = Model(with_optimizer(Cbc.Optimizer, logLevel=0));

    @variable(m, p[1:numHours] ≥ 0); # Power withdrawal from the grid (non-negative), unit: kW
    @variable(m, 0.1 ≤ l[1:numHours + 1] ≤ 0.9); # State of charge (SoC) in [0.1, 0.9]
    @variable(m, r[1:numHours] ≥ 0); # Reserve power (non-negative)
    @variable(m, u[1:numHours]); # Ice-storage charging rate
    @variable(m, u_p[1:numHours] ≥ 0, base_name = "u^+"); # Positive parts of ice-storage charging rate
    @variable(m, u_n[1:numHours] ≤ 0, base_name = "u^-"); # Negative parts of ce-storage charging rate
    @variable(m, b[1:numHours], Bin); # Ice-storage charging/discharging indicator: 1 == Charging, 0 == Discharging
    @variable(m, u_r[1:numHours], base_name = "\\tilde{u}"); # Ice-storage charging rate
    @variable(m, u_p_r[1:numHours] ≥ 0, base_name = "\\tilde{u}^+"); # Positive parts of ice-storage charging rate
    @variable(m, u_n_r[1:numHours] ≤ 0, base_name = "\\tilde{u}^-"); # Negative parts of ce-storage charging rate
    @variable(m, b_r[1:numHours], Bin, base_name = "\\tilde{b}"); # Ice-storage charging/discharging indicator: 1 == Charging, 0 == Discharging

    t_UB = [0.1, 0.6, 0.78, 0.9];
    a_UB = (u_UB.(t_UB[1:end-1]) - u_UB.(t_UB[2:end])) ./ (t_UB[1:end-1] - t_UB[2:end]);
    b_UB = u_UB.(t_UB[1:end-1]) - a_UB .* t_UB[1:end-1];
    @variable(m, l_UB[1:length(t_UB)-1, 1:numHours], base_name = "\\overline{l}"); # Piecewise linearization weights for u_UB
    @variable(m, S_UB[1:length(t_UB)-1, 1:numHours], Bin, base_name = "\\overline{S}"); # Piecewise linearization binary for u_UB

    t_LB = [0.1, 0.3, 0.5, 0.65, 0.78, 0.9];
    a_LB = (u_LB.(t_LB[1:end-1]) - u_LB.(t_LB[2:end])) ./ (t_LB[1:end-1] - t_LB[2:end]);
    b_LB = u_LB.(t_LB[1:end-1]) - a_LB .* t_LB[1:end-1];
    @variable(m, l_LB[1:length(t_LB)-1, 1:numHours], base_name = "\\underline{l}"); # Piecewise linearization weights for u_LB
    @variable(m, S_LB[1:length(t_LB)-1, 1:numHours], Bin, base_name = "\\underline{S}"); # Piecewise linearization binary for u_LB

    t_PLR = [0, 0.5, 1];
    a_Pch_p = zeros(length(t_PLR) - 1, numHours);
    for k = 1:numHours
        a_Pch_p[:, k] = 
            (
                P_chiller.(t_PLR[1:end-1] * Q_avail(T_cw_ch, T_out[k]), T_cw_ch, T_out[k]) 
                - P_chiller.(t_PLR[2:end] * Q_avail(T_cw_ch, T_out[k]), T_cw_ch, T_out[k])
            ) ./ (
                t_PLR[1:end-1] 
                - t_PLR[2:end]
            ) ./ Q_avail(T_cw_ch, T_out[k]);
    end
    b_Pch_p = zeros(length(t_PLR) - 1, numHours);
    for k = 1:numHours
        b_Pch_p[:, k] = P_chiller.(t_PLR[1:end-1] * Q_avail(T_cw_ch, T_out[k]), T_cw_ch, T_out[k]) - a_Pch_p[:, k] .* t_PLR[1:end-1] * Q_avail(T_cw_ch, T_out[k]);
    end
    a_Pch_n = zeros(length(t_PLR) - 1, numHours);
    for k = 1:numHours
        a_Pch_n[:, k] = 
            (
                P_chiller.(t_PLR[1:end-1] * Q_avail(T_cw_norm, T_out[k]), T_cw_norm, T_out[k]) 
                - P_chiller.(t_PLR[2:end] * Q_avail(T_cw_norm, T_out[k]), T_cw_norm, T_out[k])
            ) ./ (
                t_PLR[1:end-1] 
                - t_PLR[2:end]
            ) ./ Q_avail(T_cw_norm, T_out[k]);
    end
    b_Pch_n = zeros(length(t_PLR) - 1, numHours);
    for k = 1:numHours
        b_Pch_n[:, k] = P_chiller.(t_PLR[1:end-1] * Q_avail(T_cw_norm, T_out[k]), T_cw_norm, T_out[k]) - a_Pch_n[:, k] .* t_PLR[1:end-1] * Q_avail(T_cw_norm, T_out[k]);
    end
    @variable(m, q_p[1:length(t_PLR)-1, 1:numHours], base_name = "q^+");
    @variable(m, S_p[1:length(t_PLR)-1, 1:numHours], Bin, base_name = "S^+"); # Piecewise binary variable for q_p
    @variable(m, q_n[1:length(t_PLR)-1, 1:numHours], base_name = "q^-");
    @variable(m, S_n[1:length(t_PLR)-1, 1:numHours], Bin, base_name = "S^-"); # Piecewise binary variable for q_n
    @variable(m, q_p_r[1:length(t_PLR)-1, 1:numHours], base_name = "\\tilde{q}^+");
    @variable(m, S_p_r[1:length(t_PLR)-1, 1:numHours], Bin, base_name = "\\tilde{S}^+"); # Piecewise binary variable for q_p_r
    @variable(m, q_n_r[1:length(t_PLR)-1, 1:numHours], base_name = "\\tilde{q}^-");
    @variable(m, S_n_r[1:length(t_PLR)-1, 1:numHours], Bin, base_name = "\\tilde{S}^-"); # Piecewise binary variable for q_n_r

    @constraint(m, conInitSoc,  l[1] == initSoC);
    if finalSoC ≥ 0
        @constraint(m, conFinalSoc, l[numHours + 1] == finalSoC);
    end

    @constraint(m, conSocDynamics[k = 1:numHours], l[k + 1] - l[k] == u[k] / Q_stor);

    @constraint(m, conThermalCharging[k = 1:numHours], u[k] == u_p[k] + u_n[k]);

    @constraint(m, [k = 1:numHours], l[k] == sum(l_UB[:, k]));
    @constraint(m, [k = 1:numHours], sum(S_UB[:, k]) == 1);
    @constraint(m, [i = 1:length(t_UB)-1, k = 1:numHours], S_UB[i, k] * t_UB[i] ≤ l_UB[i, k]);
    @constraint(m, [i = 1:length(t_UB)-1, k = 1:numHours], l_UB[i, k] ≤ S_UB[i, k] * t_UB[i+1]);
    @constraint(m, conThermalChargingUB[k = 1:numHours], u_p[k] ≤ sum(a_UB .* l_UB[:, k] + b_UB .* S_UB[:, k])); # u_UB(l[k])
    @constraint(m, [k = 1:numHours], u_p[k] ≤ b[k] * 1e10); 

    @constraint(m, [k = 1:numHours], l[k] == sum(l_LB[:, k]));
    @constraint(m, [k = 1:numHours], sum(S_LB[:, k]) == 1);
    @constraint(m, [i = 1:length(t_LB)-1, k = 1:numHours], S_LB[i, k] * t_LB[i] ≤ l_LB[i, k]);
    @constraint(m, [i = 1:length(t_LB)-1, k = 1:numHours], l_LB[i, k] ≤ S_LB[i, k] * t_LB[i + 1]);
    @constraint(m, conThermalChargingLB[k = 1:numHours], u_n[k] ≥ -sum(a_LB .* l_LB[:, k] + b_LB .* S_LB[:, k])); # -u_LB(l[k])
    @constraint(m, [k = 1:numHours], u_n[k] == -(1 - b[k]) * Q_cool[k]); # Constraining the discharging rate u_n[k] to be either 0 or -Q_cool[k]

    # @constraint(m, [k = 1:numHours], u[k] + Q_cool[k] ≥ 0);

    @constraint(m, [k = 1:numHours], Q_cool[k] + u_p[k] == sum(q_p[:, k]));
    @constraint(m, [k = 1:numHours], sum(S_p[:, k]) == 1);
    @constraint(
        m,
        [i = 1:length(t_PLR)-1, k = 1:numHours], 
        S_p[i, k] * t_PLR[i] * Q_avail(T_cw_ch, T_out[k]) ≤ q_p[i, k]
    );
    @constraint(
        m,
        [i = 1:length(t_PLR)-1, k = 1:numHours], 
        q_p[i, k] ≤ S_p[i, k] * t_PLR[i + 1] * Q_avail(T_cw_ch, T_out[k])
    );
    @constraint(m, [k = 1:numHours], Q_cool[k] + u_n[k] == sum(q_n[:, k]));
    @constraint(m, [k = 1:numHours], sum(S_n[:, k]) == 1);
    @constraint(
        m,
        [i = 1:length(t_PLR)-1, k = 1:numHours], 
        S_n[i, k] * t_PLR[i] * Q_avail(T_cw_norm, T_out[k]) ≤ q_n[i, k]
    );
    @constraint(
        m,
        [i = 1:length(t_PLR)-1, k = 1:numHours], 
        q_n[i, k] ≤ S_n[i, k] * t_PLR[i + 1] * Q_avail(T_cw_norm, T_out[k])
    );
    @constraint(
        m, 
        conPowConsumption[k = 1:numHours], 
        p[k] == sum(a_Pch_p[:,k] .* q_p[:,k] + b_Pch_p[:,k] .* S_p[:,k]) 
                + sum(a_Pch_n[:,k] .* q_n[:,k] + b_Pch_n[:,k] .* S_n[:,k]) 
                - (1 - b[k]) * P_chiller(Q_cool[k], T_cw_ch, T_out[k]) 
                - b[k] * P_chiller(Q_cool[k], T_cw_norm, T_out[k])
    );

    @constraint(m, conReserveThermalCharging[k = 1:numHours], u_r[k] == u_p_r[k] + u_n_r[k]);
    @constraint(m, conReserveThermalChargingUB[k = 1:numHours], u_p_r[k] ≤ u_p[k]); 
    @constraint(m, [k = 1:numHours], u_p_r[k] ≤ b_r[k] * 1e10); 

    # @constraint(m, [k = 1:numHours], l[k] == sum(l_LB_r[:, k]));
    # @constraint(m, [k = 1:numHours], 1 == sum(S_LB_r[:, k]));
    # @constraint(m, [i = 1:length(t_LB)-1, k = 1:numHours], S_LB_r[i, k] * t_LB[i] ≤ l_LB_r[i, k]);
    # @constraint(m, [i = 1:length(t_LB)-1, k = 1:numHours], l_LB_r[i, k] ≤ S_LB_r[i, k] * t_LB[i + 1]);
    @constraint(m, conReserveThermalChargingLB[k = 1:numHours], u_n_r[k] ≥ -sum(a_LB .* l_LB[:, k] + b_LB .* S_LB[:, k])); # u_LB(l[k])
    @constraint(m, [k = 1:numHours], u_n_r[k] == -(1 - b_r[k]) * Q_cool[k]; # Constraining the reserve discharging rate u_n_r[k] to be either 0 or -Q_cool[k]
    # @constraint(m, [k = 1:numHours], u_n_r[k] ≥ (1 - b_r[k]) * -1e10); # u_LB(l[k])

    # @constraint(m, [k = 1:numHours], u_r[k] + Q_cool[k] ≥ 0);

    @constraint(m, [k = 1:numHours], Q_cool[k] + u_p_r[k] == sum(q_p_r[:, k]));
    @constraint(m, [k = 1:numHours], sum(S_p_r[:, k]) == 1);
    @constraint(
        m,
        [i = 1:length(t_PLR)-1, k = 1:numHours], 
        S_p_r[i, k] * t_PLR[i] * Q_avail(T_cw_ch, T_out[k]) ≤ q_p_r[i, k]
    );
    @constraint(
        m,
        [i = 1:length(t_PLR)-1, k = 1:numHours], 
        q_p_r[i, k] ≤ S_p_r[i, k] * t_PLR[i + 1] * Q_avail(T_cw_ch, T_out[k])
    );
    @constraint(m, [k = 1:numHours], Q_cool[k] + u_n_r[k] == sum(q_n_r[:, k]));
    @constraint(m, [k = 1:numHours], sum(S_n_r[:, k]) == 1);
    @constraint(
        m,
        [i = 1:length(t_PLR)-1, k = 1:numHours], 
        S_n_r[i, k] * t_PLR[i] * Q_avail(T_cw_norm, T_out[k]) ≤ q_n_r[i, k]
    );
    @constraint(
        m,
        [i = 1:length(t_PLR)-1, k = 1:numHours], 
        q_n_r[i, k] ≤ S_n_r[i, k] * t_PLR[i + 1] * Q_avail(T_cw_norm, T_out[k])
    );
    @constraint(
        m, 
        conReserveCapacity[k = 1:numHours], 
        r[k] == p[k] 
            - sum(a_Pch_p[:,k] .* q_p_r[:,k] + b_Pch_p[:,k] .* S_p_r[:,k]) 
            - sum(a_Pch_n[:,k] .* q_n_r[:,k] + b_Pch_n[:,k] .* S_n_r[:,k]) 
            + (1 - b_r[k]) * P_chiller(Q_cool[k], T_cw_ch, T_out[k]) 
            + b_r[k] * P_chiller(Q_cool[k], T_cw_norm, T_out[k])
    );

    @constraint(m, conReserveEnergyLimit[k = 1:numHours], l[k] + u_r[k] / Q_stor ≥ 0);

    @objective(m, Max, sum(-energyPrice .* p ./ 1000 + reservePrice .* r ./ 1000));
    optimize!(m);

    if termination_status(m) ≠ MOI.LOCALLY_SOLVED && termination_status(m) ≠ MOI.OPTIMAL
        println(iobuffer, "Optimization solved with status: ", termination_status(m));
        @warn String(take!(iobuffer));
    end

    powerInject = -value.(p) ./ 1000; # MW
    powerReserve = value.(r) ./ 1000; # MW
    soc = value.(l);
    therCh = value.(u);

    # df = DataFrame(
    #     Energy_Price=vcat(missing, energyPrice),
    #     Reserve_Price=vcat(missing, reservePrice), 
    #     T_out=vcat(missing, T_out),
    #     Q_cool=vcat(missing, Q_cool),
    #     SoC=soc, 
    #     thermalCharging=vcat(missing, therCh), 
    #     Power=vcat(missing, -powerInject), 
    #     Reserve=vcat(missing, powerReserve)
    # );
    # CSV.write("TESS_test.csv", df);

    return powerInject, powerReserve, soc, therCh;
end