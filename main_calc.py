from dataclasses import dataclass
import numpy as np

import aero_log_parser
import domain
import calc
from aero_constraints import violations_type4_points
from aero_log_parser import join_logs_to_points
from calc.aerodynamics import prepare_aerodynamics, calculate_weights
from calc.ballistics import calculate_ballistics
from calc.misc import reinstall_express

reinstall_express()

# ============================================================
# INPUT DATA
# ============================================================

m_0 = 705.816
m_fuel = 315
t_engine = 43.864219
m_empty = m_0 - m_fuel

d_middle = 0.4
V_0 = 375
J_1 = 2100
d_nozzle = 0.25
theta_0 = np.deg2rad(12.0)

V_target = 340

Y_min = 900
Y_max = 9800

Y_1 = Y_min
X_1 = Y_1 / np.tan(theta_0)

Y_4 = Y_min
X_4 = Y_4 / np.tan(theta_0)

mdot_A = 13.404203
mdot_B = -0.294651

aero_info = domain.AerodynamicsInfo(
    L_0=5.3351,
    L_head=0.9000,
    L_stern=0.7296,
    d_M=0.4,
    d_stern=0.3,
    L_cm=0,

    # Steering
    L_st_position=0.983100570703258,
    L_st_span=1.3194249254574384,
    L_st=0.0626025063978535,
    L_st_straight=0.07916421500336057,
    delta_st=0.02133966836451949,

    # Wing
    L_w_position=4.651368704733642,
    L_w_span=1.9256349277294273,
    L_w=0.30038120257225964,
    L_w_straight=0.29261611495574114,
    delta_w=0.03092507039939902
)

w_info = domain.WeightInfo(
    m_body=310.8413,
    a_body=665.0554,
    m_fuel=315.0,
    x_cm_fuel=3.7036,
    rho_w=2600,
    rho_st=2600,
)

print(calculate_weights(w_info, aero_info))

prepare_aerodynamics(aero_info, '3.ad')


# ============================================================
# CHECK (DETAILED)
# ============================================================

@dataclass(frozen=True)
class CheckFlags:
    mdot_ok: bool
    theta_ok: bool
    mass_ok: bool
    velocity_ok: bool

    @property
    def ok(self) -> bool:
        return self.mdot_ok and self.theta_ok and self.mass_ok and self.velocity_ok


def check_detailed(r: domain.BallisticsCalculationResult) -> CheckFlags:
    return CheckFlags(
        mdot_ok=(r.mdot_final >= 0),
        theta_ok=(r.theta_final >= 0),
        mass_ok=(r.m_final >= m_empty),
        velocity_ok=(1.5 * V_target <= r.v_final <= 2.5 * V_target),
    )


# ============================================================
# CALC
# ============================================================
# t_engine: float
#     d_nozzle: float
#     J: float
#
#     V_0: float
#     theta_0: float
#
#     p_nozzle: float = 101325

target_1 = calculate_ballistics(
    target=domain.TargetInfo(
        velocity=V_target,
        x=X_1,
        y=Y_1,
    ),
    la=domain.LAInfo(
        m_0=m_0,
        d_middle=d_middle,

        mdot_A=mdot_A,
        mdot_B=mdot_B,

        t_engine=t_engine,
        d_nozzle=d_nozzle,
        J=J_1,
        V_0=V_0,
        theta_0=theta_0,
        aerodynamics='3.ad',
    ),
    integration=domain.IntegrationInfo(),
    calculate_aerodynamics=True,
)

# target_4 = calculate_ballistics(
#     target=domain.TargetInfo(
#         velocity=V_target,
#         x = X_4,
#         y = Y_4,
#     ),
#     la=domain.LAInfo(
#         m_0=m_0,
#         d_middle=d_middle,
#
#         mdot_A=mdot_A,
#         mdot_B=mdot_B,
#
#         t_engine=t_engine,
#         d_nozzle=d_nozzle,
#         J=J_1,
#         V_0=V_0,
#         theta_0=theta_0,
#         aerodynamics='2.ad',
#     ),
#     integration=domain.IntegrationInfo(),
#     calculate_aerodynamics=False,
# )
print("TARGET 1")
print(target_1.X_final)
print(target_1.trajectory_log_position)
print()
print(target_1.trajectory_log_forces)
print()
print(target_1.trajectory_log_additional)
print(check_detailed(target_1))
print()

force = aero_log_parser.parse_forces_log(target_1.trajectory_log_forces)
position = aero_log_parser.parse_position_log(target_1.trajectory_log_position)
additional = aero_log_parser.parse_additional_log(target_1.trajectory_log_additional)
points =  join_logs_to_points(force, position, additional)
const = {
    "m_0inp": 700.0,
    "L_0": 5.3351,
    "L_head": 0.9000,
    "L_stern": 0.7296,
    "L_warhead_start": 1.2,
    "V_target": 340.0,
    "n_ymax": 12.0,

    "j_n_max": 117.72,
    "alpha_max_deg": 13.0,

    "d_M": 0.4,
    "d_stern": 0.3,
    "L_cm": 0.0,

    # scales
    "sL": 0.1,
    "sM": 1.0,
    "sV": 10.0,
    "sTheta": np.deg2rad(1.0),
    "sMdot": 0.1,
    "sAlpha": np.deg2rad(1.0),
    "sNy": 1.0,
    "sCya": 0.1,
    "sR": 0.01,

    # penalties
    "penalty_ballistics_fail": 1e3,
    "penalty_no_points": 1e3,
    "penalty_bad_denom": 1e3,
    "penalty_div0_mz": 1e3,
}

v4 = violations_type4_points(points, aero_info, const)
print(v4)


# print("TARGET 4")
# print(target_4)
# print(check_detailed(target_4))
# print()
