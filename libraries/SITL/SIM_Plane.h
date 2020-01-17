/*
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 /*
   simple plane simulator class
 */

#pragma once

#include "SIM_Aircraft.h"
#include "SIM_ICEngine.h"
#include <Filter/LowPassFilter.h>

namespace SITL {

    /*
      a very simple plane simulator
     */
    class Plane : public Aircraft {
    public:
        Plane(const char *home_str, const char *frame_str);

        /* update model by one time step */
        virtual void update(const struct sitl_input &input);

        /* static object creator */
        static Aircraft *create(const char *home_str, const char *frame_str) {
            return new Plane(home_str, frame_str);
        }

    protected:
        const float hover_throttle = 0.7f;
        const float air_density = 1.225; // kg/m^3 at sea level, ISA conditions
        float angle_of_attack;
        float beta;
        double X; //
        double X_0; // 
        double initilization_distance; //
         //

        // NEW TETHER MODEL
        double X0_max;
        double X0_min;
        double F_tether_max;
        double Elong_tether_max;
        double K_tether;
        double d_tether;
        double massperlength;
        double C_D_tether;
        double X0;
        double F_tether;


        double W_tether;
        double F_W_NED_x;
        double F_W_NED_y;
        double F_W_NED_z;

        double F_W_BODY_x;
        double F_W_BODY_y;
        double F_W_BODY_z;

        double W_e_mod;
        double Velocity_air_bf_unit_x;
        double Velocity_air_bf_unit_y;
        double Velocity_air_bf_unit_z;

        double F_drag_BODY_x;
        double F_drag_BODY_y;
        double F_drag_BODY_z;

        double Tension_2;
        double Tension_1;
        double Speed_1;
        double Speed_2;
        double m_spd_tension;
        double Tension_cntl;
        double rho_tether;





        bool updating_time;
        //AKM: is reeling out?
        bool reeling_out;
        //AKM: define reeling out time start
        bool radius_min_reached;//
        bool radius_max_reached;//
        double time_reelout_start_ms; //
        double time_reelout_elapsed_s; //
        double time_reelin_elapsed_s; //
        double time_reelin_start_ms; //
        //AKM: define reel-out speed
        double reelout_speed; //
        //AKM: reel-out distance
        double reelin_speed; //
        double radius_min; //
        double radius_max; //
        double F_tether_total;
        double F_tether_int;
        double F_tether_min;
        double F_tether_perc;
        double Winch_alpha;
        double Winch_omega;
        double update_time;
        double delta_time_s;
        double Reeling_speed;
        double Winch_radius;
        double Reeling_distance;
        double F_tether_zero;
        double Winch_omega_max;
        // for catenary tether model
        double R_sphere;
        double drag_tether;
        double weight_tether;
        double load_tether;
        double sag_tether;
        double Rx_tether;
        double E_times_area;

        double F_tether_NED_x;
        double F_tether_NED_y;
        double F_tether_NED_z;
        double F_tether_BODY_x;
        double F_tether_BODY_y;
        double F_tether_BODY_z;

        double X_0_min;
        double X_max;
        double C_4;
        double Inertia_motor;// motor inertia
        double Torque_rated; // rated torque
        double b_motor;// motor viscous friction
        double winch_gain; // multiplies omega
        double time_wait;
        double delta_time_delay;
        double m_omega;
        double limit_omega;
        double X_0_prev;
        double L_unit_x;
        double L_unit_y;
        double L_unit_z;
        double F_tether_const;





        struct {
            // from last_letter skywalker_2013/aerodynamics.yaml
            // thanks to Georacer!
            float s = 0.45;
            float b = 1.88;
            float c = 0.24;
            float c_lift_0 = 0.56;
            float c_lift_deltae = 0;
            float c_lift_a = 6.9;
            float c_lift_q = 0;
            float mcoeff = 50;
            float oswald = 0.9;
            float alpha_stall = 0.4712;
            float c_drag_q = 0;
            float c_drag_deltae = 0.0;
            float c_drag_p = 0.1;
            //float c_drag_p = 1.3; // FOr simulating the quadplane drag
            float c_y_0 = 0;
            float c_y_b = -0.98;
            float c_y_p = 0;
            float c_y_r = 0;
            float c_y_deltaa = 0;
            float c_y_deltar = -0.2;
            float c_l_0 = 0;
            float c_l_p = -1.0;
            float c_l_b = -0.12;
            float c_l_r = 0.14;
            float c_l_deltaa = 0.25;
            float c_l_deltar = -0.037;
            float c_m_0 = 0.045;
            float c_m_a = -0.7;
            float c_m_q = -20;
            float c_m_deltae = 1.0;
            float c_n_0 = 0;
            float c_n_b = 0.25;
            float c_n_p = 0.022;
            float c_n_r = -1;
            float c_n_deltaa = 0.00;
            float c_n_deltar = 0.1;
            float deltaa_max = 0.3491;
            float deltae_max = 0.3491;
            float deltar_max = 0.3491;
            // the X CoG offset should be -0.02, but that makes the plane too tail heavy
            // in manual flight. Adjusted to -0.15 gives reasonable flight
            Vector3f CGOffset{ -0.15, 0, -0.05 };
        } coefficient;

        float thrust_scale;
        bool reverse_thrust;
        bool elevons;
        bool vtail;
        bool dspoilers;
        bool reverse_elevator_rudder;
        bool ice_engine;
        bool tailsitter;
        bool have_launcher;
        float launch_accel;
        float launch_time;
        uint64_t launch_start_ms;

        ICEngine icengine{ 2, 14, 12, 13, 100 };

        float liftCoeff(float alpha) const;
        float dragCoeff(float alpha) const;
        Vector3f getForce(float inputAileron, float inputElevator, float inputRudder) const;
        Vector3f getTorque(float inputAileron, float inputElevator, float inputRudder, float inputThrust, const Vector3f &force) const;
        void calculate_forces(const struct sitl_input &input, Vector3f &rot_accel, Vector3f &body_accel);
        // AKM: adding the calculate location function 
            /* update location from position */
        //void update_position(void);
        //double distance_to_home;
    };

} // namespace SITL
