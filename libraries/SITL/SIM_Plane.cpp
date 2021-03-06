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
   very simple plane simulator class. Not aerodynamically accurate,
   just enough to be able to debug control logic for new frame types
 */

#include "SIM_Plane.h"
 // AKM
 // AKM: need to include this to calculate the position... or
 // or.. put the function definition on SIM_Plane.h



#include <stdio.h>

using namespace SITL;



Plane::Plane(const char *home_str, const char *frame_str) :
    Aircraft(home_str, frame_str)
{
    mass = 2.0f;

    /*
       scaling from motor power to Newtons. Allows the plane to hold
       vertically against gravity when the motor is at hover_throttle
    */
    thrust_scale = (mass * GRAVITY_MSS) / hover_throttle;
    frame_height = 0.1f;

    ground_behavior = GROUND_BEHAVIOR_FWD_ONLY;

    if (strstr(frame_str, "-heavy")) {
        mass = 8;
    }
    if (strstr(frame_str, "-revthrust")) {
        reverse_thrust = true;
    }
    if (strstr(frame_str, "-elevon")) {
        elevons = true;
    }
    else if (strstr(frame_str, "-vtail")) {
        vtail = true;
    }
    else if (strstr(frame_str, "-dspoilers")) {
        dspoilers = true;
    }
    if (strstr(frame_str, "-elevrev")) {
        reverse_elevator_rudder = true;
    }
    if (strstr(frame_str, "-catapult")) {
        have_launcher = true;
        launch_accel = 15;
        launch_time = 2;
    }
    if (strstr(frame_str, "-bungee")) {
        have_launcher = true;
        launch_accel = 7;
        launch_time = 4;
    }
    if (strstr(frame_str, "-throw")) {
        have_launcher = true;
        launch_accel = 10;
        launch_time = 1;
    }
    if (strstr(frame_str, "-tailsitter")) {
        tailsitter = true;
        ground_behavior = GROUND_BEHAVIOR_TAILSITTER;
        thrust_scale *= 1.5;
    }

    if (strstr(frame_str, "-ice")) {
        ice_engine = true;
    }
}

/*
  the following functions are from last_letter
  https://github.com/Georacer/last_letter/blob/master/last_letter/src/aerodynamicsLib.cpp
  many thanks to Georacer!
 */
float Plane::liftCoeff(float alpha) const
{
    const float alpha0 = coefficient.alpha_stall;
    const float M = coefficient.mcoeff;
    const float c_lift_0 = coefficient.c_lift_0;
    const float c_lift_a0 = coefficient.c_lift_a;

    // clamp the value of alpha to avoid exp(90) in calculation of sigmoid
    const float max_alpha_delta = 0.8f;
    if (alpha - alpha0 > max_alpha_delta) {
        alpha = alpha0 + max_alpha_delta;
    }
    else if (alpha0 - alpha > max_alpha_delta) {
        alpha = alpha0 - max_alpha_delta;
    }
    double sigmoid = (1 + exp(-M * (alpha - alpha0)) + exp(M*(alpha + alpha0))) / (1 + exp(-M * (alpha - alpha0))) / (1 + exp(M*(alpha + alpha0)));
    double linear = (1.0 - sigmoid) * (c_lift_0 + c_lift_a0 * alpha); //Lift at small AoA
    double flatPlate = sigmoid * (2 * copysign(1, alpha)*pow(sin(alpha), 2)*cos(alpha)); //Lift beyond stall

    float result = linear + flatPlate;
    return result;
}

float Plane::dragCoeff(float alpha) const
{
    const float b = coefficient.b;
    const float s = coefficient.s;
    const float c_drag_p = coefficient.c_drag_p;
    const float c_lift_0 = coefficient.c_lift_0;
    const float c_lift_a0 = coefficient.c_lift_a;
    const float oswald = coefficient.oswald;

    double AR = pow(b, 2) / s;
    double c_drag_a = c_drag_p + pow(c_lift_0 + c_lift_a0 * alpha, 2) / (M_PI*oswald*AR);

    return c_drag_a;
}

// Torque calculation function
Vector3f Plane::getTorque(float inputAileron, float inputElevator, float inputRudder, float inputThrust, const Vector3f &force) const
{
    float alpha = angle_of_attack;

    //calculate aerodynamic torque
    float effective_airspeed = airspeed;

    if (tailsitter) {
        /*
          tailsitters get airspeed from prop-wash
         */
        effective_airspeed += inputThrust * 20;

        // reduce effective angle of attack as thrust increases
        alpha *= constrain_float(1 - inputThrust, 0, 1);
    }

    const float s = coefficient.s;
    const float c = coefficient.c;
    const float b = coefficient.b;
    const float c_l_0 = coefficient.c_l_0;
    const float c_l_b = coefficient.c_l_b;
    const float c_l_p = coefficient.c_l_p;
    const float c_l_r = coefficient.c_l_r;
    const float c_l_deltaa = coefficient.c_l_deltaa;
    const float c_l_deltar = coefficient.c_l_deltar;
    const float c_m_0 = coefficient.c_m_0;
    const float c_m_a = coefficient.c_m_a;
    const float c_m_q = coefficient.c_m_q;
    const float c_m_deltae = coefficient.c_m_deltae;
    const float c_n_0 = coefficient.c_n_0;
    const float c_n_b = coefficient.c_n_b;
    const float c_n_p = coefficient.c_n_p;
    const float c_n_r = coefficient.c_n_r;
    const float c_n_deltaa = coefficient.c_n_deltaa;
    const float c_n_deltar = coefficient.c_n_deltar;
    const Vector3f &CGOffset = coefficient.CGOffset;

    float rho = air_density;

    //read angular rates
    double p = gyro.x;
    double q = gyro.y;
    double r = gyro.z;

    double qbar = 1.0 / 2.0*rho*pow(effective_airspeed, 2)*s; //Calculate dynamic pressure
    double la, na, ma;
    if (is_zero(effective_airspeed))
    {
        la = 0;
        ma = 0;
        na = 0;
    }
    else
    {
        la = qbar * b*(c_l_0 + c_l_b * beta + c_l_p * b*p / (2 * effective_airspeed) + c_l_r * b*r / (2 * effective_airspeed) + c_l_deltaa * inputAileron + c_l_deltar * inputRudder); // rolling moment
        ma = qbar * c*(c_m_0 + c_m_a * alpha + c_m_q * c*q / (2 * effective_airspeed) + c_m_deltae * inputElevator);                                                           // pitching moment
        na = qbar * b*(c_n_0 + c_n_b * beta + c_n_p * b*p / (2 * effective_airspeed) + c_n_r * b*r / (2 * effective_airspeed) + c_n_deltaa * inputAileron + c_n_deltar * inputRudder); // yawing moment
    }


    // Add torque to to force misalignment with CG
    // r x F, where r is the distance from CoG to CoL
    la += CGOffset.y * force.z - CGOffset.z * force.y;
    ma += -CGOffset.x * force.z + CGOffset.z * force.x;
    na += -CGOffset.y * force.x + CGOffset.x * force.y;

    return Vector3f(la, ma, na);
}

// AKM: obtaining the location of the plane
/*
void SITL::Plane::update_position(void)
{
    location = home;
    location_offset(location, position.x, position.y);

    location.alt = static_cast<int32_t>(home.alt - position.z * 100.0f);
    double distance_to_home = pow(pow(position.x, 2) + pow(position.y, 2) + pow(position.z, 2), 0.5);

}
*/

//hal.console->printf("hola esto es una prueba: ");



// Force calculation function from last_letter
Vector3f Plane::getForce(float inputAileron, float inputElevator, float inputRudder) const
{
    const float alpha = angle_of_attack;
    const float c_drag_q = coefficient.c_drag_q;
    const float c_lift_q = coefficient.c_lift_q;
    const float s = coefficient.s;
    const float c = coefficient.c;
    const float b = coefficient.b;
    const float c_drag_deltae = coefficient.c_drag_deltae;
    const float c_lift_deltae = coefficient.c_lift_deltae;
    const float c_y_0 = coefficient.c_y_0;
    const float c_y_b = coefficient.c_y_b;
    const float c_y_p = coefficient.c_y_p;
    const float c_y_r = coefficient.c_y_r;
    const float c_y_deltaa = coefficient.c_y_deltaa;
    const float c_y_deltar = coefficient.c_y_deltar;

    float rho = air_density;

    //request lift and drag alpha-coefficients from the corresponding functions
    double c_lift_a = liftCoeff(alpha);
    double c_drag_a = dragCoeff(alpha);

    //convert coefficients to the body frame
    double c_x_a = -c_drag_a * cos(alpha) + c_lift_a * sin(alpha);
    double c_x_q = -c_drag_q * cos(alpha) + c_lift_q * sin(alpha);
    double c_z_a = -c_drag_a * sin(alpha) - c_lift_a * cos(alpha);
    double c_z_q = -c_drag_q * sin(alpha) - c_lift_q * cos(alpha);

    //read angular rates
    double p = gyro.x;
    double q = gyro.y;
    double r = gyro.z;

    // AKM
    // calculate distance to home here!!!
    //double distance_to_home = pow(pow(position.x, 2) + pow(position.y, 2) + pow(position.z, 2), 0.5);
    //hal.console->printf("hola esto es una prueba: ");
    //print_distance_to_home(hal.console, distance_to_home);
    // with pre-calculated values
    //double distance_to_home = pow(pow(location_offset, 2) + pow(location.alt, 2), 0.5);



    //calculate aerodynamic force
    double qbar = 1.0 / 2.0*rho*pow(airspeed, 2)*s; //Calculate dynamic pressure
    double ax, ay, az;
    if (is_zero(airspeed))
    {
        ax = 0;
        ay = 0;
        az = 0;
    }
    else
    {
        ax = qbar * (c_x_a + c_x_q * c*q / (2 * airspeed) - c_drag_deltae * cos(alpha)*fabs(inputElevator) + c_lift_deltae * sin(alpha)*inputElevator);
        // split c_x_deltae to include "abs" term
        ay = qbar * (c_y_0 + c_y_b * beta + c_y_p * b*p / (2 * airspeed) + c_y_r * b*r / (2 * airspeed) + c_y_deltaa * inputAileron + c_y_deltar * inputRudder);
        // AKM: adding distance_to_home as force to see how it bahaves
        // the force must increase as the distance gets bigger.
        az = qbar * (c_z_a + c_z_q * c*q / (2 * airspeed) - c_drag_deltae * sin(alpha)*fabs(inputElevator) - c_lift_deltae * cos(alpha)*inputElevator);// + distance_to_home;
        // split c_z_deltae to include "abs" term
    }
    return Vector3f(ax, ay, az);
}

void Plane::calculate_forces(const struct sitl_input &input, Vector3f &rot_accel, Vector3f &body_accel)
{
    float aileron = filtered_servo_angle(input, 0);
    float elevator = filtered_servo_angle(input, 1);
    float rudder = filtered_servo_angle(input, 3);
    bool launch_triggered = input.servos[6] > 1700;
    float throttle;
    float rho = air_density;
    if (reverse_elevator_rudder) {
        elevator = -elevator;
        rudder = -rudder;
    }
    if (elevons) {
        // fake an elevon plane
        float ch1 = aileron;
        float ch2 = elevator;
        aileron = (ch2 - ch1) / 2.0f;
        // the minus does away with the need for RC2_REV=-1
        elevator = -(ch2 + ch1) / 2.0f;

        // assume no rudder
        rudder = 0;
    }
    else if (vtail) {
        // fake a vtail plane
        float ch1 = elevator;
        float ch2 = rudder;
        // this matches VTAIL_OUTPUT==2
        elevator = (ch2 - ch1) / 2.0f;
        rudder = (ch2 + ch1) / 2.0f;
    }
    else if (dspoilers) {
        // fake a differential spoiler plane. Use outputs 1, 2, 4 and 5
        float dspoiler1_left = filtered_servo_angle(input, 0);
        float dspoiler1_right = filtered_servo_angle(input, 1);
        float dspoiler2_left = filtered_servo_angle(input, 3);
        float dspoiler2_right = filtered_servo_angle(input, 4);
        float elevon_left = (dspoiler1_left + dspoiler2_left) / 2;
        float elevon_right = (dspoiler1_right + dspoiler2_right) / 2;
        aileron = (elevon_right - elevon_left) / 2;
        elevator = (elevon_left + elevon_right) / 2;
        rudder = fabsf(dspoiler1_right - dspoiler2_right) / 2 - fabsf(dspoiler1_left - dspoiler2_left) / 2;
    }
    printf("Aileron: %.1f elevator: %.1f rudder: %.1f\n", aileron, elevator, rudder);

    if (reverse_thrust) {
        throttle = filtered_servo_angle(input, 2);
    }
    else {
        throttle = filtered_servo_range(input, 2);
    }

    float thrust = throttle;

    if (ice_engine) {
        thrust = icengine.update(input);
    }

    // calculate angle of attack
    angle_of_attack = atan2f(velocity_air_bf.z, velocity_air_bf.x);
    beta = atan2f(velocity_air_bf.y, velocity_air_bf.x);

    if (tailsitter) {
        /*
          tailsitters get 4x the control surfaces
         */
        aileron *= 4;
        elevator *= 4;
        rudder *= 4;
    }

    Vector3f force = getForce(aileron, elevator, rudder);
    // =======================================================================================================================
    // AKM: add force here in z direction (must be in the home direction!)
    // create if loop that activates tether force if distance_to_home > R_sphere.
        // =======================================================================================================================

    // extract the Euler angles
    float r, p, y;
    dcm.to_euler(&r, &p, &y);


    // NEW TETHER MODEL USING WEIGHT, SPRING TENSION, AND DRAG
    // =======================================================


    // need to have an initialization distance to avoid singularities at We = 0
    //define tether specifications
    
    X0_min = 30; //m
    F_tether_max = 5800; //N
    Elong_tether_max = 0.034; // 1/100
    d_tether = 0.0016; //m
    rho_tether = 970;
    massperlength = rho_tether * 3.1416*pow(d_tether, 2)/4; //0.001; // kg/m
    C_D_tether = 1;

    X = pow(pow(position.x, 2) + pow(position.y, 2) + pow(position.z, 2), 0.5); // m
    
/*    // USE CONSTANT REEL SPEED TO CHANGE THE LENGTH OF THE TETHER (X0) THE SPRING MODEL IS THEN APPLIED
         // CONSTANT REEL SPEED MODEL
        // ^^^^^^^^^^^^^^^^^^^^^^
        // initialization distance, replaces the control_mode set initializer
    initilization_distance = 30.0; //meters
    // define the minimun and maximum radius
    radius_min = 100; // meters, make sure this is the same as Navigation.cpp

    /// CHECK THIS OUT! THE RADIUS MIN IS THE X_0 RIGHT AFTER THE INITIALIZATION DISTANCE IS REACHED, RESULTING IN dELTA X VERY BIG!! THEREFORE HIGH LOADS AT X=20

    radius_max = 200; // meters.
    // define the reeling speeds (constant)
    reelout_speed = 1.0; // m/s
    reelin_speed = 5.0; // m/s
    X = pow(pow(position.x, 2) + pow(position.y, 2) + pow(position.z, 2), 0.5); // in meters
    // set initial timing
    updating_time = true;
    // if the plane mode is farther than defined distance, then we are not updating time
    if (X >= initilization_distance) {
        updating_time = false;
    }
    //the timer stops updating when X pass the initilization distance,
    if (updating_time) {
        time_reelout_start_ms = AP_HAL::millis();
        time_reelin_start_ms = AP_HAL::millis();
        X0_max = X;  //winch is reeling out with the plane until reaching the initialization distance
    }
    // change position verification from X to X_0
   // Where are we?
   // 1) Below R_min
    if (X0_max < radius_min) {
        radius_min_reached = true;
        radius_max_reached = false;
    }
    // 2) Above R_max
    else if (X0_max >= radius_max) {
        radius_min_reached = false;
        radius_max_reached = true;
    }
    // which direction are we going?
    if (radius_min_reached == true && radius_max_reached == false) {
        // we are reeling out
        time_reelin_start_ms = AP_HAL::millis();
        time_reelout_elapsed_s = (AP_HAL::millis() - time_reelout_start_ms) * 0.001;
        X0_max = radius_min + (reelout_speed * time_reelout_elapsed_s); // in meters!
    }
    else if (radius_min_reached == false && radius_max_reached == true) {
        // we are reeling in
        time_reelout_start_ms = AP_HAL::millis();
        time_reelin_elapsed_s = (AP_HAL::millis() - time_reelin_start_ms) * 0.001;
        X0_max = radius_max - (reelin_speed * time_reelin_elapsed_s); // in meters!
    }

    // END OF CONSTANT SPEED MODEL
    */

    //COMMENT HERE FOR REMOVING THE TETHER MODEL (LATEST AND WORKING)


    // OR NOT APPLY TETHER TENSION
    X0 = X + 1;
    X0_max = 1000; //m

    // no tension, just weight and drag
    if ((X <= X0) && (X >= X0_min)) {

        // CALCULATE THE REELING SPEED
        // calculate the direction of the tether
        L_unit_x = position.x / X;
        L_unit_y = position.y / X;
        L_unit_z = -position.z / X;
        // Reeling Velocity: Dot product the Velocity in XYZ with the tether direction 
        // velocity_ef is the velocity in the earth frame. Z is pointing down (switch)
        Reeling_speed = velocity_ef.x * L_unit_x + velocity_ef.y * L_unit_y + -velocity_ef.z * L_unit_z;

        // PROPORTIONAL TENSION CONTROLLER AS A FUNCTION OF REELING SPEED
      //  Tension_2 = 30;
      //  Tension_1 = 0;
      //  Speed_2 = 10;
      //  Speed_1 = -10;
      //  m_spd_tension = (Tension_2 - Tension_1) / (Speed_2 - Speed_1);
      //  Tension_cntl = m_spd_tension * Reeling_speed;
       // if (Tension_cntl < 0) {
       //     Tension_cntl = 0;
       // }
       



        // TWO TENSION SET-POINTS FOR TRACTION AND RETRACTION
        //if (Reeling_speed < 0) {
        //    Tension_cntl = 15;
        //}
        //else {
        //    Tension_cntl = 15;
        //}
        
        Tension_cntl = 15;

        // TETHER WEIGHT
        // -------------
        W_tether = massperlength * X0 * 9.81;
        // in the Z direction of the NED coordinate system
        F_W_NED_z = W_tether;
        F_W_NED_x = 0.0;
        F_W_NED_y = 0.0;

        // rotate NED to BODY
        F_W_BODY_x = cos(y)*cos(p)*F_W_NED_x + cos(p)*sin(y)*F_W_NED_y - sin(p)*F_W_NED_z;
        F_W_BODY_y = (cos(y)*sin(r)*sin(p) - cos(r)*sin(y))*F_W_NED_x + (cos(r)*cos(y) + sin(r)*sin(y)*sin(p))*F_W_NED_y + cos(p)*sin(r)*F_W_NED_z;
        F_W_BODY_z = (sin(y)*sin(r) + cos(r)*cos(y)*sin(p))*F_W_NED_x + (cos(r)*sin(y)*sin(p) - cos(y)*sin(r))*F_W_NED_y + cos(r)*cos(p)*F_W_NED_z;

        // TETHER DRAG
        // -----------
        // calculate the module of the apparent velocity
        W_e_mod = pow(pow(velocity_air_bf.x, 2) + pow(velocity_air_bf.y, 2) + pow(velocity_air_bf.x, 2), 0.5);

        // calculate the apparent velocity unit vectors
        Velocity_air_bf_unit_x = velocity_air_bf.x / W_e_mod;
        Velocity_air_bf_unit_y = velocity_air_bf.y / W_e_mod;
        Velocity_air_bf_unit_z = velocity_air_bf.z / W_e_mod;

        // caclulate the drag force in BODY axes
        F_drag_BODY_x = rho * C_D_tether * X0*d_tether * cos(angle_of_attack) * pow(W_e_mod, 2) * Velocity_air_bf_unit_x / 8;
        F_drag_BODY_y = rho * C_D_tether * X0*d_tether * cos(angle_of_attack) * pow(W_e_mod, 2) * Velocity_air_bf_unit_y / 8;
        F_drag_BODY_z = rho * C_D_tether * X0*d_tether * cos(angle_of_attack) * pow(W_e_mod, 2) * Velocity_air_bf_unit_z / 8;


        // TETHER TENSION: CONSTANT
        // ------------------------
        F_tether_NED_x = -Tension_cntl * (position.x / X);
        F_tether_NED_y = -Tension_cntl * (position.y / X);
        F_tether_NED_z = -Tension_cntl * (position.z / X);
        // rotate the forces to body axes
        F_tether_BODY_x = cos(y)*cos(p)*F_tether_NED_x + cos(p)*sin(y)*F_tether_NED_y - sin(p)*F_tether_NED_z;
        F_tether_BODY_y = (cos(y)*sin(r)*sin(p) - cos(r)*sin(y))*F_tether_NED_x + (cos(r)*cos(y) + sin(r)*sin(y)*sin(p))*F_tether_NED_y + cos(p)*sin(r)*F_tether_NED_z;
        F_tether_BODY_z = (sin(y)*sin(r) + cos(r)*cos(y)*sin(p))*F_tether_NED_x + (cos(r)*sin(y)*sin(p) - cos(y)*sin(r))*F_tether_NED_y + cos(r)*cos(p)*F_tether_NED_z;


        // ADD FORCES TO THE BODY FORCES
        force.x = force.x + F_W_BODY_x - F_drag_BODY_x + F_tether_BODY_x;
        force.y = force.y + F_W_BODY_y - F_drag_BODY_y + F_tether_BODY_y;
        force.z = force.z + F_W_BODY_z - F_drag_BODY_z + F_tether_BODY_z;

    }

    // Tension applied, plus weight and drag
    if ((X >= X0_max) && (X >= X0_min)) {

        // CALCULATE THE REELING SPEED
        // calculate the direction of the tether
        L_unit_x = position.x / X;
        L_unit_y = position.y / X;
        L_unit_z = -position.z / X;
        // Reeling Velocity: Dot product the Velocity in XYZ with the tether direction 
        // velocity_ef is the velocity in the earth frame. Z is pointing down (switch)
        Reeling_speed = velocity_ef.x * L_unit_x + velocity_ef.y * L_unit_y + -velocity_ef.z * L_unit_z;


        // TENSION FORCE
        // -------------
        // calculate tether elasticity constant
        K_tether = F_tether_max / (Elong_tether_max * X0);

        // calculate total tether tension force 
        F_tether = K_tether * (X - X0_max);

        // calculate the NED forces by decomposing the tether tension
        F_tether_NED_x = F_tether * (position.x / X);
        F_tether_NED_y = F_tether * (position.y / X);
        F_tether_NED_z = F_tether * (position.z / X);

        // rotate tether forces to BODY axes
        F_tether_BODY_x = cos(y)*cos(p)*F_tether_NED_x + cos(p)*sin(y)*F_tether_NED_y - sin(p)*F_tether_NED_z;
        F_tether_BODY_y = (cos(y)*sin(r)*sin(p) - cos(r)*sin(y))*F_tether_NED_x + (cos(r)*cos(y) + sin(r)*sin(y)*sin(p))*F_tether_NED_y + cos(p)*sin(r)*F_tether_NED_z;
        F_tether_BODY_z = (sin(y)*sin(r) + cos(r)*cos(y)*sin(p))*F_tether_NED_x + (cos(r)*sin(y)*sin(p) - cos(y)*sin(r))*F_tether_NED_y + cos(r)*cos(p)*F_tether_NED_z;

        // TETHER WEIGHT
        // -------------
        W_tether = massperlength * X0 * 9.81;
        // in the Z direction of the NED coordinate system
        F_W_NED_z = W_tether;
        F_W_NED_x = 0.0;
        F_W_NED_y = 0.0;

        // rotate NED to BODY
        F_W_BODY_x = cos(y)*cos(p)*F_W_NED_x + cos(p)*sin(y)*F_W_NED_y - sin(p)*F_W_NED_z;
        F_W_BODY_y = (cos(y)*sin(r)*sin(p) - cos(r)*sin(y))*F_W_NED_x + (cos(r)*cos(y) + sin(r)*sin(y)*sin(p))*F_W_NED_y + cos(p)*sin(r)*F_W_NED_z;
        F_W_BODY_z = (sin(y)*sin(r) + cos(r)*cos(y)*sin(p))*F_W_NED_x + (cos(r)*sin(y)*sin(p) - cos(y)*sin(r))*F_W_NED_y + cos(r)*cos(p)*F_W_NED_z;

        // TETHER DRAG
        // -----------
        // calculate the module of the apparent velocity
        W_e_mod = pow(pow(velocity_air_bf.x, 2) + pow(velocity_air_bf.y, 2) + pow(velocity_air_bf.x, 2), 0.5);

        // calculate the apparent velocity unit vectors
        Velocity_air_bf_unit_x = velocity_air_bf.x / W_e_mod;
        Velocity_air_bf_unit_y = velocity_air_bf.y / W_e_mod;
        Velocity_air_bf_unit_z = velocity_air_bf.z / W_e_mod;

        // caclulate the drag force in BODY axes
        F_drag_BODY_x = rho * C_D_tether * X0*d_tether * cos(angle_of_attack) * pow(W_e_mod, 2) * Velocity_air_bf_unit_x / 8;
        F_drag_BODY_y = rho * C_D_tether * X0*d_tether * cos(angle_of_attack) * pow(W_e_mod, 2) * Velocity_air_bf_unit_y / 8;
        F_drag_BODY_z = rho * C_D_tether * X0*d_tether * cos(angle_of_attack) * pow(W_e_mod, 2) * Velocity_air_bf_unit_z / 8;


        // ADD FORCES TO THE BODY FORCES
        force.x = force.x - F_W_BODY_x - F_drag_BODY_x - F_tether_BODY_x;
        force.y = force.y - F_W_BODY_y - F_drag_BODY_y - F_tether_BODY_y;
        force.z = force.z - F_W_BODY_z - F_drag_BODY_z - F_tether_BODY_z;

    }

    // END OF NEW TETHER MODEL
    



    /*
     
     // ADD CONSTANT TETHER FORCE AFTER INITIALIZATION DISTANCE:
     X = pow(pow(position.x, 2) + pow(position.y, 2) + pow(position.z, 2), 0.5); // in meters
     X_0 = 60;
     X_max = 150;
     // K_tether = 950;
     F_tether_max = 4000; // [N]
     Elong_tether_max = 0.034; // [ %/100 ]

     K_tether = F_tether_max / (Elong_tether_max * X);

     if (X > X_0 && X < X_max){
     // CONSTANT FORCE for X_0_ref calculation.
     F_tether_NED_x = 0.0 * (position.x / X);
     F_tether_NED_y = 0.0 * (position.y / X);
     F_tether_NED_z = 0.0 * (position.z / X);
     // rotate the forces to body axes
     F_tether_BODY_x = cos(y)*cos(p)*F_tether_NED_x + cos(p)*sin(y)*F_tether_NED_y - sin(p)*F_tether_NED_z;
     F_tether_BODY_y = (cos(y)*sin(r)*sin(p) - cos(r)*sin(y))*F_tether_NED_x + (cos(r)*cos(y) + sin(r)*sin(y)*sin(p))*F_tether_NED_y + cos(p)*sin(r)*F_tether_NED_z;
     F_tether_BODY_z = (sin(y)*sin(r) + cos(r)*cos(y)*sin(p))*F_tether_NED_x + (cos(r)*sin(y)*sin(p) - cos(y)*sin(r))*F_tether_NED_y + cos(r)*cos(p)*F_tether_NED_z;
     }



     
     // add spring model when the max tether length is reached

     if (X > X_max) {
     // spring model after reaching X_max
     F_tether_NED_x = K_tether * (X - X_max) * (position.x / X);
     F_tether_NED_y = K_tether * (X - X_max) * (position.y / X);
     F_tether_NED_z = K_tether * (X - X_max) * (position.z / X);
     // rotate the forces to body axes
     F_tether_BODY_x = cos(y)*cos(p)*F_tether_NED_x + cos(p)*sin(y)*F_tether_NED_y - sin(p)*F_tether_NED_z;
     F_tether_BODY_y = (cos(y)*sin(r)*sin(p) - cos(r)*sin(y))*F_tether_NED_x + (cos(r)*cos(y) + sin(r)*sin(y)*sin(p))*F_tether_NED_y + cos(p)*sin(r)*F_tether_NED_z;
     F_tether_BODY_z = (sin(y)*sin(r) + cos(r)*cos(y)*sin(p))*F_tether_NED_x + (cos(r)*sin(y)*sin(p) - cos(y)*sin(r))*F_tether_NED_y + cos(r)*cos(p)*F_tether_NED_z;
     }

     // add the tether force to the plane forces
     force.x = force.x - F_tether_BODY_x;
     force.y = force.y - F_tether_BODY_y;
     force.z = force.z - F_tether_BODY_z;

     */




     /*   // =======================================================================================================================
    // AKM: create the variable tether length constraint (look in Navigation.cpp)

         // OPTION 1: PUMPING MODE, CONSTANT REEL SPEED
        // ^^^^^^^^^^^^^^^^^^^^^^
        // initialization distance, replaces the control_mode set initializer
        initilization_distance = 20.0; //meters
        // define the minimun and maximum radius
        radius_min = 50; // meters, make sure this is the same as Navigation.cpp

        /// CHECK THIS OUT! THE RADIUS MIN IS THE X_0 RIGHT AFTER THE INITIALIZATION DISTANCE IS REACHED, RESULTING IN dELTA X VERY BIG!! THEREFORE HIGH LOADS AT X=20

        radius_max = 300; // meters.
        // define the reeling speeds (constant)
        reelout_speed = 1.0; // m/s
        reelin_speed = 1.0; // m/s
        X = pow(pow(position.x, 2) + pow(position.y, 2) + pow(position.z, 2), 0.5); // in meters
        // set initial timing
        updating_time = true;
        // if the plane mode is farther than defined distance, then we are not updating time
        if (X >= initilization_distance) {
            updating_time = false;
        }
        //the timer stops updating when X pass the initilization distance,
        if (updating_time) {
            time_reelout_start_ms = AP_HAL::millis();
            time_reelin_start_ms = AP_HAL::millis();
            X_0 = X;  //winch is reeling out with the plane until reaching the initialization distance
        }
        // change position verification from X to X_0
       // Where are we?
       // 1) Below R_min
        if (X_0 < radius_min) {
            radius_min_reached = true;
            radius_max_reached = false;
        }
        // 2) Above R_max
        else if (X_0 >= radius_max) {
            radius_min_reached = false;
            radius_max_reached = true;
        }
        // which direction are we going?
        if (radius_min_reached == true && radius_max_reached == false) {
            // we are reeling out
            time_reelin_start_ms = AP_HAL::millis();
            time_reelout_elapsed_s = (AP_HAL::millis() - time_reelout_start_ms) * 0.001;
            X_0 = radius_min + (reelout_speed * time_reelout_elapsed_s); // in meters!
        }
        else if (radius_min_reached == false && radius_max_reached == true) {
            // we are reeling in
            time_reelout_start_ms = AP_HAL::millis();
            time_reelin_elapsed_s = (AP_HAL::millis() - time_reelin_start_ms) * 0.001;
            X_0 = radius_max - (reelin_speed * time_reelin_elapsed_s); // in meters!
        }

        if (X > initilization_distance && X < X_0) { // cable model
            // Tether inclined catenary model
            drag_tether = 0.5 * 1.225 * (3.1416 * pow(0.0016, 2) / 4) * 0.47 * pow(airspeed / 2, 2);
            weight_tether = 0.00194 * 9.8 * X;
            load_tether = drag_tether + (weight_tether * 0.7071);
            // only if the plane X is smaller than R_sphere, add if condition.
            sag_tether = pow((3 * X*X_0 - 3 * pow(X, 2)) / 8, 0.5);
            Rx_tether = load_tether * pow(X, 2) / (8 * sag_tether);
            // calculate the forces in NED coordinates
            F_tether_NED_x = Rx_tether * (position.x / X);
            F_tether_NED_y = Rx_tether * (position.y / X);
            F_tether_NED_z = Rx_tether * (position.z / X);
            // rotate the forces to body axes
            F_tether_BODY_x = cos(y)*cos(p)*F_tether_NED_x + cos(p)*sin(y)*F_tether_NED_y - sin(p)*F_tether_NED_z;
            F_tether_BODY_y = (cos(y)*sin(r)*sin(p) - cos(r)*sin(y))*F_tether_NED_x + (cos(r)*cos(y) + sin(r)*sin(y)*sin(p))*F_tether_NED_y + cos(p)*sin(r)*F_tether_NED_z;
            F_tether_BODY_z = (sin(y)*sin(r) + cos(r)*cos(y)*sin(p))*F_tether_NED_x + (cos(r)*sin(y)*sin(p) - cos(y)*sin(r))*F_tether_NED_y + cos(r)*cos(p)*F_tether_NED_z;
        }

        else if (X > X_0) { // spring model
            E_times_area = 2; //233230; // Modulus times the area using diameter 1.6mm
            K_tether = E_times_area / X_0; // function of the distance to home (X)
            F_tether_NED_x = K_tether * (X - X_0) * (position.x / X);
            F_tether_NED_y = K_tether * (X - X_0) * (position.y / X);
            F_tether_NED_z = K_tether * (X - X_0) * (position.z / X);
            // rotate the forces to body axes
            F_tether_BODY_x = cos(y)*cos(p)*F_tether_NED_x + cos(p)*sin(y)*F_tether_NED_y - sin(p)*F_tether_NED_z;
            F_tether_BODY_y = (cos(y)*sin(r)*sin(p) - cos(r)*sin(y))*F_tether_NED_x + (cos(r)*cos(y) + sin(r)*sin(y)*sin(p))*F_tether_NED_y + cos(p)*sin(r)*F_tether_NED_z;
            F_tether_BODY_z = (sin(y)*sin(r) + cos(r)*cos(y)*sin(p))*F_tether_NED_x + (cos(r)*sin(y)*sin(p) - cos(y)*sin(r))*F_tether_NED_y + cos(r)*cos(p)*F_tether_NED_z;
        }
        else {
            F_tether_BODY_x = 0;
            F_tether_BODY_y = 0;
            F_tether_BODY_z = 0;
        }


        force.x = force.x - F_tether_BODY_x;
        force.y = force.y - F_tether_BODY_y;
        force.z = force.z - F_tether_BODY_z;


        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         // COMMENT OUT HERE THE PUMPING MODE
     */


     //=========================================================================================================================================
 /*   // OPTION 2: FIXED TETHER LENGTH
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //R_sphere = 150; //m same as in Commands_logic.cpp
    //X = pow(pow(position.x, 2) + pow(position.y, 2) + pow(position.z, 2), 0.5); // in meters
    //X_0 = 200; // meters, this must be equal to S1_in_S2.S2_radius_cm in Commands_logic.cpp
    //X_0 = 50; // this is the R_sphere minus the oscillation of the plane distance with K=0
    //K_tether = 1; // constant
    //K_tether = 0.0087 * X; // variable K as a function of the distance to home (X)
    // the maximum is K=0.7 @ L=230m and K=0 @ L=150

    R_sphere = 150; //m length of the tether.
    X_0 = 30; // where the catenary model starts to make effect
    X = pow(pow(position.x, 2) + pow(position.y, 2) + pow(position.z, 2), 0.5); // in meters




    if (X > X_0 && X < R_sphere){
        // Tether inclined catenary model
        drag_tether = 0.5 * 1.225 * 0.0016*X_0 * 0.47 * pow(airspeed / 2, 2);
        weight_tether = 0.00194 * 9.8 * X_0;
        load_tether = drag_tether + (weight_tether * pow(pow(position.x,2)+pow(position.y,2),0.5)/X);
        // only if the plane X is smaller than R_sphere, add if condition.
        sag_tether = pow((3 * X*R_sphere - 3 * pow(X, 2)) / 8, 0.5);
        Rx_tether = load_tether * pow(X, 2) / (8 * sag_tether);
        // calculate the forces in NED coordinates
        F_tether_NED_x = Rx_tether * (position.x / X);
        F_tether_NED_y = Rx_tether * (position.y / X);
        F_tether_NED_z = Rx_tether * (position.z / X);
        // rotate the forces to body axes
        F_tether_BODY_x = cos(y)*cos(p)*F_tether_NED_x + cos(p)*sin(y)*F_tether_NED_y - sin(p)*F_tether_NED_z;
        F_tether_BODY_y = (cos(y)*sin(r)*sin(p) - cos(r)*sin(y))*F_tether_NED_x + (cos(r)*cos(y) + sin(r)*sin(y)*sin(p))*F_tether_NED_y + cos(p)*sin(r)*F_tether_NED_z;
        F_tether_BODY_z = (sin(y)*sin(r) + cos(r)*cos(y)*sin(p))*F_tether_NED_x + (cos(r)*sin(y)*sin(p) - cos(y)*sin(r))*F_tether_NED_y + cos(r)*cos(p)*F_tether_NED_z;
    }
    else if (X > R_sphere) {
        E_times_area = 233230; // Modulus times the area using diameter 1.6mm
        K_tether = E_times_area / R_sphere; // function of the distance to home (X)
        F_tether_NED_x = K_tether * (X - R_sphere) * (position.x / X);
        F_tether_NED_y = K_tether * (X - R_sphere) * (position.y / X);
        F_tether_NED_z = K_tether * (X - R_sphere) * (position.z / X);
        // rotate the forces to body axes
        F_tether_BODY_x = cos(y)*cos(p)*F_tether_NED_x + cos(p)*sin(y)*F_tether_NED_y - sin(p)*F_tether_NED_z;
        F_tether_BODY_y = (cos(y)*sin(r)*sin(p) - cos(r)*sin(y))*F_tether_NED_x + (cos(r)*cos(y) + sin(r)*sin(y)*sin(p))*F_tether_NED_y + cos(p)*sin(r)*F_tether_NED_z;
        F_tether_BODY_z = (sin(y)*sin(r) + cos(r)*cos(y)*sin(p))*F_tether_NED_x + (cos(r)*sin(y)*sin(p) - cos(y)*sin(r))*F_tether_NED_y + cos(r)*cos(p)*F_tether_NED_z;
    }


//    if (X > X_0){
//    // CONSTANT FORCE for X_0_ref calculation.
//    F_tether_NED_x = 30.0 * (position.x / X);
//    F_tether_NED_y = 30.0 * (position.y / X);
//    F_tether_NED_z = 30.0 * (position.z / X);
//    // rotate the forces to body axes
//    F_tether_BODY_x = cos(y)*cos(p)*F_tether_NED_x + cos(p)*sin(y)*F_tether_NED_y - sin(p)*F_tether_NED_z;
//    F_tether_BODY_y = (cos(y)*sin(r)*sin(p) - cos(r)*sin(y))*F_tether_NED_x + (cos(r)*cos(y) + sin(r)*sin(y)*sin(p))*F_tether_NED_y + cos(p)*sin(r)*F_tether_NED_z;
//    F_tether_BODY_z = (sin(y)*sin(r) + cos(r)*cos(y)*sin(p))*F_tether_NED_x + (cos(r)*sin(y)*sin(p) - cos(y)*sin(r))*F_tether_NED_y + cos(r)*cos(p)*F_tether_NED_z;
//    }


    force.x = force.x - F_tether_BODY_x;
    force.y = force.y - F_tether_BODY_y;
    force.z = force.z - F_tether_BODY_z;

    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // COMMENT OUT HERE THE FIXED TETHER LENGRTH
     //=========================================================================================================================================




    // how to calculate the tether constant of elasticity

R_sphere = 240; %m
E_tether = 116E9; % Pa
d_tether = 1.6E-3;  % m
A_tether = pi* d_tether^2 /4; %m^2
K_tether = E_tether * A_tether / R_sphere;


    if (X_0 >= radius_min) {
        K_tether = 116000000000 * 0.000002 / X_0;
    }
    else {
        K_tether = 116000000000 * 0.000002 / radius_min;
    }



    K_tether = 0.70;

    // define sphere radius and spring elasticity constant
    //double R_sphere = 50.0; // meters. Read from Commands_Logic.cpp

    //double K_tether = 700.0;
    //double distance_to_home = pow(pow(position.x, 2) + pow(position.y, 2) + pow(position.z, 2), 0.5); // in meters

    // option 1:
    // use tether model F= K*x only when plane distance to home is bigger than defined sphere radius
    // calculate the forces in NED coordinates

    double F_tether_NED_x = K_tether * (X - X_0) * (position.x / X);
    double F_tether_NED_y = K_tether * (X - X_0) *(position.y / X);
    double F_tether_NED_z = K_tether * (X - X_0) *(position.z / X);


    // option 2:
    // Use constant tether force pointing from the plane to the home location.
    //double F_tether_constant = 20.0; // in Newtons. Use this for constant tether force

    //double F_tether_NED_x = F_tether_constant *(position.x / distance_to_home);
    //double F_tether_NED_y = F_tether_constant *(position.y / distance_to_home);
    //double F_tether_NED_z = F_tether_constant *(position.z / distance_to_home);


    // rotate the forces to body axes
    double F_tether_BODY_x = cos(y)*cos(p)*F_tether_NED_x + cos(p)*sin(y)*F_tether_NED_y - sin(p)*F_tether_NED_z;
    double F_tether_BODY_y = (cos(y)*sin(r)*sin(p) - cos(r)*sin(y))*F_tether_NED_x + (cos(r)*cos(y) + sin(r)*sin(y)*sin(p))*F_tether_NED_y + cos(p)*sin(r)*F_tether_NED_z;
    double F_tether_BODY_z = (sin(y)*sin(r) + cos(r)*cos(y)*sin(p))*F_tether_NED_x + (cos(r)*sin(y)*sin(p) - cos(y)*sin(r))*F_tether_NED_y + cos(r)*cos(p)*F_tether_NED_z;

    if (X > X_0){
        force.x = force.x - F_tether_BODY_x;
        force.y = force.y - F_tether_BODY_y;
        force.z = force.z - F_tether_BODY_z;

    }






 */   // =======================================================================================================================
 // AKM: end

 // =======================================================================================================================
/*   // OPTION 3: WINCH CONTROL: DRAG MODE
   // =======================================================================================================================

   // set initialization distance
   initilization_distance = 5.0; //meters

   // set initial timing
   updating_time = true;

   // calculate distance to home (X)
   X = pow(pow(position.x, 2) + pow(position.y, 2) + pow(position.z, 2), 0.5); // in meters

   // if the plane mode is farther than defined distance, then we are not updating time
   if (X > initilization_distance) {
       updating_time = false;

   }
//the timer stops updating when X pass the initilization distance,
    if (updating_time) {
        update_time = AP_HAL::millis();
        X_0 = X+ 0.5*X;     //+0.05*X; // forcing to start with the cable model just when the initialization distance is passed
        X_0_prev = X_0;
                            //X_0 = X;  //winch is reeling out with the plane until reaching the initialization distance
//        Winch_omega = 0; // asume this for now, imagine tether is laid on the ground up to the initialization distance
//        Winch_alpha = 0; //
    }

    // calculate tether forces in NED coordinates
    if ( X_0 > X) { //X > initilization_distance &&: this applies the force since take-off, it should be small because airspeed= 0 then q ~= 0
        // Tether inclined cable model

        // MODIFY THE AREA!!! MUST BE THE CROSSECTION AREA OF THE CABLE LONGITUDINAL

        drag_tether = 0.5 * 1.225 * 0.0016*X_0 * 0.47 * pow(airspeed / 2, 2);
        weight_tether = 0.00194 * 9.8 * X_0;
        load_tether = drag_tether + (weight_tether * pow(pow(position.x,2)+pow(position.y,2),0.5)/X);
        // only if the plane X is smaller than X_0, add if condition.
        sag_tether = pow((3 * X*X_0 - 3 * pow(X, 2)) / 8, 0.5);
        Rx_tether = load_tether * pow(X, 2) / (8 * sag_tether);
        // calculate the forces in NED coordinates
        F_tether_NED_x = Rx_tether * (position.x / X);
        F_tether_NED_y = Rx_tether * (position.y / X);
        F_tether_NED_z = Rx_tether * (position.z / X);
        // rotate the forces to body axes
        F_tether_BODY_x = cos(y)*cos(p)*F_tether_NED_x + cos(p)*sin(y)*F_tether_NED_y - sin(p)*F_tether_NED_z;
        F_tether_BODY_y = (cos(y)*sin(r)*sin(p) - cos(r)*sin(y))*F_tether_NED_x + (cos(r)*cos(y) + sin(r)*sin(y)*sin(p))*F_tether_NED_y + cos(p)*sin(r)*F_tether_NED_z;
        F_tether_BODY_z = (sin(y)*sin(r) + cos(r)*cos(y)*sin(p))*F_tether_NED_x + (cos(r)*sin(y)*sin(p) - cos(y)*sin(r))*F_tether_NED_y + cos(r)*cos(p)*F_tether_NED_z;
    }

    else if (X > X_0) {
        // tether spring model
        E_times_area = 233230; // 116e9 Modulus times the area using diameter 1.6mm
        K_tether = E_times_area / X_0; // function of the distance to home (X)
        F_tether_NED_x = K_tether * (X - X_0) * (position.x / X);
        F_tether_NED_y = K_tether * (X - X_0) * (position.y / X);
        F_tether_NED_z = K_tether * (X - X_0) * (position.z / X);
        // rotate the forces to body axes
        F_tether_BODY_x = cos(y)*cos(p)*F_tether_NED_x + cos(p)*sin(y)*F_tether_NED_y - sin(p)*F_tether_NED_z;
        F_tether_BODY_y = (cos(y)*sin(r)*sin(p) - cos(r)*sin(y))*F_tether_NED_x + (cos(r)*cos(y) + sin(r)*sin(y)*sin(p))*F_tether_NED_y + cos(p)*sin(r)*F_tether_NED_z;
        F_tether_BODY_z = (sin(y)*sin(r) + cos(r)*cos(y)*sin(p))*F_tether_NED_x + (cos(r)*sin(y)*sin(p) - cos(y)*sin(r))*F_tether_NED_y + cos(r)*cos(p)*F_tether_NED_z;
    }


    // calculate the total tether force from the x,y,z tether forces in NED direction
    F_tether_total = pow(pow(F_tether_NED_x, 2) + pow(F_tether_NED_y, 2) + pow(F_tether_NED_z, 2), 0.5);

    // add tether forces to body axes
    force.x = force.x - F_tether_BODY_x;
    force.y = force.y - F_tether_BODY_y;
    force.z = force.z - F_tether_BODY_z;

    // calculate the tension range from the percentage (%) of maximum force
    // the maximum force should be the Lift when flying in crosswind.
    //F_tether_max = 500.0;
    //F_tether_perc = F_tether_total * 100 / F_tether_max;
    //F_tether_min = 33.3;
    //F_tether_int = 66.6;

    // MOTOR SPECS
    // the rated torque of the motor used is 7.7Nm
    // the motor inertia is 15.8 kg cm^2
    // resulting in rated alpha of 4873.4 rad/s^2


    // use a cubic curve to modify the angular acceleration of the winch
    // use points (F_tether, alpha) = (0,-4000),(15,-1000),(27.5,0) ,(40,1000),(55,4000)
    // alpha = 0.109091 F_tether^3 - 9 F_tether^2 + 310.455 F_tether - 4000
    //Winch_alpha = 0.109091*pow(F_tether_total, 3) - 9 * pow(F_tether_total, 2) + 310.455*F_tether_total - 4000; // run 9
    // Run 12: extending the poli of Run 9 to double F
    // 0.0136364 x^3 - 2.25 x^2 + 155.227 x - 4000
    // (0,-4000), (30,-1000), (55,0), (80,1000), (110,4000)
    //Winch_alpha = 0.0136364*pow(F_tether_total, 3) - 2.25*pow(F_tether_total, 2) + 155.227*F_tether_total - 4000;

    // Run 13 Compress 9th run poly by half F
    // 0.193439 x^3 - 16.9532 x^2 + 492.039 x - 3956.06
    //Winch_alpha = 0.193439*pow(F_tether_total, 3) - 16.9532*pow(F_tether_total, 2) + 492.039*F_tether_total - 3956.06;

    // Run 14: faster reel in than out
    // 0.168337 x^3 - 14.6737 x^2 + 443.617 x - 4023.29
    //Winch_alpha = 0.168337*pow(F_tether_total, 3) - 14.6737*pow(F_tether_total, 2) + 443.617*F_tether_total - 4023.29;

    //delta_time_delay = (AP_HAL::millis() - time_wait) * 0.001;
    //if (delta_time_delay > 0.02) {
    // =========================== DRAG MODE ==========================================================
    // Winch_omega = -0.0005*pow(F_tether_total, 3) - 0.0874*pow(F_tether_total, 2) + 5.6346*F_tether_total - 117.1468;
    Inertia_motor =  15.8/pow(100, 2); // motor inertia
    Torque_rated = -3.0; // rated torque
    b_motor = 0.1; // motor viscous friction
    Winch_radius = 0.075;
    delta_time_s = (AP_HAL::millis() - update_time) * 0.001;
    winch_gain = 1; // multiplies omega in the diff eqn. solution.

    F_tether_max = 50.0;
    F_tether_int = 50.0;
    F_tether_min = 25.0;

    m_omega = 5;
    limit_omega = 5;


    // linear omega as a function of tension (add Winch omega to have some kind of inertia in system.
    //Winch_omega = Winch_omega  + (0.00000007*pow(F_tether_total,3) - 0.000052189*pow(F_tether_total,2) + 0.012652*F_tether_total - 0.4297);
    if (F_tether_total > F_tether_max) {
        Winch_omega = Winch_omega + limit_omega;
    }
    else {
        Winch_omega = Winch_omega  - limit_omega;
    }
    // Holding fixed (X_0 - X) from defined value from FIX_R100_X0150_cable
    // X_0 = 1.5*X at all times.
    //X_0 = 1.2*X;
    //Winch_omega = (X_0 - X_0_prev) / (1.0/500.0*Winch_radius);
    //X_0_prev = X_0;


    // Linear alpha as a function of tension
    //Winch_alpha = m_alpha * F_tether_total - limit_alpha;


 //   if (F_tether_total > 0 && F_tether_total < F_tether_min) {
 //       Winch_alpha = -250;
 //   }
 //   else if (F_tether_total >= F_tether_min && F_tether_total < F_tether_int) {
 //       Winch_alpha = -100;
 //   }
 //   else if (F_tether_total >= F_tether_int && F_tether_total < F_tether_max) {
 //       Winch_alpha = 100;
 //   }
 //   else if (F_tether_total >= F_tether_max) {
 //       Winch_alpha = 250;
 //   }

    //Winch_omega = (Torque_rated - F_tether_total * Winch_radius - Winch_alpha * Inertia_motor) / b_motor;
    // need to set an initial update_time = now before starting to reel
    //Winch_omega = Winch_alpha * delta_time_s + Winch_omega;
//    if (updating_time) {
//        Winch_omega = (X_0 - X_0_prev) / (delta_time_s*Winch_radius);
//    }
//    else {
//        Winch_omega = 0;
//    }

    // C_4 = -Inertia_motor * exp(-b_motor * delta_time_s / Inertia_motor)*(Winch_omega - (Torque_rated - F_tether_total * Winch_radius) / b_motor) / b_motor;
    //Winch_omega = -winch_gain*((Torque_rated - F_tether_total * Winch_radius) / b_motor + exp(-2 * b_motor*delta_time_s / Inertia_motor)*(Winch_omega - (Torque_rated - F_tether_total * Winch_radius) / b_motor));
    // * b_motor*exp(-b_motor * delta_time_s / Inertia_motor) / Inertia_motor;
    // =========================== END DRAG MODE ==========================================================
    //time_wait = AP_HAL::millis();

    //}

    // ======================= PUMPING MODE ========================================================
 // Pumping mode with tether and winch.
    // separate the alpha curve into (+) only when reeling out
    //                               (-) only when reeling in
    X_0_min = 100;
    X_0_max = 350;

    // Where are we?
    // 1) Below X_0_min
    if (X_0 < X_0_min) {
        radius_min_reached = true;
        radius_max_reached = false;
    }
    // 2) Above X_0_max
    else if (X_0 >= X_0_max) {
        radius_min_reached = false;
        radius_max_reached = true;
    }

    // calculate the omega of the winch
    if (radius_min_reached == true && radius_max_reached == false) {
        // we are reeling out
    Winch_omega = 0.0002*pow(F_tether_total, 3) - 0.0691*pow(F_tether_total, 2) + 8.8606*F_tether_total;
    }
    else if (radius_min_reached == false && radius_max_reached == true) {
        // we are reeling in
        // winch alpha only negative
        Winch_omega = 0.0006*pow(F_tether_total, 3) - 0.0105*pow(F_tether_total, 2) + 1.8945*F_tether_total - 376.9911;
    }

    //Different alpha calculation for VTOL and Forward flight.
    // the Tension treshold is approx 15N

    if (F_tether_total< 15.0) {
        Winch_alpha = 0.0;
    }
        else
    {
        //0.486808 x^3 - 53.3432 x^2 + 1948.24 x - 22775.5
        Winch_alpha = 0.486808*pow(F_tether_total, 3) - 53.3432*pow(F_tether_total, 2) + 1948.24*F_tether_total - 22775.5;
    }

     // ======================= END PUMPING MODE =======================================================

    // y = 0.142424 x ^ 3 - 11.5518 x ^ 2 + 349.968 x - 4022.44
    //Winch_alpha = 0.142424* pow(F_tether_total, 3) - 11.5518*pow(F_tether_total, 2) + 349.968*F_tether_total - 4022.44; //run 10

    //Winch_alpha = -0.0030303*pow(F_tether_total, 3) + 0.25*pow(F_tether_total, 2) + 13.5985*F_tether_total - 500; // run 11


    // calculate the new output reeling speed

    //Winch_omega_max = 377.0/2;
    Winch_omega_max = 50;

    // define the constant reel-out command, Not seems to be working
    // Winch_omega = 0.1* Winch_omega_max;



    if (Winch_omega >= Winch_omega_max) {
        Winch_omega = Winch_omega_max;
    }
    else if (Winch_omega <= -Winch_omega_max) {
        Winch_omega = -Winch_omega_max;
    }

    Reeling_speed = Winch_omega * Winch_radius;
    //X_0_prev = X_0;


    // Comment here for deltaX model
    X_0 = X_0 + Reeling_speed * delta_time_s; // add the reeled length to the previous length

    // define tether length limits

    if (X_0 < 0.0) {
        X_0 = 0.0;
    }

    else if (X_0 > 200.0) {
        X_0 = 200.0;
    }

    // DIRECT CALCULATION OF X_0 SKIPPING THE ALPHA. THE ORIGINAL X_0 IS ABOVE
    //X_0 = X_0 + 1 / X + X;

    update_time = AP_HAL::millis();

    // =======================================================================================================================
 */  // END OF WINCH MODEL: DRAG MODE
 // =======================================================================================================================



    rot_accel = getTorque(aileron, elevator, rudder, thrust, force);

    if (have_launcher) {
        /*
          simple simulation of a launcher
         */
        if (launch_triggered) {
            uint64_t now = AP_HAL::millis64();
            if (launch_start_ms == 0) {
                launch_start_ms = now;
            }
            if (now - launch_start_ms < launch_time * 1000) {
                force.x += launch_accel;
                force.z += launch_accel / 3;
            }
        }
        else {
            // allow reset of catapult
            launch_start_ms = 0;
        }
    }

    // simulate engine RPM
    rpm1 = thrust * 7000;

    // scale thrust to newtons
    thrust *= thrust_scale;

    accel_body = Vector3f(thrust, 0, 0) + force;
    accel_body /= mass;

    // add some noise
    if (thrust_scale > 0) {
        add_noise(fabsf(thrust) / thrust_scale);
    }

    if (on_ground() && !tailsitter) {
        // add some ground friction
        Vector3f vel_body = dcm.transposed() * velocity_ef;
        accel_body.x -= vel_body.x * 0.3f;
    }
}


/*
  update the plane simulation by one time step
 */
void Plane::update(const struct sitl_input &input)
{
    Vector3f rot_accel;

    update_wind(input);

    calculate_forces(input, rot_accel, accel_body);

    update_dynamics(rot_accel);

    // update lat/lon/altitude
    update_position();
    time_advance();

    // update magnetic field
    update_mag_field_bf();
}