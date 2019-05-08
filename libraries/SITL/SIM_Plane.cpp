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
    } else if (strstr(frame_str, "-vtail")) {
        vtail = true;
    } else if (strstr(frame_str, "-dspoilers")) {
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
    if (alpha-alpha0 > max_alpha_delta) {
        alpha = alpha0 + max_alpha_delta;
    } else if (alpha0-alpha > max_alpha_delta) {
        alpha = alpha0 - max_alpha_delta;
    }
	double sigmoid = ( 1+exp(-M*(alpha-alpha0))+exp(M*(alpha+alpha0)) ) / (1+exp(-M*(alpha-alpha0))) / (1+exp(M*(alpha+alpha0)));
	double linear = (1.0-sigmoid) * (c_lift_0 + c_lift_a0*alpha); //Lift at small AoA
	double flatPlate = sigmoid*(2*copysign(1,alpha)*pow(sin(alpha),2)*cos(alpha)); //Lift beyond stall

	float result  = linear+flatPlate;
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
    
	double AR = pow(b,2)/s;
	double c_drag_a = c_drag_p + pow(c_lift_0+c_lift_a0*alpha,2)/(M_PI*oswald*AR);

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

	double qbar = 1.0/2.0*rho*pow(effective_airspeed,2)*s; //Calculate dynamic pressure
	double la, na, ma;
	if (is_zero(effective_airspeed))
	{
		la = 0;
		ma = 0;
		na = 0;
	}
	else
	{
		la = qbar*b*(c_l_0 + c_l_b*beta + c_l_p*b*p/(2*effective_airspeed) + c_l_r*b*r/(2*effective_airspeed) + c_l_deltaa*inputAileron + c_l_deltar*inputRudder);
		ma = qbar*c*(c_m_0 + c_m_a*alpha + c_m_q*c*q/(2*effective_airspeed) + c_m_deltae*inputElevator);
		na = qbar*b*(c_n_0 + c_n_b*beta + c_n_p*b*p/(2*effective_airspeed) + c_n_r*b*r/(2*effective_airspeed) + c_n_deltaa*inputAileron + c_n_deltar*inputRudder);
	}


	// Add torque to to force misalignment with CG
	// r x F, where r is the distance from CoG to CoL
	la +=  CGOffset.y * force.z - CGOffset.z * force.y;
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
	double c_x_a = -c_drag_a*cos(alpha)+c_lift_a*sin(alpha);
	double c_x_q = -c_drag_q*cos(alpha)+c_lift_q*sin(alpha);
	double c_z_a = -c_drag_a*sin(alpha)-c_lift_a*cos(alpha);
	double c_z_q = -c_drag_q*sin(alpha)-c_lift_q*cos(alpha);

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
	double qbar = 1.0/2.0*rho*pow(airspeed,2)*s; //Calculate dynamic pressure
	double ax, ay, az;
	if (is_zero(airspeed))
	{
		ax = 0;
		ay = 0;
		az = 0;
	}
	else
	{
		ax = qbar*(c_x_a + c_x_q*c*q/(2*airspeed) - c_drag_deltae*cos(alpha)*fabs(inputElevator) + c_lift_deltae*sin(alpha)*inputElevator);
		// split c_x_deltae to include "abs" term
		ay = qbar*(c_y_0 + c_y_b*beta + c_y_p*b*p/(2*airspeed) + c_y_r*b*r/(2*airspeed) + c_y_deltaa*inputAileron + c_y_deltar*inputRudder);
        // AKM: adding distance_to_home as force to see how it bahaves
        // the force must increase as the distance gets bigger.
        az = qbar * (c_z_a + c_z_q * c*q / (2 * airspeed) - c_drag_deltae * sin(alpha)*fabs(inputElevator) - c_lift_deltae * cos(alpha)*inputElevator);// + distance_to_home;
		// split c_z_deltae to include "abs" term
	}
    return Vector3f(ax, ay, az);
}

void Plane::calculate_forces(const struct sitl_input &input, Vector3f &rot_accel, Vector3f &body_accel)
{
    float aileron  = filtered_servo_angle(input, 0);
    float elevator = filtered_servo_angle(input, 1);
    float rudder   = filtered_servo_angle(input, 3);
    bool launch_triggered = input.servos[6] > 1700;
    float throttle;
    if (reverse_elevator_rudder) {
        elevator = -elevator;
        rudder = -rudder;
    }
    if (elevons) {
        // fake an elevon plane
        float ch1 = aileron;
        float ch2 = elevator;
        aileron  = (ch2-ch1)/2.0f;
        // the minus does away with the need for RC2_REV=-1
        elevator = -(ch2+ch1)/2.0f;

        // assume no rudder
        rudder = 0;
    } else if (vtail) {
        // fake a vtail plane
        float ch1 = elevator;
        float ch2 = rudder;
        // this matches VTAIL_OUTPUT==2
        elevator = (ch2-ch1)/2.0f;
        rudder   = (ch2+ch1)/2.0f;
    } else if (dspoilers) {
        // fake a differential spoiler plane. Use outputs 1, 2, 4 and 5
        float dspoiler1_left = filtered_servo_angle(input, 0);
        float dspoiler1_right = filtered_servo_angle(input, 1);
        float dspoiler2_left = filtered_servo_angle(input, 3);
        float dspoiler2_right = filtered_servo_angle(input, 4);
        float elevon_left  = (dspoiler1_left + dspoiler2_left)/2;
        float elevon_right = (dspoiler1_right + dspoiler2_right)/2;
        aileron  = (elevon_right-elevon_left)/2;
        elevator = (elevon_left+elevon_right)/2;
        rudder = fabsf(dspoiler1_right - dspoiler2_right)/2 - fabsf(dspoiler1_left - dspoiler2_left)/2;
    }
    printf("Aileron: %.1f elevator: %.1f rudder: %.1f\n", aileron, elevator, rudder);

    if (reverse_thrust) {
        throttle = filtered_servo_angle(input, 2);
    } else {
        throttle = filtered_servo_range(input, 2);
    }
    
    float thrust     = throttle;

    if (ice_engine) {
        thrust = icengine.update(input);
    }

    // calculate angle of attack
    angle_of_attack = atan2f(velocity_air_bf.z, velocity_air_bf.x);
    beta = atan2f(velocity_air_bf.y,velocity_air_bf.x);

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

    // =======================================================================================================================
// AKM: create the variable tether length constraint (look in Navigation.cpp)
    /*
     // OPTION 1: PUMPING MODE, CONSTANT REEL SPEED
    // ^^^^^^^^^^^^^^^^^^^^^^
    // initialization distance, replaces the control_mode set initializer 
    initilization_distance = 20.0; //meters
    // define the minimun and maximum radius
    radius_min = 50; // meters, make sure this is the same as Navigation.cpp
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
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     // COMMENT OUT HERE THE PUMPING MODE
    */

    /*
    // OPTION 2: FIXED TETHER LENGTH
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    X = pow(pow(position.x, 2) + pow(position.y, 2) + pow(position.z, 2), 0.5); // in meters
    //X_0 = 200; // meters, this must be equal to S1_in_S2.S2_radius_cm in Commands_logic.cpp
    X_0 = 200 - 80; // this is the R_sphere minus the oscillation of the plane distance with K=0
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // COMMENT OUT HERE THE FIXED TETHER LENGRTH
    */



    // how to calculate the tether constant of elasticity
/*
R_sphere = 240; %m
E_tether = 116E9; % Pa
d_tether = 1.6E-3;  % m
A_tether = pi* d_tether^2 /4; %m^2
K_tether = E_tether * A_tether / R_sphere;
    */
    /*
    if (X_0 >= radius_min) {
        K_tether = 116000000000 * 0.000002 / X_0;
    }
    else {
        K_tether = 116000000000 * 0.000002 / radius_min;
    }
    */

    /*
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

    /
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






    // =======================================================================================================================
    // AKM: end
    */

    // OPTION 3: WINCH CONTROL: DRAG MODE

    // set initial timing
    updating_time = true;

    // calculate distance to home (X)
    X = pow(pow(position.x, 2) + pow(position.y, 2) + pow(position.z, 2), 0.5); // in meters

    // if the plane mode is farther than defined distance, then we are not updating time
    if (X >= initilization_distance) {
        updating_time = false;
    }
//the timer stops updating when X pass the initilization distance,
    if (updating_time) {
        update_time = AP_HAL::millis();
        X_0 = X;  //winch is reeling out with the plane until reaching the initialization distance
        Winch_omega = 0;
    }

    // calculate tether forces in NED coordinates
    K_tether = 0.7;
    double F_tether_NED_x = K_tether * (X - X_0) * (position.x / X);
    double F_tether_NED_y = K_tether * (X - X_0) *(position.y / X);
    double F_tether_NED_z = K_tether * (X - X_0) *(position.z / X);

    // rotate the forces to body axes
    double F_tether_BODY_x = cos(y)*cos(p)*F_tether_NED_x + cos(p)*sin(y)*F_tether_NED_y - sin(p)*F_tether_NED_z;
    double F_tether_BODY_y = (cos(y)*sin(r)*sin(p) - cos(r)*sin(y))*F_tether_NED_x + (cos(r)*cos(y) + sin(r)*sin(y)*sin(p))*F_tether_NED_y + cos(p)*sin(r)*F_tether_NED_z;
    double F_tether_BODY_z = (sin(y)*sin(r) + cos(r)*cos(y)*sin(p))*F_tether_NED_x + (cos(r)*sin(y)*sin(p) - cos(y)*sin(r))*F_tether_NED_y + cos(r)*cos(p)*F_tether_NED_z;


    // calculate the total tether force from the x,y,z tether forces in NED direction
if (X > X_0) {
    F_tether_total = pow(pow(F_tether_NED_x, 2) + pow(F_tether_NED_y, 2) + pow(F_tether_NED_z, 2), 0.5);
    // add tether forces to body axes
    force.x = force.x - F_tether_BODY_x;
    force.y = force.y - F_tether_BODY_y;
    force.z = force.z - F_tether_BODY_z;
}
else
{
    F_tether_total = 0;
}
 
    // calculate the tension range from the percentage (%) of maximum force
    F_tether_max = 40.0;
    F_tether_perc = F_tether_total * 100 / F_tether_max;
    F_tether_zero = 0.0;

    // determine the reeling angular acceleration
    if (F_tether_perc <= F_tether_zero) {
        Winch_alpha = -750; // change this for realistic value
    }
    else if (F_tether_perc > 0.0 && F_tether_perc <= 35.0) {
        Winch_alpha = 0;
    }
    else if (F_tether_perc > 35.0 && F_tether_perc <= 70.0) {
        Winch_alpha = 750.0;
    }
    else if (F_tether_perc > 70.0) {
        Winch_alpha = 1500.0;
    }


    // calculate the new output reeling speed
    Winch_radius = 0.075;
    Winch_omega_max = 377.0;
    delta_time_s = (AP_HAL::millis() - update_time) * 0.001; // need to set an initial update_time = now before starting to reel
    Winch_omega = Winch_alpha * delta_time_s + Winch_omega;

    if (Winch_omega >= Winch_omega_max) {
        Winch_omega = Winch_omega_max;
    }
    else if (Winch_omega <= -Winch_omega_max) {
        Winch_omega = -Winch_omega_max;
    }


    Reeling_speed = Winch_omega * Winch_radius;
    X_0 = X_0 + Reeling_speed * delta_time_s;

    update_time = AP_HAL::millis();


    // END OF WINCH MODEL: DRAG MODE
    //=============================================================================




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
            if (now - launch_start_ms < launch_time*1000) {
                force.x += launch_accel;
                force.z += launch_accel/3;
            }
        } else {
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