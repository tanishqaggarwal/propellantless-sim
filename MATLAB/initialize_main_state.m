function main_state = initialize_main_state(seed,condition)
%initialize_main_state Sample an initial main state given a seed and condition
%   To set the seed based on current time, use seed = 'shuffle'
global const

main_state=struct();
rng(seed,'threefry');
main_state.leader=struct();
main_state.follower=struct();
[dynamics,actuators,sensors] = initialize_states(condition);
main_state.leader.dynamics=dynamics;
main_state.leader.actuators=actuators;
main_state.leader.sensors=sensors;
[dynamics,actuators,sensors] = initialize_states(condition);
main_state.follower.dynamics=dynamics;
main_state.follower.actuators=actuators;
main_state.follower.sensors=sensors;

end

function [dynamics,actuators,sensors] = initialize_states(condition)
%initialize_states Samples the initial satellite state of on satellite. 
% dynamics,actuators, and sensors are described in the MATLAB readme
global const
dynamics=struct();
actuators=struct();
sensors=struct();

% dynamics
dynamics.time_ns = int64(0);% int64
dynamics.time= double(dynamics.time_ns)*1E-9;
a  = 6860636.6;  % Semimajor axis                        (m)
e  = 0.001;      % Eccentricity                          (unitless)
i  = 45*pi/180;  % Inclination angle                     (rad)
O  = 0.0;        % Right ascension of the ascending node (rad)
o  = 0.0;        % Argument of perigee                   (rad)
nu = 0*pi/180;   % True anamoly                          (rad)

[   r,...  % Position (m)   [eci]
    v,...  % Velocity (m/s) [eci]
] = utl_orb2rv(a*(1-e), e, i, O, o, nu, const.mu);
dynamics.position_eci= r;
dynamics.velocity_eci= v;


switch condition
    case 'detumbled'
        dynamics.angular_rate_body= [0;0;0;];
    otherwise
        dynamics.angular_rate_body=randn(3,1)*5*pi/180;
end
dynamics.quat_body_eci=randn(4,1);
dynamics.quat_body_eci=dynamics.quat_body_eci/norm(dynamics.quat_body_eci);
dynamics.wheel_rate_body=[0;0;0;];
dynamics.fuel_net_angular_momentum_eci=[0;0;0;];
dynamics.fuel_mass=0.16;

% actuators
actuators= actuators_off_state();

% sensors
sensors.gyro_bias= zeros(3,1);
sensors.magnetometer_bias= zeros(3,1);
sensors.sunsensor_real_normals= transpose([ 0.9397	0.3420      0
                                    0.9397	-0.3420 	0
                                    0.9397	0           0.3420
                                    0.9397	0           -0.3420
                                    -0.9397	0.3420      0
                                    -0.9397	-0.3420     0
                                    -0.9397	0           0.3420
                                    -0.9397	0           -0.3420
                                    0.3420	0.9397      0
                                    -0.3420 0.9397      0
                                    0       0.9397      0.3420
                                    0       0.9397      -0.3420
                                    0.3420	-0.9397      0
                                    -0.3420	-0.9397      0
                                    0       -0.9397      0.3420
                                    0       -0.9397      -0.3420
                                    0.3420	0           0.9397
                                    -0.3420	0           0.9397
                                    0       0.3420      0.9397
                                    0       -0.3420     0.9397]);
sensors.sunsensor_measured_normals= sensors.sunsensor_real_normals;
sensors.gps_bias= zeros(6,1);
sensors.gps_time_till_lock= const.GPS_LOCK_TIME;
sensors.cdgps_bias= zeros(6,1);
sensors.cdgps_time_till_lock= const.CDGPS_LOCK_TIME;
end

