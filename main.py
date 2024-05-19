import numpy as np
from matplotlib import pyplot as plt
import math
import networkx

#traction effort curve: gives the maximum tractive effort that can be applied at a speed
#resistance is a function of velocity and curvature (which we can neglect, as we assume there is none). So, resistance
#Will be completely dependent on speed i.e, not influenced immediately by instantaneous tractive effort.
#Assumptions: No loss due to friction
def train_force_resistance_curve(speed):
    #approximation based on Fig2 of "A Power-Management Strategy for Multiple Unit Railroad Vehicles"
    #x-coordinate is speed in km/hr
    #y-coordinate is kN

    tractive_effort = 0
    if speed <90:
        tractive_effort = 200

    if 90 <= speed and speed <= 200:
        tractive_effort = -1*speed + 290

    if speed>200:
        tractive_effort = 90

    #Crude approximation of speed and resistance, according to figure 2
    resistance = 0.2*speed + 5

    return [tractive_effort,resistance]

def motor_energy_efficiency_curve(motor_energy_output):
    if motor_energy_output >= 0 and motor_energy_output <= 1000:
        motor_efficiency = 0.8*np.exp((-1*(motor_energy_output +300-1000)**2)/(2*(460)**2))
    else:
        motor_efficiency = 0.8*np.exp((-1*(1000)**2)/(2*(460)**2))
    return motor_efficiency

#Returns the braking rate information, based on speed. Speed of 100 km/hr chosen for the cutoff.
def electric_maximum_breaking_rate(train_speed):
    if train_speed < 90:
        return 1

    if train_speed >= 90:
        electric_braking_rate = -1*(1/110)*(train_speed - 90)
        return electric_braking_rate

def convert_vehicle_speed_into_ang_speed(vehicle_speed, wheel_radius):
    return vehicle_speed/wheel_radius

#kinematic equations
def final_velocity(initial_velocity, acceleration,time):
    return initial_velocity + acceleration*time

def position_in_1_dimension(initial_position, initial_vel, time, acceleration):
    return initial_position + initial_vel*time + 0.5*acceleration*(time)**2




class train():
    #all units are metric
    def __init__(self):
        self.tare_mass = 1000
        self.load_mass = 500
        self.effective_mass_constant = 0.1
        self.effective_mass = self.tare_mass*(1 + 0.1) + self.load_mass
        self.tangential_velocity = 0
        self.current_acceleration = 0
        self.gravitational_constant_newtons = 6.6743*(10**-11)
        self.angle_of_incline = 0
        self.wheel_radius =  0.762/2
        self.angular_speed = 0
        self.current_tractive_effort = 0
        self.current_train_resistance = 0
        self.vehicle_speed = 0
        self.powertrain_efficiency_coefficient = 0.9
        #kilowats is max power output number
        self.max_power_output = 3400
        self.time_increment = 10
        self.total_wheel_power = 0
        self.total_current_engine_output = 0
        self.engine_idle_power_cost = 500
        self.train_acceleration = 0
        self.total_power_expended = 0
        self.train_position = 0
        self.state_building_speed = True
        self.time_elapsed = 0
        self.total_energy_in_system = 0
        self.energy_added_in_step = 0
        self.total_engine_energy_consumption = 0
        self.current_effective_tractive_effort = 0
        self.state_braking = False
        self.force_total_braking = 0
        self.current_breaking_force = 0
        self.proportion_electrical_break = 0
        self.proportion_mechanical_break = 0 


    #Uses fig 2, traction force curve/speed in Km/hr and compares with
    # power output curve, to estimate motor power expenditure.
    def estimate_current_diesel_motor_energy(self, vehicle_speed):
        #incorrect
        #if vehicle speed is above 0, set vehicle speed to slightly above zero, to calculate power using the power calculation curve (to avoid a divide by zero error)
        if vehicle_speed == 0:
            vehicle_speed = 2
        self.angular_speed = convert_vehicle_speed_into_ang_speed(vehicle_speed, self.wheel_radius)
        self.total_wheel_power = self.angular_speed*self.current_effective_tractive_effort*self.wheel_radius
        self.energy_added_in_step = np.abs(self.total_energy_in_system - self.total_wheel_power)
        self.total_energy_in_system = self.total_wheel_power
        #ouput power efficiency graph is related to rates
        return #self.energy_added_in_step/motor_energy_efficiency_curve(self.energy_added_in_step)

    #Equation of motion for electric trains, returns acceleration.
    def lomonoffs_equation_of_motion_acc(self, effective_mass,tare_mass, slope_angle, vehicle_resistance, tractive_effort, acceleration_from_g):
        return (tractive_effort - vehicle_resistance - tare_mass*acceleration_from_g*np.sin(slope_angle))/effective_mass

    #in our simulation, we assume the train utilizes the maximum tractive effort possible, according to the tractive effort curve, finally
    #achieving the top speed.
    #Then, when it comes time to break, the train idles and uses either the regenerative braking or the forward breaking.
    def simulation(self):
        if self.state_building_speed:
            result = train_force_resistance_curve(self.vehicle_speed)
            self.current_tractive_effort = result[0]
            self.current_train_resistance = result[1]
            self.train_acceleration = self.lomonoffs_equation_of_motion_acc(self.effective_mass,self.tare_mass,self.angle_of_incline,self.current_train_resistance,self.current_tractive_effort,self.gravitational_constant_newtons)
            self.current_effective_tractive_effort = self.effective_mass*self.train_acceleration
            self.current_kinetic_energy = 0.5*self.effective_mass*(self.vehicle_speed)**2

            #update vehicle speed, based on traction effort
            #use accel to update position, speed, with kinematic equations
            self.train_position = position_in_1_dimension(self.train_position,self.vehicle_speed,self.time_increment,self.train_acceleration)
            self.vehicle_speed = final_velocity(self.vehicle_speed, self.train_acceleration, self.time_increment)

            self.time_elapsed += self.time_increment

        if self.state_braking:
            self.current_train_resistance = train_force_resistance_curve(self.vehicle_speed)[1]
            self.current_tractive_effort = 0
            self.current_tractive_effort = -1*self.current_breaking_force


        return [self.time_elapsed, self.train_acceleration, self.train_position, self.vehicle_speed]
    
    #sets the proportion of electrical and mechanical breaking
    def set_breaking_combination(self):
        return

    def run_simulation(self):
        curr_time = []
        train_acc_arr = []
        train_pos_arr = []
        train_speed_arr = []
        #bring vehicle up to 150 km/hr
        while self.vehicle_speed <= 150:
            sim_variables = self.simulation()
            curr_time += [sim_variables[0]]
            train_acc_arr += [sim_variables[1]]
            train_pos_arr += [sim_variables[2]]
            train_speed_arr += [sim_variables[3]]

        print(train_acc_arr)
        #set traction effort to 0, to represent engine idling, change to braking state
        self.current_tractive_effort = 0
        self.state_building_speed = False
        self.state_braking = True
        return

    # Generate curves associated with braking speeds for electric braking/forward braking.
    def generate_braking_curves(self):

        return

    #routine for testing the breaking physics.
    #To determine the range, we need to use the maximum breaking rate.
    #Forwards and backwards calculations refer to the fastest you can stop and the latest you can stop. I.e., what path does the train
    #take when you start maximal breaking immediately
    #and what path does the train take when you start your maximal breaking at the last moment possible to stop on time

    def breaking_routine(self):

        return


if __name__ == '__main__':
    electric_train = train()
    electric_train.run_simulation()
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
