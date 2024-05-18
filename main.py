import numpy as np
from matplotlib import pyplot as plt
import math

#traction effort curve: gives the maximum tractive effort that can be applied at a speed
#resistance is a function of velocity and curvature (which we can neglect, as we assume there is none). So, resistance
#Will be completely dependent on speed i.e, not influenced immediately by instantaneous tractive effort.
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

def motor_power_efficiency_curve(motor_power_output):
    if motor_power_output >= 0 and motor_power_output <= 2000:
        motor_efficiency = 0.3*np.exp((-1*(motor_power_output-1000)**2)/(2*(460)**2))
    return motor_efficiency

def convert_vehicle_speed_into_ang_speed(vehicle_speed, wheel_radius):
    return vehicle_speed/wheel_radius



class train():
    #all units are metric
    def __init__(self):
        self.tare_mass = 4000
        self.load_mass = 2000
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
        self.time_increment = 1
        self.total_wheel_power = 0
        self.total_current_engine_output = 0
        self.engine_idle_power_cost = 500
        self.train_acceleration = 0
        self.total_power_expended = 0
        self.train_position = 0
        self.state_building_speed = True


    #Uses fig 2, traction force curve/speed in Km/hr and compares with
    # power output curve, to estimate motor power expenditure.
    def estimate_current_diesel_motor_power(self, vehicle_speed):
        #if vehicle speed is above 0, set vehicle speed to slightly above zero, to calculate power using the power calculation curve (to avoid a divide by zero error)
        if vehicle_speed == 0:
            vehicle_speed = 2
        self.angular_speed = convert_vehicle_speed_into_ang_speed(vehicle_speed, self.wheel_radius)
        self.total_wheel_power = self.angular_speed*self.current_tractive_effort
        self.total_current_engine_output = self.total_wheel_power/self.powertrain_efficiency_coefficient
        print(self.total_current_engine_output)
        return self.total_current_engine_output/motor_power_efficiency_curve(self.total_current_engine_output)

    #Equation of motion for electric trains, returns acceleration.
    def lomonoffs_equation_of_motion(self, effective_mass,tare_mass, slope_angle, vehicle_resistance, tractive_effort, acceleration_from_g):
        return (tractive_effort - vehicle_resistance - tare_mass*acceleration_from_g*np.sin(slope_angle))/effective_mass

    #in our simulation, we assume the train utilizes the maximum tractive effort possible, according to the tractive effort curve, finally
    #achieving the top speed.
    #Then, when it comes time to break, the train idles and uses either the regenerative braking or the forward breaking.
    def simulation(self):
        if self.state_building_speed:
            result = train_force_resistance_curve(self.vehicle_speed)
            self.current_tractive_effort = result[0]
            self.current_train_resistance = result[1]
            self.current_engine_power_consumption = self.estimate_current_diesel_motor_power(self.vehicle_speed)
            self.total_power_expended += self. current_engine_power_consumption*self.time_increment
            #update vehicle speed, based on traction effort
            self.train_acceleration = self.lomonoffs_equation_of_motion(self.effective_mass,self.tare_mass,self.angle_of_incline,self.current_train_resistance,self.current_tractive_effort,self.gravitational_constant_newtons)
            #use equation of motion to update position, speed, with kinematic equations
            
        return


if __name__ == '__main__':
    electric_train = train()
    electric_train.simulation()
