import random
import numpy as np
from matplotlib import pyplot as plt
import math
import networkx as nx

#Force returned is in newtons.
#speed in km/hr
def tractive_force_curve(speed):
    if speed < 90:
        tractive_effort = 200000
    else:
        tractive_effort = 200000*np.exp(-1.9*(speed -90)/137.77)

    return tractive_effort

def resistance_curve(speed):
    if speed < 0:
        resistance = 0
    else:
        resistance = (0.2*speed + 5)*1000

    return resistance

def electric_braking_rate_curve(speed):
    if speed < 60:
        regen_brake_force = 350000
    else:
        regen_brake_force = 350000*np.exp(-1*(speed - 60 )/30)

    return regen_brake_force

def final_velocity(initial_velocity, acceleration, time_step):
    return initial_velocity + acceleration*time_step

def position_in_1_dimension(initial_position,initial_vel, time_step, acceleration):
    return initial_position + initial_vel*(time_step) + 0.5*acceleration*(time_step**2)

def plot_curve(function, lower_bound,upper_bound, partitions):
    array_func = []
    array_step = []
    step_size = (upper_bound-lower_bound)/partitions
    step = lower_bound
    while step <= upper_bound:
        array_func += [function(step)]
        array_step += [step]
        step += step_size

    plt.plot(array_step,array_func)
    return

class train():
    def __init__(self):
        #masses here are in kg
        self.tare_mass = 50000
        self.rotary_allowance = 0.1
        self.freight_load = 200000
        self.effective_mass = (1+self.rotary_allowance)*self.tare_mass + self.freight_load
        self.current_velocity = 0
        self.current_distance_from_origin = 0
        self.current_acceleration = 0
        self.current_max_avail_elc_brak_forc = 0
        self.constant_braking_force = 350000
        self.current_max_attainable_tract_effort = 0
        self.current_resistance = 0
        self.time_step = 5
        self.accel_array_curve_vel = []
        self.accel_array_curve_pos = []
        self.coast_array_curve_vel = []
        self.coast_array_curve_pos = []
        self.max_decel_array_curve_vel = []
        self.max_decel_array_curve_pos = []
        self.elec_decel_array_curve_vel = []
        self.elec_decel_array_curve_pos = []
        self.max_decel_array_curve_vel2 = []
        self.max_decel_array_curve_pos2 = []
        self.braking_starting_vel = 0
        self.braking_starting_post = 0
        self.coasting_range = 23000

    #assumes the train operates on a 1-dimensional flat track, hence, the last term, with sin(theta), is not present
    def lomonosoffs_equation(self, tractive_force_app, resistance):
        return (tractive_force_app - resistance)/self.effective_mass

    def simulation_step(self, train_state):
        self.current_max_attainable_tract_effort = tractive_force_curve(self.current_velocity)
        self.current_resistance = resistance_curve(self.current_velocity)
        self.current_max_avail_elc_brak_forc = electric_braking_rate_curve(self.current_velocity)

        if train_state == "acceleration":
            self.current_acceleration = self.lomonosoffs_equation(self.current_max_attainable_tract_effort, self.current_resistance)
        if train_state == "braking_full":
            total_force_wo_resist = (self.constant_braking_force + self.current_max_avail_elc_brak_forc)*-1
            self.current_acceleration = self.lomonosoffs_equation(total_force_wo_resist, self.current_resistance)
        if train_state == "coasting":
            self.current_acceleration = self.lomonosoffs_equation(0, self.current_resistance)
        if train_state == "braking_elec":
            total_force_wo_resist = self.current_max_avail_elc_brak_forc*-1
            self.current_acceleration = self.lomonosoffs_equation(total_force_wo_resist, self.current_resistance)

        self.current_distance_from_origin = position_in_1_dimension(self.current_distance_from_origin, self.current_velocity, self.time_step, self.current_acceleration)
        self.current_velocity = final_velocity(self.current_velocity,self.current_acceleration, self.time_step)

        return

    def simulation(self):
        self.accel_array_curve_pos += [self.current_distance_from_origin]
        self.accel_array_curve_vel += [self.current_velocity]
        while self.current_velocity <= 150:
            self.simulation_step("acceleration")
            self.accel_array_curve_pos += [self.current_distance_from_origin]
            self.accel_array_curve_vel += [self.current_velocity]

        plt.plot(self.accel_array_curve_pos,self.accel_array_curve_vel)

        self.braking_starting_vel = self.current_velocity
        self.braking_starting_post = self.current_distance_from_origin

        self.max_decel_array_curve_pos += [self.braking_starting_post]
        self.max_decel_array_curve_vel += [self.braking_starting_vel]
        while self.current_velocity > 0:
            self.simulation_step("braking_full")
            self.max_decel_array_curve_pos += [self.current_distance_from_origin]
            self.max_decel_array_curve_vel += [self.current_velocity]

        plt.plot(self.max_decel_array_curve_pos,self.max_decel_array_curve_vel)
        self.current_velocity = self.braking_starting_vel
        self.current_distance_from_origin = self.braking_starting_post

        self.elec_decel_array_curve_pos += [self.current_distance_from_origin]
        self.elec_decel_array_curve_vel += [self.current_velocity]

        while self.current_velocity > 0:
            self.simulation_step("braking_elec")
            self.elec_decel_array_curve_pos += [self.current_distance_from_origin]
            self.elec_decel_array_curve_vel += [self.current_velocity]
        plt.plot(self.elec_decel_array_curve_pos, self.elec_decel_array_curve_vel)

        self.current_velocity = self.braking_starting_vel
        self.current_distance_from_origin = self.braking_starting_post

        self.coast_array_curve_pos += [self.braking_starting_post]
        self.coast_array_curve_vel += [self.braking_starting_vel]

        while self.current_distance_from_origin <= self.braking_starting_post + self.coasting_range:
            self.simulation_step("coasting")
            self.coast_array_curve_pos += [self.current_distance_from_origin]
            self.coast_array_curve_vel += [self.current_velocity]

        plt.plot(self.coast_array_curve_pos,self.coast_array_curve_vel)

        self.max_decel_array_curve_pos2 += [self.current_distance_from_origin]
        self.max_decel_array_curve_vel2 += [self.current_velocity]

        while self.current_velocity > 0:
            self.simulation_step("braking_full")
            self.max_decel_array_curve_pos2 += [self.current_distance_from_origin]
            self.max_decel_array_curve_vel2 += [self.current_velocity]

        plt.plot(self.max_decel_array_curve_pos2,self.max_decel_array_curve_vel2)
        plt.show()




        return


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    electric_train = train()
    plot_curve(tractive_force_curve,0,150,400)
    plot_curve(resistance_curve,0,150,400)
    plot_curve(electric_braking_rate_curve,0,150,400)
    plt.close()
    electric_train.simulation()
    print("steve")


