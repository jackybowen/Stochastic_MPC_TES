#!/usr/bin/env python3

import torch
import numpy as np
import os
import yaml
import math
import json
import pandas as pd

from tcilc_model import simple_dnn as DNN

def get_prediction(u,
            model_ckpt='atlanta_model.ckpt',
            torch_config='atlanta_config.yaml'):

    '''
    The main function called for receiving a TCILC model prediction. Loads model and sets up configuration.

        Parameters:
            u (dict): {
                outT_arr (float, length:n): [...], Outdoor temperature for the future n steps,
                Tset_arr (float, length: [1,n]): [...], Current zone cooling temperature setpoint for choice length of timesteps
                T_arr (float, length: [1,n]): [...], Current zone temperature setpoint for choice length of timesteps
                outT (float, scalar): Current outdoor temperature (currently not being used)
                half_hour (int, scalar): Index of half hour through a week. Range = [0,240).
                }
            model_ckpt (str): Name of pytorch checkpoint file (.ckpt). Should be in directory with this file.,
            config (str): Name of yaml config file associated with the pytorch checkpoint (.yaml). Should be in directory with this file.,
            power_mean (float): Mean power flexibility value used for standardization
            power_std (float): Standard deviation power flexibility value used for standardization
            temp_mean (float): Mean outdoor air temp value used for standardization
            temp_std (float): Standard deviation outdoor air temp value used for standardization

        Returns:
            y (dict): {
                base_power_arr (float, length:n): [...], Baseline power (zone temperature setpoints are kept uncahnged) for the future n steps,
                power_upper_arr (float, length:n): [...], Power flexibility upper bounds for the future n steps,
                power_lower_arr (float, length: n): [...], Power flexibility lower bounds for the future n steps,
                }
    '''

    dirname = os.path.dirname(__file__) # Absolute path to the directory of tcilc_model.py
    if 'atlanta' in model_ckpt:
        power_mean, power_std, temp_mean, temp_std = get_norm_vals('atlanta')
    elif 'denver' in model_ckpt:
        power_mean, power_std, temp_mean, temp_std = get_norm_vals('denver')
    else:
        raise Exception('Expected "denver" or "atlanta" as a parameter to get_norm_vals(). Check that you have correct ckpt and yaml files.')

    checkpoint = os.path.join(dirname, model_ckpt)
    conf = os.path.join(dirname, torch_config)
    config = convert_conf(conf)

    model, device = load_model(checkpoint, config)

    return predict(model,
            device,
            u,
            power_mean,
            power_std,
            temp_mean,
            temp_std)

def convert_conf(config):
    '''
    Converts pytorch config file to dictionary for easier access.

        Parameters:
            config (str): Path to pytorch configuartion file for the model.

        Returns:
            results (dict): Various hyperparameters.
    '''

    with open(config, "r") as yamlfile:
        data = yaml.load(yamlfile, Loader=yaml.FullLoader)

    results = {
            'model': data['model/_target_']['value'],
            'which_city': data['dataset/which_city']['value'],
            'layer_specs': data['model/layer_specs']['value'],
            'lr': data['optimizer/lr']['value'],
            'optimizer': data['optimizer/optimizer']['value'],
            'input_size': data['model/input_size']['value'],
            'output_size': data['model/output_size']['value'],
            'batch_size': data['dataloader/batch_size']['value'],
            'criterion': data['model/criterion']['value'],
            }

    return results

def get_norm_vals(city):
    '''
    Finds the values used to normalize the data for the given city.

        Parameters:
            city (str): The city that predictions will be made for.

        Returns:
            power_mean (float): Mean of total power in the "city" dataset.
            power_std (float): Standard deviation of total power in the "city" dataset.
            temp_mean (float): Mean of outdoor air temperature in the "city" dataset.
            temp_std (float): Standard deviation of outdoor air temperature in the "city" dataset.
    '''

    if city == 'atlanta':
        power_mean = 1316421.4548172045
        power_std = 315130.72659125103
        temp_mean = 25.44871527777778
        temp_std = 3.234439481087087
    elif city == 'denver':
        power_mean = 1272671.4671975863
        power_std = 357222.34846292157
        temp_mean = 21.418472222222224
        temp_std = 5.880239074106885

    return power_mean, power_std, temp_mean, temp_std


def predict(model,
            device,
            u,
            power_mean,
            power_std,
            temp_mean,
            temp_std):

    '''
    Run model and format predictions for output.

        Parameters:
            model (pl.LightningModule): An instance of the DNN from checkpoint.
            device (torch.device): Either GPU or CPU.
            u (dict): {
                outT_arr (float, length: n): [...], Outdoor temperature for the future n steps,
                Tset_arr (float, length: [1,n]): [...], Current zone cooling temperature set points for choice length of timesteps

                T_arr (float, length: [1,n]): [...], Current zone temperature set points for choice length of timesteps

                outT (float, scalar): Current outdoor temperature (currently not being used)
                half_hour (float, scalar): Index of half hour through a week. Range = [0,240)
                }
            power_mean (float): Mean power flexibility value used for standardization
            power_std (float): Standard deviation power flexibility value used for standardization
            temp_mean (float): Mean outdoor air temp value used for standardization
            temp_std (float): Standard deviation outdoor air temp value used for standardization

        Returns:
            y (dict): {
                base_power_arr (float, length:n): [...], Baseline power (zone temperature setpoints are kept unchanged) for the future n steps,
                power_upper_arr (float, length:n): [...], Power flexibility upper bounds for the future n steps,
                power_lower_arr (float, length: n): [...], Power flexibility lower bounds for the future n steps,
                }
    '''
    y = {'base_power_arr': [],
        'power_upper_arr': [],
        'power_lower_arr': []}

    # Scalar for outdoor air temp at the time step of interest, normalized for model comprehension
    temp = (torch.tensor(u['outT_arr']) - temp_mean) / temp_std # We are only making one prediction at the first index of this array

    # Whole building set point (assuming single zone for now)
    set_point = interpret_set_point(torch.tensor(u['Tset_arr']), len(temp))

    output_list = ['power_upper_arr', 'base_power_arr', 'power_lower_arr']

    # Disable pytorch gradient calculation
    # Get model predictions and convert back to watt space
    with torch.no_grad():
        for k in range(len(temp)): # For each half hour in the future
            half_hour = int(u['half_hour']) + k

            x = create_input(set_point[k], temp[k], pos_encodings(half_hour)).to(device)

            pred = model(x)
            for j, prediction in enumerate(output_list): # Let j be the index in model output tuple
                # Appending to python list to support future time step predictions
                y[prediction].append(pred[j].cpu().item()) # Take tensor off of GPU and convert to numpy array

        #Turn them all into numpy arrays and convert back to watt space
        for prediction in output_list:
            y[prediction] = np.asarray(y[prediction])
            y[prediction] *= power_std
            y[prediction] += power_mean

    # For the power flexibility, subtract base power to get flexibility
    y['power_lower_arr'] -= y['base_power_arr']
    y['power_upper_arr'] -= y['base_power_arr']
    for key, value in y.items():
        y[key] = value.tolist()
    return json.dumps(y)

def load_model(checkpoint, config):
    '''
    Loads pytorch lightning module from checkpoint.

        Parameters:
            checkpoint (str): Path to pytorch model checkpoint.
            config (dict): Dictionary of model hyperparameters.

        Returns:
            model (pl.LightningModule): An instance of the DNN from checkpoint.
            device (torch.device): Either GPU or CPU.
    '''

    # Load model on local GPU, if possible
    if torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')

    model = DNN.load_from_checkpoint(checkpoint, **config)
    model.to(device)
    model.eval()

    return model, device

def interpret_set_point(T_arr, n):
    '''
    Creates pytorch tensor of the correct shape for correct input shape

        Parameters:
            T_arr (pytorch tensor): T_arr as is from "u" dictionary
            n (int): # of time steps the model is set to predict

        Returns:
            T_arr (pytorch tensor): A modified 1D tensor with appropriate set points for each time step
    '''

    if len(T_arr.shape) > 1:
        raise Exception('Until we start predicting per zone, T_arr should be a 1D vector')
    if len(T_arr) < 1 or len(T_arr) > n:
        raise Exception('T_arr needs to describe a single set point across the prediction horizon \
        or any amount < n of set points for each time step')

    # If len(T_arr) < n, all time steps from index len(T_arr)-1 to n-1 should be equal to the last known set point.
    #   That is, the final set point in the input vector T_arr.
    temp = torch.ones((n - len(T_arr)))
    temp *= T_arr[-1]
    temp = torch.cat((T_arr, temp)) # Assume all set points after index len(T_arr) - 1 are broadcast out to n time steps
    return temp

def pos_encodings(half_hour):
    '''
    Generates positional encodings for all 240 half hours, returns the half hour of interest.

        Parameters:
            half_hour (float): Half hour index. Integer in the range [1,240)

        Returns:
            hour_sin_encoding[half_hour] (pytorch float32): A sine unit circle representation of hour in day.
            hour_cos_encoding[half_hour] (pytorch float32): A cosine unit circle representation of hour in day.
            day_sin_encoding[half_hour] (pytorch float32): A sine unit circle representation of day of week.
            day_cos_encoding[half_hour] (pytorch float32): A cosine unit circle representation of day of week.
    '''

    general_indexed_array = torch.arange(240)
    hour_sin_encoding = torch.sin(general_indexed_array/48 * (2*math.pi))
    hour_cos_encoding = torch.cos((general_indexed_array/48 * (2*math.pi)))
    day_sin_encoding = torch.sin((general_indexed_array//48) / 7 * (2*math.pi))
    day_cos_encoding = torch.cos((general_indexed_array//48) / 7 * (2*math.pi))

    return hour_sin_encoding[half_hour],hour_cos_encoding[half_hour],day_sin_encoding[half_hour],day_cos_encoding[half_hour]

def create_input(sp, temp, positions):
    '''
    Create input tensor for pytorch model.

        Parameters:
            sp (float): Setpoint(s) of previous time step for m zones.
            temp (float): Outdoor air temperature of current time step.
            positions (pytorch tensor): Positional encodings for time of day and day of week.

        Returns:
            x (pytorch tensor): Concatenated input features.
    '''
    x = torch.cat((
        sp.unsqueeze(0).clone().detach(), # (1,1)
        temp.unsqueeze(0).clone().detach(), # (1,1)
        torch.tensor(positions) # (4,1)
        )).type(torch.float32)
    return x
