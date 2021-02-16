#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 21:25:37 2021
@author: rga
@contact: roy_gonzalez@fq.uh.cu, roy.gonzalez-aleman@u-psud.fr
"""
import pytest
import pickle


def unpickle_from_file(file_name):
    ''' Unserialize a **pickle** file.

    Args:
        file_name (str): file to unserialize.
    Returns:
        (object): an unserialized object.
    '''
    with open(file_name, 'rb') as file:
        data = pickle.load(file)
    return data


golden = 'aligned_tau_5ef8f444.pick'
target = '/home/rga/BSProject/runners/BitQT/bitQT_outputs/aligned_tau.pick'
final_array = unpickle_from_file(golden)
tinal_array = unpickle_from_file(target)


def test_final_clustering():
    assert (final_array == tinal_array).sum() == final_array.size
