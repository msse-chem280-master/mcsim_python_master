"""
Tests for Python Standard Library implementation of MC code.
"""

import math
import random

import pytest

from monte_carlo import calculate_distance, accept_or_reject

def test_calculate_distance():

    point1 = [0, 0, 0]
    point2 = [0, 0, 1]

    expected_value = 1
    observed_value = calculate_distance(point1, point2)
    
    assert expected_value == observed_value

def test_calculate_distance2():
    point_3 = [0, 0, 0]
    point_4 = [0, 1, 1]

    dist2 = calculate_distance(point_3, point_4)
    assert math.sqrt(2) == dist2

def test_calculate_distance_periodic():
    point_1 = [0, 0, 0]
    point_2 = [0, 0, 8]
    box_length = 10

    expected_distance = 2
    observed_distance = calculate_distance(point_1, point_2, box_length=box_length)

    assert expected_distance == observed_distance

def test_accept_or_reject_positive_false():
    try:
        random.seed(0)
        delta_e = 1
        beta = 1
        assert accept_or_reject(delta_e, beta) is False
    finally:
        random.seed()

"""
@pytest.mark.parametrize('variable_name1, variable_name2, ... , variable_nameN', [
    (variable_value1, variable_value2, ... , variable_valueN),
    (variable_value1, variable_value2, ... , variable_valueN),
]
)
def test_function(variable_name1, variable_name2, ... , variable_nameN):
    *** TEST CODE HERE ***
"""

@pytest.mark.parametrize("point1, point2, expected_distance, box_length", [
    ([0, 0, 0], [1, 0, 0], 1, None ),
    ([0, 0, 0], [8, 0, 0], 8, None ),
    ([0, 0, 0], [8, 0, 0], 2, 10 ),
])
def test_calculate_distance_many(point1, point2, expected_distance, box_length):

    observed_distance = calculate_distance(point1, point2, box_length)

    assert expected_distance == observed_distance
