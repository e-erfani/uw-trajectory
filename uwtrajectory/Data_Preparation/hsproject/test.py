# -*- coding: utf-8 -*-
"""
test.py: unittesting for other files in this folder.
"""
from __future__ import division, unicode_literals
import unittest
import numpy as np
import sys
from util import find_nearest_index, find_index_limits, get_region
from util import get_trajectory_data, trajectory_stats
from util import interpolate_to_regular_grid, closest_index


class FindNearestIndexTests(unittest.TestCase):

    def setUp(self):
        self.integers = range(10)
        self.floats = [float(i) for i in range(10)]
        self.np_floats = np.array(self.floats, dtype=np.float)
        self.single_list = [1., ]

    def tearDown(self):
        self.integers, self.floats, self.np_floats = None, None, None

    def test_returns_integer(self):
        """
        Make sure the function returns an integer.
        """
        self.assertIsInstance(find_nearest_index(5, self.integers), int,
                              'The function should return an integer')

    def test_int_identical(self):
        """
        Test case with integers where the number being searched for is in
        the iterable.
        """
        self.assertEqual(find_nearest_index(5, self.integers), 5,
                         'The nearest index is 5 because 5 is our search' +
                         'value and there is an equal value in the data')
        self.assertEqual(find_nearest_index(0, self.integers), 0,
                         'The nearest index is 0 because 0 is our search' +
                         'value and there is an equal value in the data')
        self.assertEqual(find_nearest_index(9, self.integers), 9,
                         'The nearest index is 9 because 9 is our search' +
                         'value and there is an equal value in the data')

    def test_float_identical(self):
        """
        Test case with floats where the number being searched for is in
        the iterable.
        """
        self.assertEqual(find_nearest_index(5., self.floats), 5.,
                         'The nearest index is 5 because 5 is our search' +
                         'value and there is an equal value in the data')

    def test_float_identical_from_int(self):
        """
        Test case with an integer search value where the float that is equal
        to the integer is in the iterable.
        """
        self.assertEqual(find_nearest_index(5, self.floats), 5.,
                         'Nearest index is 5 since our search term is 5 and' +
                         'there is a float in the data with an equal value')

    def test_int_above(self):
        """
        Test case with a float search value where the float is greater
        than the closest int in the iterable consisting of ints
        """
        self.assertEqual(find_nearest_index(4.1, self.integers), 4,
                         'The nearest index to the float 4.1 is the integer 4')

    def test_int_below(self):
        """
        Test case with a float search value where the float is less
        than the closest int in the iterable consisting of ints
        """
        self.assertEqual(find_nearest_index(5.6, self.integers), 6,
                         'The nearest index to the float 5.6 is the integer 6')

    def test_float_above(self):
        """
        Test case with a float search value where the float is greater
        than the closest float in the iterable consisting of floats
        """
        self.assertEqual(find_nearest_index(4.1, self.floats), 4,
                         'The nearest index to the float 4.1 is the float 4')

    def test_float_below(self):
        """
        Test case with a float search value where the float is less
        than the closest float in the iterable consisting of floats
        """
        self.assertEqual(find_nearest_index(5.6, self.floats), 6,
                         'The nearest index to the float 5.6 is the float 6')

    def test_float_out_of_range(self):
        """
        Test case with a float search value where the float is out of the range
        of the iterable consisting of floats
        """
        self.assertEqual(find_nearest_index(-10., self.floats), 0.,
                         'Nearest index is 0 because our search value is' +
                         'below the minimum index and 0 is the lowest index' +
                         'value')

    def test_single_above(self):
        """
        Test case with a single number in the iterable with a search value that
        is higher than the number in the iterable
        """
        self.assertEqual(find_nearest_index(5., self.single_list), 0,
                         'The nearest index must be 0 because 0 is the only' +
                         'index)')

    def test_single_below(self):
        """
        Test case with a single number in the iterable with search value that
        is lower than the number in the iterable
        """
        self.assertEqual(find_nearest_index(-5., self.single_list), 0,
                         'The nearest index must be 0 because 0 is the only' +
                         'index)')

    def test_single_equal(self):
        """
        Test case with a single number in the iterable with search value
        equal to the number in the iterable
        """
        self.assertEqual(find_nearest_index(1., self.single_list), 0,
                         'The nearest index must be 0 because 0 is the only' +
                         'index)')


class FindIndexLimitTests(unittest.TestCase):

    def setUp(self):
        self.axis = np.array([1., 2, 3, 4, 5, 6, 7, 8, 9, 10])

    def tearDown(self):
        self.axis = None

    def test_is_min_an_int(self):
        """
        Tests that the min_index is an integer because otherwise it would be
        impossible for it to be an index
        """
        min_index, max_index = find_index_limits(self.axis, [self.axis[2],
                                                 self.axis[5]])
        self.assertIsInstance(min_index, int,
                              'The min index should be an integer')

    def test_is_max_an_int(self):
        """
        Tests that the max_index is an integer because otherwise it would be
        impossible for it to be an index
        """
        # see if this is necessary or can this assert statement be put in the
        # previous method
        min_index, max_index = find_index_limits(self.axis, [self.axis[2],
                                                 self.axis[5]])
        self.assertTrue(type(max_index) is int,
                        'The max index should be an integer')

    def test_if_inclusive(self):
        """
        Tests if the method is inclusive by setting the limits equal to an
        item in the array
        """
        min_index, max_index = find_index_limits(self.axis, [self.axis[2],
                                                 self.axis[5]])
        self.assertTrue(min_index == 2 and max_index == 5,
                        'Min/max must be 2 and 5 since these are our set' +
                        'limits and there are indices at these values.')

    def test_single_index(self):
        min_index, max_index = find_index_limits(self.axis, (3, 3))
        self.assertTrue(min_index == 2 and max_index == 2,
                        'Min and max must both be 2 because 2 is the only' +
                        'index between the set limits')

    def test_in_between_numbers(self):
        min_index, max_index = find_index_limits(self.axis, (1.5, 5.5))
        self.assertTrue(min_index == 1 and max_index == 4)

    def test_single_index_in_between_numbers(self):
        min_index, max_index = find_index_limits(self.axis, (3.5, 4.5))
        assert(min_index == 3 and max_index == 3)


class ClosestIndexTests(unittest.TestCase):
    def setUp(self):
        self.lat_traj = 53
        self.lon_traj = 52
        self.regular_lat = np.array([[50, 51, 52, 53, 54],
                                     [50, 51, 52, 53, 54],
                                     [50, 51, 52, 53, 54],
                                     [50, 51, 52, 53, 54],
                                     [50, 51, 52, 53, 54]])
        self.regular_lon = np.array([[50, 50, 50, 50, 50],
                                     [51, 51, 51, 51, 51],
                                     [52, 52, 52, 52, 52],
                                     [53, 53, 53, 53, 53],
                                     [54, 54, 54, 54, 54]])
        self.offset_lat = np.array([[50, 51, 52, 53, 54],
                                    [51, 52, 53, 54, 55],
                                    [50, 51, 52, 53, 54],
                                    [50, 51, 52, 53, 54],
                                    [50, 51, 52, 53, 54]])
        self.offset_lon = np.array([[50, 51, 52, 53, 54],
                                    [51, 52, 53, 54, 55],
                                    [52, 53, 54, 55, 56],
                                    [53, 54, 55, 56, 57],
                                    [54, 55, 56, 57, 58]])

    def tearDown(self):
        self.offset_lon, self.offset_lat, self.regular_lon = None, None, None
        self.regular_lat, self.lon_traj, self.lat_traj = None, None, None

    def test_regular_index_exact_value(self):
        index = closest_index(self.lat_traj, self.lon_traj, self.regular_lat,
                              self.regular_lon)
        self.assertTrue((index == np.array([2, 3])).all(), 'the function ' +
                        'should return the index that contains the two ' +
                        'values from lat and lon_traj')


class GetRegionTests(unittest.TestCase):

    def setUp(self):
        self.latsInt = np.array([1, 2, 3, 4, 5])
        self.longsInt = np.array([1, 3, 5, 7, 9])
        self.D2Ints = np.array([[1, 2, 3, 4, 5],
                                [2, 3, 4, 5, 6],
                                [3, 4, 5, 6, 7],
                                [4, 5, 6, 7, 8],
                                [5, 6, 7, 8, 9]])

    def tearDown(self):
        self.latsInt, self.longsInt, self.D2Ints = None, None, None

    def test_returns_correct_type(self):
        """Test that get_region returns a numpy array"""
        array, lat, lon = get_region(self.D2Ints, self.latsInt,
                                     self.longsInt, (0, 3), (0, 3))
        assert(type(array) is np.ndarray and type(lat) is np.ndarray and
               type(lon) is np.ndarray)

    def test_returns_correct_dtype(self):
        """
        Test that giving an int64 array as input returns an int64 array
        as output.
        """
        int64_array = np.zeros((4,), dtype=np.int64)
        array, lat, lon = get_region(self.D2Ints, int64_array, self.longsInt,
                                     (0, 0), (0, 10))
        assert lat.dtype == np.int64

    def test_region_outside_bounds(self):
        """
        Test that asking for a region outside the bounds of the input array
        returns an empty 2D array with the same dtype as the input array.
        """
        array, lat, lon = get_region(self.D2Ints, self.latsInt, self.longsInt,
                                     (11, 15), (-5, -2))
        assert len(array) == 0 and len(lat) == 0 and len(lon) == 0

    def test_only_max(self):
        new_array, lat, lon = get_region(self.D2Ints, self.latsInt,
                                         self.longsInt, [0, 3], [0, 3])
        assert((new_array == np.array([[1, 2], [2, 3], [3, 4]])).all() and
               (lat == np.array([1, 2, 3])).all() and
               (lon == np.array([1, 3]))).all()


class GetTrajectoryDataTests(unittest.TestCase):
    def setUp(self):
        self.array = np.array([[2, 3, 4, 3, 2],
                               [1, 3, 4, 5, 2],
                               [4, 2, 1, 3, 2],
                               [3, 6, 4, 2, 3],
                               [2, 4, 3, 6, 7]])
        self.lat = np.array([[31, 32, 33, 34, 35],
                             [31, 32, 33, 34, 35],
                             [31, 32, 33, 34, 35],
                             [31, 32, 33, 34, 35],
                             [31, 32, 33, 34, 35]])
        self.lon = np.array([[12, 12, 12, 12, 12],
                             [13, 13, 13, 13, 13],
                             [14, 14, 14, 14, 14],
                             [15, 15, 15, 15, 15],
                             [16, 16, 16, 16, 16]])
        self.lat_traj = np.array([31, 32, 33, 34])
        self.lon_traj = np.array([12, 15, 14, 16])
        self.dist_max = 0

    def tearDown(self):
        self.array, self.lon, self.lat = None, None, None
        self.lat_traj, self.lon_traj = None, None

    def test_correct_length_of_trajectory(self):
        """
        Tests if the second index which is the amount of grids taken is equal
        to the amount of elements in lat_traj and lon_traj
        """
        test_array = get_trajectory_data(self.array, self.lat, self.lon,
                                         self.lat_traj, self.lon_traj, 0)
        assert(len(test_array) == len(self.lat_traj))

    def test_non_equal_trajectory(self):
        """
        tests if a trajectory that is not equal to the lat and lon values will
        still return correct input by having one value be higher than the point
        it is suppost to represent and one value lower
        """
        test_array = get_trajectory_data(self.array, self.lat, self.lon,
                                         np.array([32.7]), np.array([13.3]), 0)
        assert(test_array[0] == 4)

    def test_point(self):
        """
        Tests to see if the correct set of points are picked if there is
        a window of 1 value.
        """
        test_array = get_trajectory_data(self.array, self.lat, self.lon,
                                         self.lat_traj, self.lon_traj, 0)
        # use if testing with delta_index
        assert((test_array == np.array([[2], [6], [1], [6]])).all())
        # use if testing with a distance rather than a delta_index
        # assert((test_array == np.array([[2], [6], [1], [6]])).all())

    def test_nonzero_dist_max(self):
        """
        Tests the output of a grid that has more data than 1 value
        """
        # for lat index: 1, 3, 2
        # for lon index: 1, 2, 3
        test_array = get_trajectory_data(self.array, self.lat, self.lon,
                                         np.array([32.3, 34.1, 33]),
                                         np.array([12.7, 14.2, 14.8]), 1)
        # use if testing with delta_index
        # use if testing with a distance rather than a delta_index
        # assert((test_array == np.array([[2], [6], [1], [6]])).all())
        # assert((test_array == np.array([[3, 3, 4, 4],
        #                                  [3, 2, 2],
        #                                  [4, 1]])).any())


class TrajectoryStatsTests(unittest.TestCase):
    def setUp(self):
        self.array = np.array([[2, 3, 4, 5],
                               [4, 6, 3, 2],
                               [3, 7, 5, 4],
                               [4, 8, 9, 1]])
        self.lat = np.array([34, 34.5, 35, 35.5])
        self.lon = np.array([140, 140.5, 141, 141.5])
        self.lat_traj = np.array([34.6, 34, 35.4, 35])
        self.lon_traj = np.array([140.3, 140, 141.6, 141.1])

    def tearDown(self):
        self.array, self.lat, self.lon, self.lat_traj = None, None, None, None
        self.lon_traj = None

    def test_single_point(self):
        """
        tests if the standard deviation is equal to 0 when there is only one
        point and if the mean is equal to the single value obtained
        """
        mean, std = trajectory_stats(self.array, self.lat, self.lon,
                                     self.lat_traj, self.lon_traj,
                                     delta_index=0)
        assert((std == np.array([0, 0, 0, 0])).all() and
               (mean == np.array([6, 2, 1, 5])).all())

    def test_more_than_single_point(self):
        """
        tests if the standard deviation and mean is equal to what is expected
        for a grid that has a side_length greater than one.
        """
        mean, std = trajectory_stats(self.array, self.lat, self.lon,
                                     np.array([34.5, 35, 34.5, 35]),
                                     np.array([140.5, 140.5, 141, 141]),
                                     delta_index=1)
        assert((mean == np.array([37/9, 49/9, 39/9, 45/9])).all() and
               (.0001 > np.abs(std[0] - 1.52348)))


class InterpolateToRegularGridTests(unittest.TestCase):
    def setUp(self):
        self.array = np.array([[1, 2, 3, 4, 5],
                               [2, 3, 4, 5, 6],
                               [6, 4, 3, 2, 7],
                               [5, 3, 6, 4, 7],
                               [5, 3, 6, 4, 7]])
        self.lat_in = np.array([[1, 2, 3, 4, 5],
                                [1, 2, 3, 4, 5],
                                [2, 3, 4, 5, 6],
                                [1, 2, 3, 4, 5],
                                [1, 2, 3, 4, 5]])
        self.lon_in = np.array([[1, 1, 1, 1, 1],
                                [2, 2, 2, 2, 2],
                                [3, 3, 3, 3, 3],
                                [4, 4, 4, 4, 4],
                                [5, 5, 5, 5, 5]])
        self.lat_out = np.array([3, 3.5, 4])
        self.lon_out = np.array([2, 2.5, 3])

    def tearDown(self):
        self.array, self.lat_in, self.lon_in = None, None, None
        self.lat_out, self.lat_in = None, None

    def test_with_regular_grid_with_slight_offset(self):
        """
        tests the output with a regular grid with one set of lat values offset
        by 1.
        """
        easy_array = np.array([[0, 0, 0, 0, 0],
                               [1, 1, 1, 1, 1],
                               [2, 2, 2, 2, 2],
                               [3, 3, 3, 3, 3],
                               [4, 4, 4, 4, 4]])
        function = interpolate_to_regular_grid(easy_array, self.lat_in,
                                               self.lon_in,
                                               self.lat_out, self.lon_out)
        self.assertAlmostEqual(function[1][0], 1.5)
        # if that test is too restrictive then use this one:
        # assert(np.abs(1.5-function[0][1]) < .1)

    def test_grid_with_exact_values(self):
        """
        test with lat_out and lon_out representing values that are equal to
        values in lat_in and lon_in
        """
        exact_lat_out = np.array([3, 4])
        exact_lon_out = np.array([2, 3, 4])
        function = interpolate_to_regular_grid(self.array, self.lat_in,
                                               self.lon_in,
                                               exact_lat_out, exact_lon_out)
        self.assertAlmostEqual(function[2][1], 4)


if __name__ == '__main__':
    # This command runs all tests in the module
    # unittest.main()

    # These commands let you select specific tests to run

    # initialize a test suite and loader
    suite = unittest.TestSuite()
    loader = unittest.TestLoader()

    # add test cases
    # test cases for find_nearest_index
    # suite.addTest(loader.loadTestsFromTestCase(FindNearestIndexTests))
    # test cases for find_index_limits
    # suite.addTest(loader.loadTestsFromTestCase(FindIndexLimitTests))
    # test cases for closest_index
    # suite.addTest(loader.loadTestsFromTestCase(ClosestIndexTests))
    # test cases for get_region
    # suite.addTest(loader.loadTestsFromTestCase(GetRegionTests))
    # test cases for get_trajectory_data
    suite.addTest(loader.loadTestsFromTestCase(GetTrajectoryDataTests))
    # test cases for trajectory_stats
    # suite.addTest(loader.loadTestsFromTestCase(TrajectoryStatsTests))
    # test cases for interpolate_to_regular_grid
    # suite.addTest(loader.loadTestsFromTestCase(InterpolateToRegularGridTests))

    # initialize a test runner and run the test suite
    runner = unittest.TextTestRunner(stream=sys.stdout, verbosity=2)
    result = runner.run(suite)
