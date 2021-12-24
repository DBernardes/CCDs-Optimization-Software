#!/usr/bin/env python
# coding: utf-8
# 24/10/2019. Denis Varise Bernardes.


class Operation_Modes:
    """Operation Modes class.

    This class store a list of python dictionaries, where each dic has
    a CCD operation mode.
    """

    def __init__(self):

        self.operation_modes = []

    def write_operation_mode(self, dic):
        """Append the CCD operation mode intto the class."""
        self.operation_modes.append(dic)

    def write_list_operation_modes(self, list):
        """Write an entire list of modes into the class."""
        self.operation_modes = list

    def get_list_operation_modes(self):
        """Returns the list of modes."""
        return self.operation_modes

    def clear_list_operation_modes(self):
        """Clear the list of operation modes."""
        self.operation_modes = []
