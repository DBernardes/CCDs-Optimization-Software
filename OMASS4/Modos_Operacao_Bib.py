#!/usr/bin/env python
# coding: utf-8
# 24/10/2019. Denis Varise Bernardes.


class ModosOperacao:
    def __init__(self):
        self.em_mode = []
        self.em_gain = []
        self.hss = []
        self.preamp = []
        self.binn = []
        self.sub_img = []
        self.t_exp = []

        self.modos_operacao = []
        self.modo_atual = {}

    def write_mode(
        self, em_mode, em_gain, hss, preamp, binn, sub_img, max_t_exp, min_t_exp=0.00001
    ):
        # Write the operation mode to the class in the dictionaru format
        dic = {
            "em_mode": em_mode,
            "em_gain": em_gain,
            "hss": hss,
            "preamp": preamp,
            "binn": binn,
            "sub_img": sub_img,
            "max_t_exp": max_t_exp,
            "min_t_exp": min_t_exp,
        }
        self.modos_operacao.append(dic)

    def write_list_of_modes(self, lista):
        # Write an entire list of modes to the class
        self.modos_operacao = lista

    def get_list_of_modes(self):
        # returns the list of modes
        return self.modos_operacao

    def clear_list_of_modes(self):
        # clear the list of modes
        self.modos_operacao = []
