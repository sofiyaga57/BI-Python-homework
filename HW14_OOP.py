import numpy as np
import re
from abc import ABC, abstractmethod


class User:
    def __init__(self, name: str):
        self.name = name


class Booking:
    def __init__(self, equipment: str, user: str, start_time, end_time):
        self.equipment = equipment
        self.user = user
        self.start_time = start_time
        self.end_time = end_time

    def is_intersect(self, other) -> bool:
        if self.end_time <= other.start_time or self.start_time >= other.end_time:
            return False
        return True


class LabEquipment:
    def __init__(self, equipment):
        self.equipment = equipment
        self.booking_history = []

    def is_available(self, equipment_name: str, user: str, start_time, end_time) -> bool:
        for booking in self.booking_history:
            if booking.equipment == equipment_name:
                if (booking.user != user and
                        booking.is_intersect(Booking(equipment_name, user, start_time, end_time))):
                    return False
        return True

    def book(self, equipment_name: str, user: str, start_time, end_time):
        if equipment_name not in self.equipment:
            raise ValueError(f"{equipment_name} is not available in the lab.")

        new_booking = Booking(equipment_name, user, start_time, end_time)
        if self.is_available(equipment_name, user, start_time, end_time):
            self.booking_history.append(new_booking)
            self.equipment.append(equipment_name)
            print("Successful booking!")
        else:
            raise ValueError("Booking not available")


class WrongAlphabetError(ValueError):
    pass


class WrongInputError(ValueError):
    pass


class GenCodeInterpreter:
    def __init__(self):
        self.memory = np.zeros(500)
        self.pointer = 0
        self.buffer = ''

    def eval(self, string: str):
        pattern = r'^[GCTAN]+$'
        if not bool(re.match(pattern, string)):
            raise WrongAlphabetError('You can only use ATGC-alphabet + N character')
        if not string.endswith('N'):
            raise WrongInputError('The last character should be N')
        self.buffer = ''
        for char in list(string):
            match char:
                case 'A':
                    self.pointer += 1
                case 'T':
                    self.pointer -= 1
                case 'G':
                    self.memory[self.pointer] += 1
                case 'C':
                    self.memory[self.pointer] -= 1
                case 'N':
                    self.buffer += chr(int(self.memory[self.pointer]))
        return self.buffer


def meet_the_dunders():
    res = 0

    matrix = []
    for idx in range(0, 100, 10):
        matrix.__iadd__([list(range(idx, idx + 10))])

    def func_1(x):
        return x.__eq__(1) or x.__eq__(3)

    def func_2(x):
        return [x.__getitem__(col) for col in selected_columns_indices]

    selected_columns_indices = list(filter(func_1.__call__, range(len(matrix))))
    selected_columns = map(func_2.__call__, matrix)

    arr = np.array(list(selected_columns))

    mask = arr.__getitem__((slice(None), 1)).__mod__(3).__eq__(0)
    new_arr = arr.__getitem__(mask)

    product = new_arr.__matmul__(new_arr.__getattribute__('T'))

    if (product.__getitem__(0).__lt__(1000)).all() and (product.__getitem__(2).__gt__(1000)).any():
        res = int(product.__getattribute__('mean')().__floordiv__(10).__mod__(100))
    return res


class BiologicalSequence(ABC):
    @abstractmethod
    def check_alphabet(self, sequence):
        pass

    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, position):
        pass

    @abstractmethod
    def __str__(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    valid_nucleotides = "ACGTU"

    def __init__(self, sequence):
        if not self.check_alphabet(sequence):
            raise ValueError("Invalid DNA nucleotide alphabet.")
        self.sequence = sequence.upper()

    def check_alphabet(self, sequence):
        return set(sequence.upper()).issubset(set(self.valid_nucleotides))

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, position):
        return self.sequence[position]

    def __str__(self):
        return self.sequence

    def complement(self):
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        complement_sequence = ''.join([complement_dict[nuc] for nuc in self.sequence])
        return NucleicAcidSequence(complement_sequence)

    def gc_content(self):
        return ((self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence)) * 100
