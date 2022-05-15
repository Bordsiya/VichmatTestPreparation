from dataclasses import dataclass


@dataclass
class TableFunction:
    n: int
    x_arr: list[float]
    y_arr: list[float]

