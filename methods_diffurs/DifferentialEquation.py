from collections import Callable
from dataclasses import dataclass


@dataclass
class DifferentialEquation:
    f: Callable[[float, float], float]
    x_0: float
    y_0: float
    a: float
    b: float
    h: float
    e: float
